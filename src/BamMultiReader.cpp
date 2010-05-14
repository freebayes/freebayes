// ***************************************************************************
// BamMultiReader.cpp (c) 2010 Erik Garrison
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ---------------------------------------------------------------------------
// Last modified: 23 Februrary 2010 (EG)
// ---------------------------------------------------------------------------
// Uses BGZF routines were adapted from the bgzf.c code developed at the Broad
// Institute.
// ---------------------------------------------------------------------------
// Functionality for simultaneously reading multiple BAM files.
//
// This functionality allows applications to work on very large sets of files
// without requiring intermediate merge, sort, and index steps for each file
// subset.  It also improves the performance of our merge system as it
// precludes the need to sort merged files.
// ***************************************************************************

// C++ includes
#include <algorithm>
#include <iterator>
#include <string>
#include <vector>
#include <iostream>
#include <sstream>

// BamTools includes
#include "BGZF.h"
#include "BamMultiReader.h"
using namespace BamTools;
using namespace std;

// -----------------------------------------------------
// BamMultiReader implementation
// -----------------------------------------------------

// constructor
BamMultiReader::BamMultiReader(void)
    : CurrentRefID(0)
    , CurrentLeft(0)
{ }

// destructor
BamMultiReader::~BamMultiReader(void) {
    Close(); // close the bam files
    // clean up reader objects
    for (vector<BamReader*>::iterator br = readers.begin(); br != readers.end(); ++br) {
        delete *br;
    }
}

// close the BAM files
void BamMultiReader::Close(void) {
    int index = 0;
    for (vector<BamReader*>::iterator br = readers.begin(); br != readers.end(); ++br) {
        BamReader* reader = *br;
        reader->Close();  // close the reader
        readerStates.at(index++) = CLOSED;  // track its state
    }
}

// updates the reference id stored in the BamMultiReader
// to reflect the current state of the readers
void BamMultiReader::UpdateReferenceID(void) {
    bool updateRefID = true;
    int i = 0;
    for (vector<BamAlignment*>::iterator it = alignments.begin(); it != alignments.end(); ++it) {
        BamAlignment* alignment = *it;
        if (readerStates.at(i++) != END && alignment->RefID == CurrentRefID) {
            updateRefID = false;
            break;
        }
    }
    if (updateRefID) {
        // get the next reference id
        // while there aren't any readers at the next ref id
        // increment the ref id
        int nextRefID = CurrentRefID;
        bool ok = false;
        while (!ok) {
            ++nextRefID;
            int i = 0;
            for (vector<BamAlignment*>::iterator it = alignments.begin(); it != alignments.end(); ++it) {
                BamAlignment* alignment = *it;
                if (readerStates.at(i++) != END && alignment->RefID == nextRefID) {
                    ok = true;
                    break;
                }
            }
        }
        cerr << "updating reference id from " << CurrentRefID << " to " << nextRefID << endl;
        CurrentRefID = nextRefID;
    }
}

// checks if any readers still have alignments
bool BamMultiReader::HasOpenReaders() {
    for (vector<BamReaderState>::iterator it = readerStates.begin(); it != readerStates.end(); ++it) {
        BamReaderState readerState = *it;
        if (readerState != END)
            return true;
    }
    return false;
}

// get next alignment among all files (from specified region, if given)
bool BamMultiReader::GetNextAlignment(BamAlignment& nextAlignment) {

    // bail out if we are at EOF in all files, means no more alignments to process
    if (!HasOpenReaders())
        return false;

    // when all alignments have stepped into a new target sequence, update our
    // current reference sequence id
    UpdateReferenceID();

    // find the alignment with the lowest position that is also in the current
    // reference sequence, getnextalignment on its bam reader to step it
    // forward, and return it
    int i = 0, lowestPosition = -1, lowestAlignmentIndex = -1;
    BamAlignment* lowestAlignment = NULL;
    for (vector<BamAlignment*>::const_iterator it = alignments.begin(); it != alignments.end(); ++it) {
        BamAlignment* ba = *it;
        if (readerStates.at(i) != END && 
                ba->RefID == CurrentRefID && 
                (lowestAlignment == NULL || ba->Position < lowestPosition)) {
            lowestPosition = ba->Position;
            lowestAlignmentIndex = i;
            lowestAlignment = ba;
        }
        ++i;
    }

    // now that we have the lowest alignment in the set, save it by copy to our argument
    nextAlignment = BamAlignment(*lowestAlignment);

    // else continue and load the next alignment
    bool r = (readers.at(lowestAlignmentIndex))->GetNextAlignment(*alignments.at(lowestAlignmentIndex));
    if (!r) {
        cerr << "reached end of file " << readers.at(lowestAlignmentIndex)->GetFilename() << endl;
        readerStates.at(lowestAlignmentIndex) = END;  // set flag for end of file
    }

    return true;
}

// jumps to specified region(refID, leftBound) in BAM files, returns success/fail
bool BamMultiReader::Jump(int refID, int position) {

    //if ( References.at(refID).RefHasAlignments && (position <= References.at(refID).RefLength) ) {
    CurrentRefID = refID;
    CurrentLeft  = position;

    bool result = true;
    for (vector<BamReader*>::iterator br = readers.begin(); br != readers.end(); ++br) {
        BamReader* reader = *br;
        result &= reader->Jump(refID, position);
    }
    if (result)
        UpdateAlignments();
    return result;
}

// opens BAM files
void BamMultiReader::Open(const vector<string> filenames) {
    // for filename in filenames
    fileNames = filenames; // save filenames in our multireader
    for (vector<string>::const_iterator it = filenames.begin(); it != filenames.end(); ++it) {
        string filename = *it;
        BamReader* reader = new BamReader;
        reader->Open(filename, filename + ".bai");
        BamAlignment* alignment = new BamAlignment;
        reader->GetNextAlignment(*alignment);
        readers.push_back(reader); // tracks readers
        alignments.push_back(alignment); // tracks current alignments of the readers
        BamReaderState readerState = START;
        readerStates.push_back(readerState);
    }
}

// Runs GetNextAlignment for all BAM readers; used during jumps
void BamMultiReader::UpdateAlignments(void) {
    int i = 0;
    for (vector<BamAlignment*>::iterator it = alignments.begin(); it != alignments.end(); ++it) {
        BamAlignment* ba = *it;
        readers.at(i++)->GetNextAlignment(*ba);
    }
}

void BamMultiReader::PrintFilenames(void) {
    for (vector<BamReader*>::iterator br = readers.begin(); br != readers.end(); ++br) {
        BamReader* reader = *br;
        cout << reader->GetFilename() << endl;
    }
}

// returns BAM file pointers to beginning of alignment data
bool BamMultiReader::Rewind(void) { 
    bool result = true;
    int index = 0;
    for (vector<BamReader*>::iterator br = readers.begin(); br != readers.end(); ++br) {
        BamReader* reader = *br;
        result &= reader->Rewind();
        readerStates.at(index++) = START;
    }
    return result;
}

// saves index data to BAM index files (".bai") where necessary, returns success/fail
bool BamMultiReader::CreateIndexes(void) {
    bool result = true;
    for (vector<BamReader*>::iterator br = readers.begin(); br != readers.end(); ++br) {
        BamReader* reader = *br;
        result &= reader->CreateIndex();
    }
    return result;
}

// makes a virtual, unified header for all the bam files in the multireader
const string BamMultiReader::GetHeaderText(void) const {

    string mergedHeader = "";

    // foreach extraction entry (each BAM file)
    bool isFirstTime = true;
    for (vector<BamReader*>::const_iterator it = readers.begin(); it != readers.end(); ++it) {

        BamReader* reader = *it;

        stringstream header(reader->GetHeaderText());
        vector<string> lines;
        string item;
        while (getline(header, item))
            lines.push_back(item);

        for (vector<string>::const_iterator it = lines.begin(); it != lines.end(); ++it) {

            // get next line from header, skip if empty
            string headerLine = *it;
            if ( headerLine.empty() ) { continue; }

            // if first file, save HD & SQ entries
            if ( isFirstTime ) {
                if ( headerLine.find("@HD") == 0 || headerLine.find("@SQ") == 0) {
                    mergedHeader.append(headerLine.c_str());
                    mergedHeader.append(1, '\n');
                }
            }

            // (for all files) append RG entries
            if ( headerLine.find("@RG") == 0 ) {
                mergedHeader.append(headerLine.c_str() );
                mergedHeader.append(1, '\n');
            }

        }

        // set iteration flag
        isFirstTime = false;
    }

    // return merged header text
    return mergedHeader;
}


// TODO XXX insecure, we should check that we have the same reference data for all files!!!!
const int BamMultiReader::GetReferenceCount(void) const {
    return readers.front()->GetReferenceCount();
}

// returns vector of reference objects
const BamTools::RefVector BamMultiReader::GetReferenceData(void) const {
    return readers.front()->GetReferenceData();
}

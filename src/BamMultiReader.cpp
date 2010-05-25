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
        //cerr << "updating reference id from " << CurrentRefID << " to " << nextRefID << endl;
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
        //cerr << "reached end of file " << readers.at(lowestAlignmentIndex)->GetFilename() << endl;
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
void BamMultiReader::Open(const vector<string> filenames, bool openIndexes) {
    // for filename in filenames
    fileNames = filenames; // save filenames in our multireader
    for (vector<string>::const_iterator it = filenames.begin(); it != filenames.end(); ++it) {
        string filename = *it;
        BamReader* reader = new BamReader;
        if (openIndexes) {
            reader->Open(filename, filename + ".bai");
        } else {
            reader->Open(filename); // for merging, jumping is disallowed
        }
        BamAlignment* alignment = new BamAlignment;
        reader->GetNextAlignment(*alignment);
        readers.push_back(reader); // tracks readers
        alignments.push_back(alignment); // tracks current alignments of the readers
        BamReaderState readerState = START;
        readerStates.push_back(readerState);
    }
    ValidateReaders();
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

// ValidateReaders checks that all the readers point to BAM files representing
// alignments against the same set of reference sequences, and that the
// sequences are identically ordered.  If these checks fail the operation of
// the multireader is undefined, so we force program exit.
void BamMultiReader::ValidateReaders(void) const {
    int firstRefCount = readers.front()->GetReferenceCount();
    BamTools::RefVector firstRefData = readers.front()->GetReferenceData();
    for (vector<BamReader*>::const_iterator it = readers.begin(); it != readers.end(); ++it) {
        BamTools::RefVector currentRefData = (*it)->GetReferenceData();
        BamTools::RefVector::const_iterator f = firstRefData.begin();
        BamTools::RefVector::const_iterator c = currentRefData.begin();
        if ((*it)->GetReferenceCount() != firstRefCount || firstRefData.size() != currentRefData.size()) {
            cerr << "ERROR: mismatched number of references in " << (*it)->GetFilename()
                      << " expected " << firstRefCount 
                      << " reference sequences but only found " << (*it)->GetReferenceCount() << endl;
            exit(1);
        }
        // this will be ok; we just checked above that we have identically-sized sets of references
        // here we simply check if they are all, in fact, equal in content
        while (f != firstRefData.end()) {
            if (f->RefName != c->RefName || f->RefLength != c->RefLength) {
                cerr << "ERROR: mismatched references found in " << (*it)->GetFilename()
                          << " expected: " << endl;
                for (BamTools::RefVector::const_iterator a = firstRefData.begin(); a != firstRefData.end(); ++a)
                    cerr << a->RefName << " " << a->RefLength << endl;
                cerr << "but found: " << endl;
                for (BamTools::RefVector::const_iterator a = currentRefData.begin(); a != currentRefData.end(); ++a)
                    cerr << a->RefName << " " << a->RefLength << endl;
                exit(1);
            }
            ++f; ++c;
        }
    }
}

// NB: The following functions assume that we have identical references for all
// BAM files.  We enforce this by invoking the above validation function
// (ValidateReaders) to verify that our reference data is the same across all
// files on Open, so we will not encounter a situation in which there is a
// mismatch and we are still live.

// returns the number of reference sequences
const int BamMultiReader::GetReferenceCount(void) const {
    return readers.front()->GetReferenceCount();
}

// returns vector of reference objects
const BamTools::RefVector BamMultiReader::GetReferenceData(void) const {
    return readers.front()->GetReferenceData();
}

const int BamMultiReader::GetReferenceID(const string& refName) const { 
    return readers.front()->GetReferenceID(refName);
}

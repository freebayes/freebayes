// ***************************************************************************
// BamMultiReader.h (c) 2010 Erik Garrison
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ---------------------------------------------------------------------------
// Last modified: 22 February 2010 (EG)
// ---------------------------------------------------------------------------
// Functionality for simultaneously reading multiple BAM files
// ***************************************************************************

#ifndef BAMMULTIREADER_H
#define BAMMULTIREADER_H

// C++ includes
#include <string>

using namespace std;

// BamTools includes
#include "BamAux.h"
#include "BamReader.h"

namespace BamTools {

enum BamReaderState { START, END, CLOSED };

class BamMultiReader {

    // constructor / destructor
    public:
        BamMultiReader(void);
        ~BamMultiReader(void);

    // public interface
    public:

        // positioning
        int CurrentRefID;
        int CurrentLeft;

        // ----------------------
        // BAM file operations
        // ----------------------

        // close BAM files
        void Close(void);
        // performs random-access jump to reference, position
        bool Jump(int refID, int position = 0);
        // opens BAM files (and optional BAM index files, if provided)
        //void Open(const vector<std::string&> filenames, const vector<std::string&> indexFilenames);
        void Open(const vector<string> filenames, bool openIndexes = true);
        // returns file pointers to beginning of alignments
        bool Rewind(void);

        // ----------------------
        // access alignment data
        // ----------------------
        // updates the reference id marker to match the lower limit of our readers
        void UpdateReferenceID(void);

        // retrieves next available alignment (returns success/fail) from all files
        bool GetNextAlignment(BamAlignment&);
        // ... should this be private?
        bool HasOpenReaders(void);

        // ----------------------
        // access auxiliary data
        // ----------------------

        // returns unified SAM header text for all files
        const string GetHeaderText(void) const;
        // returns number of reference sequences
        const int GetReferenceCount(void) const;
        // returns vector of reference objects
        const BamTools::RefVector GetReferenceData(void) const;
        // returns reference id (used for BamMultiReader::Jump()) for the given reference name
        const int GetReferenceID(const std::string& refName) const;
        // validates that we have a congruent set of BAM files that are aligned against the same reference sequences
        void ValidateReaders() const;

        // ----------------------
        // BAM index operations
        // ----------------------

        // creates index for BAM files which lack them, saves to files (default = bamFilename + ".bai")
        bool CreateIndexes(void);

        //const int GetReferenceID(const string& refName) const;

        // utility
        void PrintFilenames(void);
        void UpdateAlignments(void);


    // private implementation
    private:
        // TODO perhaps, for legibility, I should use a struct to wrap them all up
        //      But this may actually make things more confusing, as I'm only
        //      operating on them all simultaneously during GetNextAlignment
        //      calls.
        // all these vectors are ordered the same
        // readers.at(N) refers to the same reader as alignments.at(N) and readerStates.at(N)
        vector<BamReader*> readers; // the set of readers which we operate on
        vector<BamAlignment*> alignments; // the equivalent set of alignments we use to step through the files
        vector<BamReaderState> readerStates; // states of the various readers
        // alignment position?
        vector<string> fileNames;
};

} // namespace BamTools

#endif // BAMMULTIREADER_H

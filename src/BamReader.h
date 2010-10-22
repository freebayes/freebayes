// ***************************************************************************
// BamReader.h (c) 2009 Derek Barnett, Michael Str�mberg
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ---------------------------------------------------------------------------
// Last modified: 13 October 2010 (DB)
// ---------------------------------------------------------------------------
// Uses BGZF routines were adapted from the bgzf.c code developed at the Broad
// Institute.
// ---------------------------------------------------------------------------
// Provides the basic functionality for reading BAM files
// ***************************************************************************

#ifndef BAMREADER_H
#define BAMREADER_H

#include <string>
#include "BamAlignment.h"
#include "BamIndex.h"

namespace BamTools {
  
class BamReader {

    // constructor / destructor
    public:
        BamReader(void);
        ~BamReader(void);

    // public interface
    public:

        // ----------------------
        // BAM file operations
        // ----------------------

        // close BAM file
        void Close(void);
        // returns whether reader is open for reading or not
        bool IsOpen(void) const;
        // performs random-access jump using (reference, position) as a left-bound
        bool Jump(int refID, int position = 0);
        // opens BAM file (and optional BAM index file, if provided)
        // @lookForIndex - if no indexFilename provided, look in BAM file's directory for an existing index file
	//   default behavior is to skip index file search if no index filename given
        // @preferStandardIndex - if true, give priority in index file searching to standard BAM index (*.bai)
	//   default behavior is to prefer the BamToolsIndex (*.bti) if both are available
        bool Open(const std::string& filename, 
                  const std::string& indexFilename = "", 
                  const bool lookForIndex = false, 
                  const bool preferStandardIndex = false);
        // returns file pointer to beginning of alignments
        bool Rewind(void);
        // sets a region of interest (with left & right bound reference/position)
        // returns success/failure of seeking to left bound of region
        bool SetRegion(const BamRegion& region);
        bool SetRegion(const int& leftRefID, const int& leftBound, const int& rightRefID, const int& rightBound);

        // ----------------------
        // access alignment data
        // ----------------------

        // retrieves next available alignment (returns success/fail)
        bool GetNextAlignment(BamAlignment& bAlignment);
        // retrieves next available alignment core data (returns success/fail)
        // ** DOES NOT parse any character data (read name, bases, qualities, tag data) **
        // useful for operations requiring ONLY aligner-related information 
	// (refId/position, alignment flags, CIGAR, mapQuality, etc)
        bool GetNextAlignmentCore(BamAlignment& bAlignment);

        // ----------------------
        // access auxiliary data
        // ----------------------

        // returns SAM header text
        const std::string GetHeaderText(void) const;
        // returns number of reference sequences
        int GetReferenceCount(void) const;
        // returns vector of reference objects
        const BamTools::RefVector& GetReferenceData(void) const;
        // returns reference id (used for BamReader::Jump()) for the given reference name
        int GetReferenceID(const std::string& refName) const;
        // returns the name of the file associated with this BamReader
        const std::string GetFilename(void) const;

        // ----------------------
        // BAM index operations
        // ----------------------

        // creates index for BAM file, saves to file
        // default behavior is to create the BAM standard index (".bai")
        // set flag to false to create the BamTools-specific index (".bti")
        bool CreateIndex(bool useStandardIndex = true);
	// returns whether index data is available for reading 
	// (e.g. if true, BamReader should be able to seek to a region)
	bool HasIndex(void) const;
	// change the index caching behavior
	// default BamReader/Index mode is LimitedIndexCaching
	// @mode - can be either FullIndexCaching, LimitedIndexCaching, 
	//   or NoIndexCaching. See BamIndex.h for more details
        void SetIndexCacheMode(const BamIndex::BamIndexCacheMode mode);
	
    // deprecated methods
    public:
	
	// deprecated (but still available): prefer HasIndex() instead
	//
	// Deprecated purely for API semantic clarity - HasIndex() should be clearer 
	// than IsIndexLoaded() in light of the new caching modes that may clear the 
	// index data from memory, but leave the index file open for later random access 
	// seeks.
	//
	// For example, what would (IsIndexLoaded() == true) mean when cacheMode has been 
	// explicitly set to NoIndexCaching? This is confusing at best, misleading about 
	// current memory behavior at worst.
	//
	// returns whether index data is available 
	// (e.g. if true, BamReader should be able to seek to a region)
        bool IsIndexLoaded(void) const;
        
    // private implementation
    private:
        struct BamReaderPrivate;
        BamReaderPrivate* d;
};

} // namespace BamTools

#endif // BAMREADER_H

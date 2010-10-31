// ***************************************************************************
// BamIndex.cpp (c) 2009 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ---------------------------------------------------------------------------
// Last modified: 21 October 2010 (DB)
// ---------------------------------------------------------------------------
// Provides index functionality - both for the default (standardized) BAM 
// index format (.bai) as well as a BamTools-specific (nonstandard) index 
// format (.bti).
// ***************************************************************************

#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <iostream>
#include <map>
#include "BamIndex.h"
#include "BamReader.h"
#include "BGZF.h"
using namespace std;
using namespace BamTools;

// -------------------------------
// BamIndex implementation

// ctor
BamIndex::BamIndex(BamTools::BgzfData* bgzf, BamTools::BamReader* reader)
    : m_BGZF(bgzf)
    , m_reader(reader)
    , m_cacheMode(BamIndex::LimitedIndexCaching)
    , m_indexStream(0)
{ 
    if ( m_reader && m_reader->IsOpen() ) 
        m_references = m_reader->GetReferenceData();
}

// dtor
BamIndex::~BamIndex(void) {
    if ( IsOpen() )
        fclose(m_indexStream);
}

// return true if FILE* is open
bool BamIndex::IsOpen(void) const {
    return ( m_indexStream != 0 );
}

// loads existing data from file into memory
bool BamIndex::Load(const string& filename)  {

    // open index file, abort on error
    if ( !OpenIndexFile(filename, "rb") ) {
        fprintf(stderr, "ERROR: Unable to open the BAM index file %s for reading.\n", filename.c_str());
        return false;
    }

    // check magic number
    if ( !LoadHeader() ) {
        fclose(m_indexStream);
        return false;
    }

    // load reference data (but only keep in memory if full caching requested)
    bool saveInitialLoad = ( m_cacheMode == BamIndex::FullIndexCaching );
    if ( !LoadAllReferences(saveInitialLoad) ) {
        fclose(m_indexStream);
        return false;
    }

    // update index cache based on selected mode
    UpdateCache();

    // return success
    return true;
}

// opens index file for reading/writing, return true if opened OK
bool BamIndex::OpenIndexFile(const string& filename, const string& mode) {
    m_indexStream = fopen(filename.c_str(), mode.c_str());
    return ( m_indexStream != 0 );
}

// rewind index file to beginning of index data, return true if rewound OK
bool BamIndex::Rewind(void) {
    return ( fseek64(m_indexStream, DataBeginOffset(), SEEK_SET) == 0 );
}

// change the index caching behavior
void BamIndex::SetCacheMode(const BamIndexCacheMode mode) {
    if ( mode != m_cacheMode ) {
        m_cacheMode = mode;
        UpdateCache();
    }
}

// updates in-memory cache of index data, depending on current cache mode
void BamIndex::UpdateCache(void) {

    // skip if file not open
    if ( !IsOpen() ) return;

    // reflect requested cache mode behavior
    switch ( m_cacheMode ) {

        case (BamIndex::FullIndexCaching) :
            Rewind();
            LoadAllReferences(true);
            break;

        case (BamIndex::LimitedIndexCaching) :
            if ( HasFullDataCache() )
                KeepOnlyFirstReferenceOffsets();
            else {
                ClearAllData();
                SkipToFirstReference();
                LoadFirstReference(true);
            }
            break;
        case(BamIndex::NoIndexCaching) :
            ClearAllData();
            break;
        default :
            // unreachable
            ;
    }
}

// writes in-memory index data out to file
bool BamIndex::Write(const string& bamFilename) {

    // open index file for writing
    string indexFilename = bamFilename + Extension();
    if ( !OpenIndexFile(indexFilename, "wb") ) {
        fprintf(stderr, "ERROR: Could not open file to save index.\n");
        return false;
    }

    // write index header data
    if ( !WriteHeader() ) {
        fprintf(stderr, "ERROR: There was a problem writing index metadata to new index file.\n");
        fflush(m_indexStream);
        fclose(m_indexStream);
        exit(1);
    }

    // write main index data
    if ( !WriteAllReferences() ) {
        fprintf(stderr, "ERROR: There was a problem writing index data to new index file.\n");
        fflush(m_indexStream);
        fclose(m_indexStream);
        exit(1);
    }

    // flush any remaining output, rewind file, and return success
    fflush(m_indexStream);
    fclose(m_indexStream);

    // re-open index file for later reading
    if ( !OpenIndexFile(indexFilename, "rb") ) {
        fprintf(stderr, "ERROR: Could not open newly created index file for reading.\n");
        return false;
    }

    // return success/failure of write
    return true;
}

// #########################################################################################
// #########################################################################################

// -------------------------------
// BamStandardIndex structs & typedefs 
 
namespace BamTools { 

// BAM index constants
const int MAX_BIN        = 37450;    // =(8^6-1)/7+1
const int BAM_LIDX_SHIFT = 14;
  
// --------------------------------------------------
// BamStandardIndex data structures & typedefs
struct Chunk {

    // data members
    uint64_t Start;
    uint64_t Stop;

    // constructor
    Chunk(const uint64_t& start = 0, 
          const uint64_t& stop = 0)
        : Start(start)
        , Stop(stop)
    { }
};

bool ChunkLessThan(const Chunk& lhs, const Chunk& rhs) {
    return lhs.Start < rhs.Start;
}

typedef vector<Chunk> ChunkVector;
typedef map<uint32_t, ChunkVector> BamBinMap;
typedef vector<uint64_t> LinearOffsetVector;

struct ReferenceIndex {
    
    // data members
    BamBinMap Bins;
    LinearOffsetVector Offsets;
    bool HasAlignments;
    
    // constructor
    ReferenceIndex(const BamBinMap& binMap           = BamBinMap(),
                   const LinearOffsetVector& offsets = LinearOffsetVector(),
                   const bool hasAlignments          = false)
        : Bins(binMap)
        , Offsets(offsets)
        , HasAlignments(hasAlignments)
    { }
};

typedef map<int32_t, ReferenceIndex> BamStandardIndexData;

} // namespace BamTools
 
// -------------------------------
// BamStandardIndexPrivate implementation
  
struct BamStandardIndex::BamStandardIndexPrivate { 
  
    // parent object
    BamStandardIndex* m_parent;
    
    // data members
    BamStandardIndexData m_indexData;
    off_t m_dataBeginOffset;
    bool  m_hasFullDataCache;
    bool  m_isBigEndian;

    // ctor & dtor    
    BamStandardIndexPrivate(BamStandardIndex* parent);
    ~BamStandardIndexPrivate(void);
    
    // parent interface methods
    public:

        // creates index data (in-memory) from current reader data
        bool Build(void);
        // clear all current index offset data in memory
        void ClearAllData(void);
        // return file position after header metadata
        const off_t DataBeginOffset(void) const;
        // returns whether reference has alignments or no
        bool HasAlignments(const int& referenceID) const;
        // return true if all index data is cached
        bool HasFullDataCache(void) const;
        // attempts to use index to jump to region; returns success/fail
        bool Jump(const BamTools::BamRegion& region, bool* hasAlignmentsInRegion);
        // clears index data from all references except the first reference
        void KeepOnlyFirstReferenceOffsets(void);
        // load index data for all references, return true if loaded OK
        // @saveData - save data in memory if true, just read & discard if false
        bool LoadAllReferences(bool saveData = true);
        // load first reference from file, return true if loaded OK
        // @saveData - save data in memory if true, just read & discard if false
        bool LoadFirstReference(bool saveData = true);
        // load header data from index file, return true if loaded OK
        bool LoadHeader(void);
        // position file pointer to first reference begin, return true if skipped OK
        bool SkipToFirstReference(void);
        // write header to new index file
        bool WriteHeader(void);
        // write index data for all references to new index file
        bool WriteAllReferences(void);
    
    // internal methods
    private:

        // -----------------------
        // index file operations

        // check index file magic number, return true if OK
        bool CheckMagicNumber(void);
        // check index file version, return true if OK
        bool CheckVersion(void);
        // load a single index bin entry from file, return true if loaded OK
        // @saveData - save data in memory if true, just read & discard if false
        bool LoadBin(ReferenceIndex& refEntry, bool saveData = true);
        bool LoadBins(ReferenceIndex& refEntry, bool saveData = true);
        // load a single index bin entry from file, return true if loaded OK
        // @saveData - save data in memory if true, just read & discard if false
        bool LoadChunk(ChunkVector& chunks, bool saveData = true);
        bool LoadChunks(ChunkVector& chunks, bool saveData = true);
        // load a single index linear offset entry from file, return true if loaded OK
        // @saveData - save data in memory if true, just read & discard if false
        bool LoadLinearOffsets(ReferenceIndex& refEntry, bool saveData = true);
        // load a single reference from file, return true if loaded OK
        // @saveData - save data in memory if true, just read & discard if false
        bool LoadReference(const int& refId, bool saveData = true);
        // loads number of references, return true if loaded OK
        bool LoadReferenceCount(int& numReferences);
        // position file pointer to desired reference begin, return true if skipped OK
        bool SkipToReference(const int& refId);
        // write index data for bin to new index file
        bool WriteBin(const uint32_t& binId, const ChunkVector& chunks);
        // write index data for bins to new index file
        bool WriteBins(const BamBinMap& bins);
        // write index data for chunk entry to new index file
        bool WriteChunk(const Chunk& chunk);
        // write index data for chunk entry to new index file
        bool WriteChunks(const ChunkVector& chunks);
        // write index data for linear offsets entry to new index file
        bool WriteLinearOffsets(const LinearOffsetVector& offsets);
        // write index data single reference to new index file
        bool WriteReference(const ReferenceIndex& refEntry);

        // -----------------------
        // index data operations

        // calculate bins that overlap region
        int BinsFromRegion(const BamRegion& region,
                           const bool isRightBoundSpecified,
                           uint16_t bins[BamTools::MAX_BIN]);
        // clear all index offset data for desired reference
        void ClearReferenceOffsets(const int& refId);
        // calculates offset(s) for a given region
        bool GetOffsets(const BamRegion& region,
                        const bool isRightBoundSpecified,
                        vector<int64_t>& offsets,
                        bool* hasAlignmentsInRegion);
        // returns true if index cache has data for desired reference
        bool IsDataLoaded(const int& refId) const;
        // clears index data from all references except the one specified
        void KeepOnlyReferenceOffsets(const int& refId);
        // simplifies index by merging 'chunks'
        void MergeChunks(void);
        // saves BAM bin entry for index
        void SaveBinEntry(BamBinMap& binMap,
                          const uint32_t& saveBin,
                          const uint64_t& saveOffset,
                          const uint64_t& lastOffset);
        // saves linear offset entry for index
        void SaveLinearOffset(LinearOffsetVector& offsets,
                              const BamAlignment& bAlignment,
                              const uint64_t& lastOffset);
        // initializes index data structure to hold @count references
        void SetReferenceCount(const int& count);

};

BamStandardIndex::BamStandardIndexPrivate::BamStandardIndexPrivate(BamStandardIndex* parent)
    : m_parent(parent)
    , m_dataBeginOffset(0)
    , m_hasFullDataCache(false)
{
    m_isBigEndian = BamTools::SystemIsBigEndian();
}

BamStandardIndex::BamStandardIndexPrivate::~BamStandardIndexPrivate(void) {
    ClearAllData();
}

// calculate bins that overlap region
int BamStandardIndex::BamStandardIndexPrivate::BinsFromRegion(const BamRegion& region,
                                                              const bool isRightBoundSpecified,
                                                              uint16_t bins[MAX_BIN])
{
    // get region boundaries
    uint32_t begin = (unsigned int)region.LeftPosition;
    uint32_t end;
    
    // if right bound specified AND left&right bounds are on same reference
    // OK to use right bound position
    if ( isRightBoundSpecified && ( region.LeftRefID == region.RightRefID ) )
        end = (unsigned int)region.RightPosition;
    
    // otherwise, use end of left bound reference as cutoff
    else
        end = (unsigned int)m_parent->m_references.at(region.LeftRefID).RefLength - 1;
    
    // initialize list, bin '0' always a valid bin
    int i = 0;
    bins[i++] = 0;

    // get rest of bins that contain this region
    unsigned int k;
    for (k =    1 + (begin>>26); k <=    1 + (end>>26); ++k) { bins[i++] = k; }
    for (k =    9 + (begin>>23); k <=    9 + (end>>23); ++k) { bins[i++] = k; }
    for (k =   73 + (begin>>20); k <=   73 + (end>>20); ++k) { bins[i++] = k; }
    for (k =  585 + (begin>>17); k <=  585 + (end>>17); ++k) { bins[i++] = k; }
    for (k = 4681 + (begin>>14); k <= 4681 + (end>>14); ++k) { bins[i++] = k; }

    // return number of bins stored
    return i;
}

// creates index data (in-memory) from current reader data
bool BamStandardIndex::BamStandardIndexPrivate::Build(void) { 
  
    // localize parent data
    if ( m_parent == 0 ) return false;
    BamReader* reader = m_parent->m_reader;
    BgzfData*  mBGZF  = m_parent->m_BGZF;
  
    // be sure reader & BGZF file are valid & open for reading
    if ( reader == 0 || mBGZF == 0 || !mBGZF->IsOpen ) 
        return false;

    // move file pointer to beginning of alignments
    reader->Rewind();

    // get reference count, reserve index space
    const int numReferences = (int)m_parent->m_references.size();
    m_indexData.clear();
    m_hasFullDataCache = false;
    SetReferenceCount(numReferences);
    
    // sets default constant for bin, ID, offset, coordinate variables
    const uint32_t defaultValue = 0xffffffffu;

    // bin data
    uint32_t saveBin(defaultValue);
    uint32_t lastBin(defaultValue);

    // reference ID data
    int32_t saveRefID(defaultValue);
    int32_t lastRefID(defaultValue);

    // offset data
    uint64_t saveOffset = mBGZF->Tell();
    uint64_t lastOffset = saveOffset;

    // coordinate data
    int32_t lastCoordinate = defaultValue;

    BamAlignment bAlignment;
    while ( reader->GetNextAlignmentCore(bAlignment) ) {

        // change of chromosome, save ID, reset bin
        if ( lastRefID != bAlignment.RefID ) {
            lastRefID = bAlignment.RefID;
            lastBin   = defaultValue;
        }

        // if lastCoordinate greater than BAM position - file not sorted properly
        else if ( lastCoordinate > bAlignment.Position ) {
            fprintf(stderr, "BAM file not properly sorted:\n");
            fprintf(stderr, "Alignment %s : %d > %d on reference (id = %d)", bAlignment.Name.c_str(),
                    lastCoordinate, bAlignment.Position, bAlignment.RefID);
            exit(1);
        }

        // if valid reference && BAM bin spans some minimum cutoff (smaller bin ids span larger regions)
        if ( (bAlignment.RefID >= 0) && (bAlignment.Bin < 4681) ) {

            // save linear offset entry (matched to BAM entry refID)
	    BamStandardIndexData::iterator indexIter = m_indexData.find(bAlignment.RefID);
	    if ( indexIter == m_indexData.end() ) return false; // error
	    ReferenceIndex& refIndex = (*indexIter).second;
            LinearOffsetVector& offsets = refIndex.Offsets;
            SaveLinearOffset(offsets, bAlignment, lastOffset);
        }

        // if current BamAlignment bin != lastBin, "then possibly write the binning index"
        if ( bAlignment.Bin != lastBin ) {

            // if not first time through
            if ( saveBin != defaultValue ) {

                // save Bam bin entry
		BamStandardIndexData::iterator indexIter = m_indexData.find(saveRefID);
		if ( indexIter == m_indexData.end() ) return false; // error
		ReferenceIndex& refIndex = (*indexIter).second;
                BamBinMap& binMap = refIndex.Bins;
                SaveBinEntry(binMap, saveBin, saveOffset, lastOffset);
            }

            // update saveOffset
            saveOffset = lastOffset;

            // update bin values
            saveBin = bAlignment.Bin;
            lastBin = bAlignment.Bin;

            // update saveRefID
            saveRefID = bAlignment.RefID;

            // if invalid RefID, break out 
            if ( saveRefID < 0 ) break; 
        }

        // make sure that current file pointer is beyond lastOffset
        if ( mBGZF->Tell() <= (int64_t)lastOffset ) {
            fprintf(stderr, "Error in BGZF offsets.\n");
            exit(1);
        }

        // update lastOffset
        lastOffset = mBGZF->Tell();

        // update lastCoordinate
        lastCoordinate = bAlignment.Position;
    }

    // save any leftover BAM data (as long as refID is valid)
    if ( saveRefID >= 0 ) {
        // save Bam bin entry
	BamStandardIndexData::iterator indexIter = m_indexData.find(saveRefID);
	if ( indexIter == m_indexData.end() ) return false; // error
	ReferenceIndex& refIndex = (*indexIter).second;
        BamBinMap& binMap = refIndex.Bins;
        SaveBinEntry(binMap, saveBin, saveOffset, lastOffset);
    }

    // simplify index by merging chunks
    MergeChunks();

    // iterate through references in index
    // sort offsets in linear offset vector
    BamStandardIndexData::iterator indexIter = m_indexData.begin();
    BamStandardIndexData::iterator indexEnd  = m_indexData.end();
    for ( int i = 0; indexIter != indexEnd; ++indexIter, ++i ) {

        // get reference index data
        ReferenceIndex& refIndex = (*indexIter).second;
        LinearOffsetVector& offsets = refIndex.Offsets;

        // sort linear offsets
        sort(offsets.begin(), offsets.end());
    }

    // rewind file pointer to beginning of alignments, return success/fail
    return reader->Rewind();
}

// check index file magic number, return true if OK
bool BamStandardIndex::BamStandardIndexPrivate::CheckMagicNumber(void) {

    // read in magic number
    char magic[4];
    size_t elementsRead = fread(magic, sizeof(char), 4, m_parent->m_indexStream);

    // compare to expected value
    if ( strncmp(magic, "BAI\1", 4) != 0 ) {
        fprintf(stderr, "Problem with index file - invalid format.\n");
        fclose(m_parent->m_indexStream);
        return false;
    }

    // return success/failure of load
    return (elementsRead == 4);
}

// clear all current index offset data in memory
void BamStandardIndex::BamStandardIndexPrivate::ClearAllData(void) {
    BamStandardIndexData::const_iterator indexIter = m_indexData.begin();
    BamStandardIndexData::const_iterator indexEnd  = m_indexData.end();
    for ( ; indexIter != indexEnd; ++indexIter ) {
        const int& refId = (*indexIter).first;
        ClearReferenceOffsets(refId);
    }
}

// clear all index offset data for desired reference
void BamStandardIndex::BamStandardIndexPrivate::ClearReferenceOffsets(const int& refId) {

    // look up refId, skip if not found
    BamStandardIndexData::iterator indexIter = m_indexData.find(refId);
    if ( indexIter == m_indexData.end() ) return ;

    // clear reference data
    ReferenceIndex& refEntry = (*indexIter).second;
    refEntry.Bins.clear();
    refEntry.Offsets.clear();

    // set flag
    m_hasFullDataCache = false;
}

// return file position after header metadata
const off_t BamStandardIndex::BamStandardIndexPrivate::DataBeginOffset(void) const {
    return m_dataBeginOffset;
}

// calculates offset(s) for a given region
bool BamStandardIndex::BamStandardIndexPrivate::GetOffsets(const BamRegion& region,
                                                           const bool isRightBoundSpecified,
                                                           vector<int64_t>& offsets,
                                                           bool* hasAlignmentsInRegion)
{
    // return false if leftBound refID is not found in index data
    if ( m_indexData.find(region.LeftRefID) == m_indexData.end() )
        return false;

    // load index data for region if not already cached
    if ( !IsDataLoaded(region.LeftRefID) ) {
        bool loadedOk = true;
        loadedOk &= SkipToReference(region.LeftRefID);
        loadedOk &= LoadReference(region.LeftRefID);
        if ( !loadedOk ) return false;
    }

    // calculate which bins overlap this region
    uint16_t* bins = (uint16_t*)calloc(MAX_BIN, 2);
    int numBins = BinsFromRegion(region, isRightBoundSpecified, bins);

    // get bins for this reference
    BamStandardIndexData::const_iterator indexIter = m_indexData.find(region.LeftRefID);
    if ( indexIter == m_indexData.end() ) return false; // error
    const ReferenceIndex& refIndex = (*indexIter).second;
    const BamBinMap& binMap        = refIndex.Bins;

    // get minimum offset to consider
    const LinearOffsetVector& linearOffsets = refIndex.Offsets;
    const uint64_t minOffset = ( (unsigned int)(region.LeftPosition>>BAM_LIDX_SHIFT) >= linearOffsets.size() )
                               ? 0 : linearOffsets.at(region.LeftPosition>>BAM_LIDX_SHIFT);

    // store all alignment 'chunk' starts (file offsets) for bins in this region
    for ( int i = 0; i < numBins; ++i ) {
      
        const uint16_t binKey = bins[i];
        map<uint32_t, ChunkVector>::const_iterator binIter = binMap.find(binKey);
        if ( (binIter != binMap.end()) && ((*binIter).first == binKey) ) {

            // iterate over chunks
            const ChunkVector& chunks = (*binIter).second;
            std::vector<Chunk>::const_iterator chunksIter = chunks.begin();
            std::vector<Chunk>::const_iterator chunksEnd  = chunks.end();
            for ( ; chunksIter != chunksEnd; ++chunksIter) {
              
                // if valid chunk found, store its file offset
                const Chunk& chunk = (*chunksIter);
                if ( chunk.Stop > minOffset )
                    offsets.push_back( chunk.Start );
            }
        }
    }

    // clean up memory
    free(bins);

    // sort the offsets before returning
    sort(offsets.begin(), offsets.end());
    
    // set flag & return success
    *hasAlignmentsInRegion = (offsets.size() != 0 );

    // if cache mode set to none, dump the data we just loaded
    if (m_parent->m_cacheMode == BamIndex::NoIndexCaching )
        ClearReferenceOffsets(region.LeftRefID);

    // return succes
    return true;
}

// returns whether reference has alignments or no
bool BamStandardIndex::BamStandardIndexPrivate::HasAlignments(const int& refId) const {
    BamStandardIndexData::const_iterator indexIter = m_indexData.find(refId);
    if ( indexIter == m_indexData.end() ) return false; // error
    const ReferenceIndex& refEntry = (*indexIter).second;
    return refEntry.HasAlignments;
}

// return true if all index data is cached
bool BamStandardIndex::BamStandardIndexPrivate::HasFullDataCache(void) const {
    return m_hasFullDataCache;
}

// returns true if index cache has data for desired reference
bool BamStandardIndex::BamStandardIndexPrivate::IsDataLoaded(const int& refId) const {

    // look up refId, return false if not found
    BamStandardIndexData::const_iterator indexIter = m_indexData.find(refId);
    if ( indexIter == m_indexData.end() ) return false;

    // see if reference has alignments
    // if not, it's not a problem to have no offset data
    const ReferenceIndex& refEntry = (*indexIter).second;
    if ( !refEntry.HasAlignments ) return true;

    // return whether bin map contains data
    return ( !refEntry.Bins.empty() );
}

// attempts to use index to jump to region; returns success/fail
bool BamStandardIndex::BamStandardIndexPrivate::Jump(const BamRegion& region, bool* hasAlignmentsInRegion) {
  
    // localize parent data
    if ( m_parent == 0 ) return false;
    BamReader* reader     = m_parent->m_reader;
    BgzfData*  mBGZF      = m_parent->m_BGZF;
    RefVector& references = m_parent->m_references;
  
    // be sure reader & BGZF file are valid & open for reading
    if ( reader == 0 || mBGZF == 0 || !mBGZF->IsOpen ) 
        return false;
    
    // make sure left-bound position is valid
    if ( region.LeftPosition > references.at(region.LeftRefID).RefLength )
        return false;
        
    // calculate offsets for this region
    // if failed, print message, set flag, and return failure
    vector<int64_t> offsets;
    if ( !GetOffsets(region, region.isRightBoundSpecified(), offsets, hasAlignmentsInRegion) ) {
        fprintf(stderr, "ERROR: Could not jump: unable to calculate offset(s) for specified region.\n");
        *hasAlignmentsInRegion = false;
        return false;
    }
  
    // iterate through offsets
    BamAlignment bAlignment;
    bool result = true;
    for ( vector<int64_t>::const_iterator o = offsets.begin(); o != offsets.end(); ++o) {
        
        // attempt seek & load first available alignment
        // set flag to true if data exists
        result &= mBGZF->Seek(*o);
        *hasAlignmentsInRegion = reader->GetNextAlignmentCore(bAlignment);
        
        // if this alignment corresponds to desired position
        // return success of seeking back to the offset before the 'current offset' (to cover overlaps)
        if ( ((bAlignment.RefID == region.LeftRefID) &&
               ((bAlignment.Position + bAlignment.Length) > region.LeftPosition)) ||
             (bAlignment.RefID > region.LeftRefID) )
        {
            if ( o != offsets.begin() ) --o;
            return mBGZF->Seek(*o);
        }
    }
    
    // if error in jumping, print message & set flag
    if ( !result ) {
        fprintf(stderr, "ERROR: Could not jump: unable to determine correct offset for specified region.\n");
        *hasAlignmentsInRegion = false;
    }
    
    // return success/failure
    return result;
}

// clears index data from all references except the first
void BamStandardIndex::BamStandardIndexPrivate::KeepOnlyFirstReferenceOffsets(void) {
    BamStandardIndexData::const_iterator indexBegin = m_indexData.begin();
    KeepOnlyReferenceOffsets((*indexBegin).first);
}

// clears index data from all references except the one specified
void BamStandardIndex::BamStandardIndexPrivate::KeepOnlyReferenceOffsets(const int& refId) {
    BamStandardIndexData::iterator mapIter = m_indexData.begin();
    BamStandardIndexData::iterator mapEnd  = m_indexData.end();
    for ( ; mapIter != mapEnd; ++mapIter ) {
        const int entryRefId = (*mapIter).first;
        if ( entryRefId != refId )
            ClearReferenceOffsets(entryRefId);
    }
}

bool BamStandardIndex::BamStandardIndexPrivate::LoadAllReferences(bool saveData) {

    // skip if data already loaded
    if ( m_hasFullDataCache ) return true;

    // get number of reference sequences
    uint32_t numReferences;
    if ( !LoadReferenceCount((int&)numReferences) )
        return false;

    // iterate over reference entries
    bool loadedOk = true;
    for ( int i = 0; i < (int)numReferences; ++i )
        loadedOk &= LoadReference(i, saveData);

    // set flag
    if ( loadedOk && saveData )
        m_hasFullDataCache = true;

    // return success/failure of loading references
    return loadedOk;
}

// load header data from index file, return true if loaded OK
bool BamStandardIndex::BamStandardIndexPrivate::LoadHeader(void) {

    bool loadedOk = CheckMagicNumber();

    // store offset of beginning of data
    m_dataBeginOffset = ftell64(m_parent->m_indexStream);

    // return success/failure of load
    return loadedOk;
}

// load a single index bin entry from file, return true if loaded OK
// @saveData - save data in memory if true, just read & discard if false
bool BamStandardIndex::BamStandardIndexPrivate::LoadBin(ReferenceIndex& refEntry, bool saveData) {

    size_t elementsRead = 0;

    // get bin ID
    uint32_t binId;
    elementsRead += fread(&binId, sizeof(binId), 1, m_parent->m_indexStream);
    if ( m_isBigEndian ) SwapEndian_32(binId);

    // load alignment chunks for this bin
    ChunkVector chunks;
    bool chunksOk = LoadChunks(chunks, saveData);

    // store bin entry
    if ( chunksOk && saveData )
        refEntry.Bins.insert(pair<uint32_t, ChunkVector>(binId, chunks));

    // return success/failure of load
    return ( (elementsRead == 1) && chunksOk );
}

bool BamStandardIndex::BamStandardIndexPrivate::LoadBins(ReferenceIndex& refEntry, bool saveData) {

    size_t elementsRead = 0;

    // get number of bins
    int32_t numBins;
    elementsRead += fread(&numBins, sizeof(numBins), 1, m_parent->m_indexStream);
    if ( m_isBigEndian ) SwapEndian_32(numBins);

    // set flag
    refEntry.HasAlignments = ( numBins != 0 );

    // iterate over bins
    bool binsOk = true;
    for ( int i = 0; i < numBins; ++i )
        binsOk &= LoadBin(refEntry, saveData);

    // return success/failure of load
    return ( (elementsRead == 1) && binsOk );
}

// load a single index bin entry from file, return true if loaded OK
// @saveData - save data in memory if true, just read & discard if false
bool BamStandardIndex::BamStandardIndexPrivate::LoadChunk(ChunkVector& chunks, bool saveData) {

    size_t elementsRead = 0;

    // read in chunk data
    uint64_t start;
    uint64_t stop;
    elementsRead += fread(&start, sizeof(start), 1, m_parent->m_indexStream);
    elementsRead += fread(&stop,  sizeof(stop),  1, m_parent->m_indexStream);

    // swap endian-ness if necessary
    if ( m_isBigEndian ) {
        SwapEndian_64(start);
        SwapEndian_64(stop);
    }

    // save data if requested
    if ( saveData ) chunks.push_back( Chunk(start, stop) );

    // return success/failure of load
    return ( elementsRead == 2 );
}

bool BamStandardIndex::BamStandardIndexPrivate::LoadChunks(ChunkVector& chunks, bool saveData) {

    size_t elementsRead = 0;

    // read in number of chunks
    uint32_t numChunks;
    elementsRead += fread(&numChunks, sizeof(numChunks), 1, m_parent->m_indexStream);
    if ( m_isBigEndian ) SwapEndian_32(numChunks);

    // initialize space for chunks if we're storing this data
    if ( saveData ) chunks.reserve(numChunks);

    // iterate over chunks
    bool chunksOk = true;
    for ( int i = 0; i < (int)numChunks; ++i )
        chunksOk &= LoadChunk(chunks, saveData);

    // sort chunk vector
    sort( chunks.begin(), chunks.end(), ChunkLessThan );

    // return success/failure of load
    return ( (elementsRead == 1) && chunksOk );
}

// load a single index linear offset entry from file, return true if loaded OK
// @saveData - save data in memory if true, just read & discard if false
bool BamStandardIndex::BamStandardIndexPrivate::LoadLinearOffsets(ReferenceIndex& refEntry, bool saveData) {

    size_t elementsRead = 0;

    // read in number of linear offsets
    int32_t numLinearOffsets;
    elementsRead += fread(&numLinearOffsets, sizeof(numLinearOffsets), 1, m_parent->m_indexStream);
    if ( m_isBigEndian ) SwapEndian_32(numLinearOffsets);

    // set up destination vector (if we're saving the data)
    LinearOffsetVector linearOffsets;
    if ( saveData ) linearOffsets.reserve(numLinearOffsets);

    // iterate over linear offsets
    uint64_t linearOffset;
    for ( int i = 0; i < numLinearOffsets; ++i ) {
        elementsRead += fread(&linearOffset, sizeof(linearOffset), 1, m_parent->m_indexStream);
        if ( m_isBigEndian ) SwapEndian_64(linearOffset);
        if ( saveData ) linearOffsets.push_back(linearOffset);
    }

    // sort linear offsets
    sort ( linearOffsets.begin(), linearOffsets.end() );

    // save in reference index entry if desired
    if ( saveData ) refEntry.Offsets = linearOffsets;

    // return success/failure of load
    return ( elementsRead == (size_t)(numLinearOffsets + 1) );
}

bool BamStandardIndex::BamStandardIndexPrivate::LoadFirstReference(bool saveData) {
    BamStandardIndexData::const_iterator indexBegin = m_indexData.begin();
    return LoadReference((*indexBegin).first, saveData);
}

// load a single reference from file, return true if loaded OK
// @saveData - save data in memory if true, just read & discard if false
bool BamStandardIndex::BamStandardIndexPrivate::LoadReference(const int& refId, bool saveData) {    

    // look up refId
    BamStandardIndexData::iterator indexIter = m_indexData.find(refId);

    // if reference not previously loaded, create new entry
    if ( indexIter == m_indexData.end() ) {
        ReferenceIndex newEntry;
        newEntry.HasAlignments = false;
        m_indexData.insert( pair<int32_t, ReferenceIndex>(refId, newEntry) );
    }

    // load reference data
    indexIter = m_indexData.find(refId);
    ReferenceIndex& entry = (*indexIter).second;
    bool loadedOk = true;
    loadedOk &= LoadBins(entry, saveData);
    loadedOk &= LoadLinearOffsets(entry, saveData);
    return loadedOk;
}

// loads number of references, return true if loaded OK
bool BamStandardIndex::BamStandardIndexPrivate::LoadReferenceCount(int& numReferences) {

    size_t elementsRead = 0;

    // read reference count
    elementsRead += fread(&numReferences, sizeof(numReferences), 1, m_parent->m_indexStream);
    if ( m_isBigEndian ) SwapEndian_32(numReferences);

    // return success/failure of load
    return ( elementsRead == 1 );
}

// merges 'alignment chunks' in BAM bin (used for index building)
void BamStandardIndex::BamStandardIndexPrivate::MergeChunks(void) {

    // iterate over reference enties
    BamStandardIndexData::iterator indexIter = m_indexData.begin();
    BamStandardIndexData::iterator indexEnd  = m_indexData.end();
    for ( ; indexIter != indexEnd; ++indexIter ) {

        // get BAM bin map for this reference
        ReferenceIndex& refIndex = (*indexIter).second;
        BamBinMap& bamBinMap = refIndex.Bins;

        // iterate over BAM bins
        BamBinMap::iterator binIter = bamBinMap.begin();
        BamBinMap::iterator binEnd  = bamBinMap.end();
        for ( ; binIter != binEnd; ++binIter ) {

            // get chunk vector for this bin
            ChunkVector& binChunks = (*binIter).second;
            if ( binChunks.size() == 0 ) continue; 

            ChunkVector mergedChunks;
            mergedChunks.push_back( binChunks[0] );

            // iterate over chunks
            int i = 0;
            ChunkVector::iterator chunkIter = binChunks.begin();
            ChunkVector::iterator chunkEnd  = binChunks.end();
            for ( ++chunkIter; chunkIter != chunkEnd; ++chunkIter) {

                // get 'currentChunk' based on numeric index
                Chunk& currentChunk = mergedChunks[i];

                // get iteratorChunk based on vector iterator
                Chunk& iteratorChunk = (*chunkIter);

                // if chunk ends where (iterator) chunk starts, then merge
                if ( currentChunk.Stop>>16 == iteratorChunk.Start>>16 )
                    currentChunk.Stop = iteratorChunk.Stop;

                // otherwise
                else {
                    // set currentChunk + 1 to iteratorChunk
                    mergedChunks.push_back(iteratorChunk);
                    ++i;
                }
            }

            // saved merged chunk vector
            (*binIter).second = mergedChunks;
        }
    }
}

// saves BAM bin entry for index
void BamStandardIndex::BamStandardIndexPrivate::SaveBinEntry(BamBinMap& binMap,
                                                             const uint32_t& saveBin,
                                                             const uint64_t& saveOffset,
                                                             const uint64_t& lastOffset)
{
    // look up saveBin
    BamBinMap::iterator binIter = binMap.find(saveBin);

    // create new chunk
    Chunk newChunk(saveOffset, lastOffset);

    // if entry doesn't exist
    if ( binIter == binMap.end() ) {
        ChunkVector newChunks;
        newChunks.push_back(newChunk);
        binMap.insert( pair<uint32_t, ChunkVector>(saveBin, newChunks));
    }

    // otherwise
    else {
        ChunkVector& binChunks = (*binIter).second;
        binChunks.push_back( newChunk );
    }
}

// saves linear offset entry for index
void BamStandardIndex::BamStandardIndexPrivate::SaveLinearOffset(LinearOffsetVector& offsets,
                                                                 const BamAlignment& bAlignment,
                                                                 const uint64_t&     lastOffset)
{
    // get converted offsets
    int beginOffset = bAlignment.Position >> BAM_LIDX_SHIFT;
    int endOffset   = (bAlignment.GetEndPosition() - 1) >> BAM_LIDX_SHIFT;

    // resize vector if necessary
    int oldSize = offsets.size();
    int newSize = endOffset + 1;
    if ( oldSize < newSize )
        offsets.resize(newSize, 0);

    // store offset
    for( int i = beginOffset + 1; i <= endOffset; ++i ) {
        if ( offsets[i] == 0 )
            offsets[i] = lastOffset;
    }
}

// initializes index data structure to hold @count references
void BamStandardIndex::BamStandardIndexPrivate::SetReferenceCount(const int& count) {
    for ( int i = 0; i < count; ++i )
        m_indexData[i].HasAlignments = false;
}

bool BamStandardIndex::BamStandardIndexPrivate::SkipToFirstReference(void) {
    BamStandardIndexData::const_iterator indexBegin = m_indexData.begin();
    return SkipToReference( (*indexBegin).first );
}

// position file pointer to desired reference begin, return true if skipped OK
bool BamStandardIndex::BamStandardIndexPrivate::SkipToReference(const int& refId) {

    // attempt rewind
    if ( !m_parent->Rewind() ) return false;

    // read in number of references
    uint32_t numReferences;
    size_t elementsRead = fread(&numReferences, sizeof(numReferences), 1, m_parent->m_indexStream);
    if ( elementsRead != 1 ) return false;
    if ( m_isBigEndian ) SwapEndian_32(numReferences);

    // iterate over reference entries
    bool skippedOk = true;
    int currentRefId = 0;
    while (currentRefId != refId) {
        skippedOk &= LoadReference(currentRefId, false);
        ++currentRefId;
    }

    // return success
    return skippedOk;
}
 
// write header to new index file
bool BamStandardIndex::BamStandardIndexPrivate::WriteHeader(void) {

    size_t elementsWritten = 0;

    // write magic number
    elementsWritten += fwrite("BAI\1", sizeof(char), 4, m_parent->m_indexStream);

    // store offset of beginning of data
    m_dataBeginOffset = ftell64(m_parent->m_indexStream);

    // return success/failure of write
    return (elementsWritten == 4);
}

// write index data for all references to new index file
bool BamStandardIndex::BamStandardIndexPrivate::WriteAllReferences(void) {

    size_t elementsWritten = 0;

    // write number of reference sequences
    int32_t numReferenceSeqs = m_indexData.size();
    if ( m_isBigEndian ) SwapEndian_32(numReferenceSeqs);
    elementsWritten += fwrite(&numReferenceSeqs, sizeof(numReferenceSeqs), 1, m_parent->m_indexStream);

    // iterate over reference sequences
    bool refsOk = true;
    BamStandardIndexData::const_iterator indexIter = m_indexData.begin();
    BamStandardIndexData::const_iterator indexEnd  = m_indexData.end();
    for ( ; indexIter != indexEnd; ++ indexIter )
        refsOk &= WriteReference( (*indexIter).second );

    // return success/failure of write
    return ( (elementsWritten == 1) && refsOk );
}

// write index data for bin to new index file
bool BamStandardIndex::BamStandardIndexPrivate::WriteBin(const uint32_t& binId, const ChunkVector& chunks) {

    size_t elementsWritten = 0;

    // write BAM bin ID
    uint32_t binKey = binId;
    if ( m_isBigEndian ) SwapEndian_32(binKey);
    elementsWritten += fwrite(&binKey, sizeof(binKey), 1, m_parent->m_indexStream);

    // write chunks
    bool chunksOk = WriteChunks(chunks);

    // return success/failure of write
    return ( (elementsWritten == 1) && chunksOk );
}

// write index data for bins to new index file
bool BamStandardIndex::BamStandardIndexPrivate::WriteBins(const BamBinMap& bins) {

    size_t elementsWritten = 0;

    // write number of bins
    int32_t binCount = bins.size();
    if ( m_isBigEndian ) SwapEndian_32(binCount);
    elementsWritten += fwrite(&binCount, sizeof(binCount), 1, m_parent->m_indexStream);

    // iterate over bins
    bool binsOk = true;
    BamBinMap::const_iterator binIter = bins.begin();
    BamBinMap::const_iterator binEnd  = bins.end();
    for ( ; binIter != binEnd; ++binIter )
        binsOk &= WriteBin( (*binIter).first, (*binIter).second );

    // return success/failure of write
    return ( (elementsWritten == 1) && binsOk );
}

// write index data for chunk entry to new index file
bool BamStandardIndex::BamStandardIndexPrivate::WriteChunk(const Chunk& chunk) {

    size_t elementsWritten = 0;

    // localize alignment chunk offsets
    uint64_t start = chunk.Start;
    uint64_t stop  = chunk.Stop;

    // swap endian-ness if necessary
    if ( m_isBigEndian ) {
        SwapEndian_64(start);
        SwapEndian_64(stop);
    }

    // write to index file
    elementsWritten += fwrite(&start, sizeof(start), 1, m_parent->m_indexStream);
    elementsWritten += fwrite(&stop,  sizeof(stop),  1, m_parent->m_indexStream);

    // return success/failure of write
    return ( elementsWritten == 2 );
}

// write index data for chunk entry to new index file
bool BamStandardIndex::BamStandardIndexPrivate::WriteChunks(const ChunkVector& chunks) {

    size_t elementsWritten = 0;

    // write chunks
    int32_t chunkCount = chunks.size();
    if ( m_isBigEndian ) SwapEndian_32(chunkCount);
    elementsWritten += fwrite(&chunkCount, sizeof(chunkCount), 1, m_parent->m_indexStream);

    // iterate over chunks
    bool chunksOk = true;
    ChunkVector::const_iterator chunkIter = chunks.begin();
    ChunkVector::const_iterator chunkEnd  = chunks.end();
    for ( ; chunkIter != chunkEnd; ++chunkIter )
        chunksOk &= WriteChunk( (*chunkIter) );

    // return success/failure of write
    return ( (elementsWritten == 1) && chunksOk );
}

// write index data for linear offsets entry to new index file
bool BamStandardIndex::BamStandardIndexPrivate::WriteLinearOffsets(const LinearOffsetVector& offsets) {

    size_t elementsWritten = 0;

    // write number of linear offsets
    int32_t offsetCount = offsets.size();
    if ( m_isBigEndian ) SwapEndian_32(offsetCount);
    elementsWritten += fwrite(&offsetCount, sizeof(offsetCount), 1, m_parent->m_indexStream);

    // iterate over linear offsets
    LinearOffsetVector::const_iterator offsetIter = offsets.begin();
    LinearOffsetVector::const_iterator offsetEnd  = offsets.end();
    for ( ; offsetIter != offsetEnd; ++offsetIter ) {

        // write linear offset
        uint64_t linearOffset = (*offsetIter);
        if ( m_isBigEndian ) SwapEndian_64(linearOffset);
        elementsWritten += fwrite(&linearOffset, sizeof(linearOffset), 1, m_parent->m_indexStream);
    }

    // return success/failure of write
    return ( elementsWritten == (size_t)(offsetCount + 1) );
}

// write index data for a single reference to new index file
bool BamStandardIndex::BamStandardIndexPrivate::WriteReference(const ReferenceIndex& refEntry) {
    bool refOk = true;
    refOk &= WriteBins(refEntry.Bins);
    refOk &= WriteLinearOffsets(refEntry.Offsets);
    return refOk;
}

// ---------------------------------------------------------------
// BamStandardIndex implementation
 
BamStandardIndex::BamStandardIndex(BgzfData* bgzf, BamReader* reader)
    : BamIndex(bgzf, reader)
{
    d = new BamStandardIndexPrivate(this);
}    

BamStandardIndex::~BamStandardIndex(void) {
    delete d;
    d = 0;
}

// BamIndex interface implementation
bool BamStandardIndex::Build(void) { return d->Build(); }
void BamStandardIndex::ClearAllData(void) { d->ClearAllData(); }
const off_t BamStandardIndex::DataBeginOffset(void) const { return d->DataBeginOffset(); }
bool BamStandardIndex::HasAlignments(const int& referenceID) const { return d->HasAlignments(referenceID); }
bool BamStandardIndex::HasFullDataCache(void) const { return d->HasFullDataCache(); }
bool BamStandardIndex::Jump(const BamRegion& region, bool* hasAlignmentsInRegion) { return d->Jump(region, hasAlignmentsInRegion); }
void BamStandardIndex::KeepOnlyFirstReferenceOffsets(void) { d->KeepOnlyFirstReferenceOffsets(); }
bool BamStandardIndex::LoadAllReferences(bool saveData) { return d->LoadAllReferences(saveData); }
bool BamStandardIndex::LoadFirstReference(bool saveData) { return d->LoadFirstReference(saveData); }
bool BamStandardIndex::LoadHeader(void) { return d->LoadHeader(); }
bool BamStandardIndex::SkipToFirstReference(void) { return d->SkipToFirstReference(); }
bool BamStandardIndex::WriteAllReferences(void) { return d->WriteAllReferences(); }
bool BamStandardIndex::WriteHeader(void) { return d->WriteHeader(); }

// #########################################################################################
// #########################################################################################

// ---------------------------------------------------
// BamToolsIndex structs & typedefs

namespace BamTools {
  
// individual index offset entry
struct BamToolsIndexEntry {
    
    // data members
    int32_t MaxEndPosition;
    int64_t StartOffset;
    int32_t StartPosition;
    
    // ctor
    BamToolsIndexEntry(const int32_t& maxEndPosition = 0,
                       const int64_t& startOffset    = 0,
                       const int32_t& startPosition  = 0)
        : MaxEndPosition(maxEndPosition)
        , StartOffset(startOffset)
        , StartPosition(startPosition)
    { }
};

// reference index entry
struct BamToolsReferenceEntry {

    // data members
    bool HasAlignments;
    vector<BamToolsIndexEntry> Offsets;

    // ctor
    BamToolsReferenceEntry(void)
        : HasAlignments(false)
    { }
};

// the actual index data structure
typedef map<int, BamToolsReferenceEntry> BamToolsIndexData;
  
} // namespace BamTools

// ---------------------------------------------------
// BamToolsIndexPrivate implementation

struct BamToolsIndex::BamToolsIndexPrivate {
    
    // keep a list of any supported versions here 
    // (might be useful later to handle any 'legacy' versions if the format changes)
    // listed for example like: BTI_1_0 = 1, BTI_1_1 = 2, BTI_1_2 = 3, BTI_2_0 = 4, and so on
    //
    // so a change introduced in (hypothetical) BTI_1_2 would be handled from then on by: 
    //
    // if ( indexVersion >= BTI_1_2 ) 
    //   do something new 
    // else 
    //   do the old thing
    enum Version { BTI_1_0 = 1 
                 , BTI_1_1
                 , BTI_1_2
                 };  
  
    // parent object
    BamToolsIndex* m_parent;
    
    // data members
    int32_t           m_blockSize;
    BamToolsIndexData m_indexData;
    off_t             m_dataBeginOffset;
    bool              m_hasFullDataCache;
    bool              m_isBigEndian;
    int32_t           m_inputVersion; // Version is serialized as int
    Version           m_outputVersion;
    
    // ctor & dtor    
    BamToolsIndexPrivate(BamToolsIndex* parent);
    ~BamToolsIndexPrivate(void);
    
    // parent interface methods
    public:

        // creates index data (in-memory) from current reader data
        bool Build(void);
        // clear all current index offset data in memory
        void ClearAllData(void);
        // return file position after header metadata
        const off_t DataBeginOffset(void) const;
        // returns whether reference has alignments or no
        bool HasAlignments(const int& referenceID) const;
        // return true if all index data is cached
        bool HasFullDataCache(void) const;
        // attempts to use index to jump to region; returns success/fail
        bool Jump(const BamTools::BamRegion& region, bool* hasAlignmentsInRegion);
        // clears index data from all references except the first
        void KeepOnlyFirstReferenceOffsets(void);
        // load index data for all references, return true if loaded OK
        // @saveData - save data in memory if true, just read & discard if false
        bool LoadAllReferences(bool saveData = true);
        // load first reference from file, return true if loaded OK
        // @saveData - save data in memory if true, just read & discard if false
        bool LoadFirstReference(bool saveData = true);
        // load header data from index file, return true if loaded OK
        bool LoadHeader(void);
        // position file pointer to desired reference begin, return true if skipped OK
        bool SkipToFirstReference(void);
        // write header to new index file
        bool WriteHeader(void);
        // write index data for all references to new index file
        bool WriteAllReferences(void);
    
    // internal methods
    private:

        // -----------------------
        // index file operations

        // check index file magic number, return true if OK
        bool CheckMagicNumber(void);
        // check index file version, return true if OK
        bool CheckVersion(void);
        // return true if FILE* is open
        bool IsOpen(void) const;
        // load a single index entry from file, return true if loaded OK
        // @saveData - save data in memory if true, just read & discard if false
        bool LoadIndexEntry(const int& refId, bool saveData = true);
        // load a single reference from file, return true if loaded OK
        // @saveData - save data in memory if true, just read & discard if false
        bool LoadReference(const int& refId, bool saveData = true);
        // loads number of references, return true if loaded OK
        bool LoadReferenceCount(int& numReferences);
        // position file pointer to desired reference begin, return true if skipped OK
        bool SkipToReference(const int& refId);
        // write current reference index data to new index file
        bool WriteReferenceEntry(const BamToolsReferenceEntry& refEntry);
        // write current index offset entry to new index file
        bool WriteIndexEntry(const BamToolsIndexEntry& entry);

        // -----------------------
        // index data operations

        // clear all index offset data for desired reference
        void ClearReferenceOffsets(const int& refId);
        // calculate BAM file offset for desired region
        // return true if no error (*NOT* equivalent to "has alignments or valid offset")
        //   check @hasAlignmentsInRegion to determine this status
        // @region - target region
        // @offset - resulting seek target
        // @hasAlignmentsInRegion - sometimes a file just lacks data in region, this flag indicates that status
        bool GetOffset(const BamRegion& region, int64_t& offset, bool* hasAlignmentsInRegion);
        // returns true if index cache has data for desired reference
        bool IsDataLoaded(const int& refId) const;
        // clears index data from all references except the one specified
        void KeepOnlyReferenceOffsets(const int& refId);
        // saves an index offset entry in memory
        void SaveOffsetEntry(const int& refId, const BamToolsIndexEntry& entry);
        // pre-allocates size for offset vector
        void SetOffsetCount(const int& refId, const int& offsetCount);
        // initializes index data structure to hold @count references
        void SetReferenceCount(const int& count);
};

// ctor
BamToolsIndex::BamToolsIndexPrivate::BamToolsIndexPrivate(BamToolsIndex* parent)
    : m_parent(parent)
    , m_blockSize(1000)
    , m_dataBeginOffset(0)
    , m_hasFullDataCache(false)
    , m_inputVersion(0)
    , m_outputVersion(BTI_1_2) // latest version - used for writing new index files
{
    m_isBigEndian = BamTools::SystemIsBigEndian();
}

// dtor
BamToolsIndex::BamToolsIndexPrivate::~BamToolsIndexPrivate(void) {
    ClearAllData();
}

// creates index data (in-memory) from current reader data
bool BamToolsIndex::BamToolsIndexPrivate::Build(void) { 
  
    // localize parent data
    if ( m_parent == 0 ) return false;
    BamReader* reader = m_parent->m_reader;
    BgzfData*  mBGZF  = m_parent->m_BGZF;
  
    // be sure reader & BGZF file are valid & open for reading
    if ( reader == 0 || mBGZF == 0 || !mBGZF->IsOpen ) 
        return false;

    // move file pointer to beginning of alignments
    if ( !reader->Rewind() ) return false;
    
    // initialize index data structure with space for all references
    const int numReferences = (int)m_parent->m_references.size();
    m_indexData.clear();
    m_hasFullDataCache = false;
    SetReferenceCount(numReferences);

    // set up counters and markers
    int32_t currentBlockCount      = 0;
    int64_t currentAlignmentOffset = mBGZF->Tell();
    int32_t blockRefId             = 0;
    int32_t blockMaxEndPosition    = 0;
    int64_t blockStartOffset       = currentAlignmentOffset;
    int32_t blockStartPosition     = -1;
    
    // plow through alignments, storing index entries
    BamAlignment al;
    while ( reader->GetNextAlignmentCore(al) ) {
        
        // if block contains data (not the first time through) AND alignment is on a new reference
        if ( currentBlockCount > 0 && al.RefID != blockRefId ) {
          
	    // store previous data
            BamToolsIndexEntry entry(blockMaxEndPosition, blockStartOffset, blockStartPosition);
            SaveOffsetEntry(blockRefId, entry);

            // intialize new block for current alignment's reference
            currentBlockCount   = 0;
            blockMaxEndPosition = al.GetEndPosition();
            blockStartOffset    = currentAlignmentOffset;
        }
        
        // if beginning of block, save first alignment's refID & position
        if ( currentBlockCount == 0 ) {
            blockRefId         = al.RefID;
            blockStartPosition = al.Position;
        }
      
        // increment block counter
        ++currentBlockCount;

        // check end position
        int32_t alignmentEndPosition = al.GetEndPosition();
        if ( alignmentEndPosition > blockMaxEndPosition ) 
            blockMaxEndPosition = alignmentEndPosition;
        
        // if block is full, get offset for next block, reset currentBlockCount
        if ( currentBlockCount == m_blockSize ) {
            BamToolsIndexEntry entry(blockMaxEndPosition, blockStartOffset, blockStartPosition);
            SaveOffsetEntry(blockRefId, entry);
            blockStartOffset  = mBGZF->Tell();
            currentBlockCount = 0;
        }
        
        // not the best name, but for the next iteration, this value will be the offset of the *current* alignment
        // necessary because we won't know if this next alignment is on a new reference until we actually read it
        currentAlignmentOffset = mBGZF->Tell();  
    }
    
    // store final block with data
    BamToolsIndexEntry entry(blockMaxEndPosition, blockStartOffset, blockStartPosition);
    SaveOffsetEntry(blockRefId, entry);
    
    // set flag
    m_hasFullDataCache = true;

    // return success/failure of rewind
    return reader->Rewind();
}

// check index file magic number, return true if OK
bool BamToolsIndex::BamToolsIndexPrivate::CheckMagicNumber(void) {

    // see if index is valid BAM index
    char magic[4];
    size_t elementsRead = fread(magic, 1, 4, m_parent->m_indexStream);
    if ( elementsRead != 4 ) return false;
    if ( strncmp(magic, "BTI\1", 4) != 0 ) {
        fprintf(stderr, "Problem with index file - invalid format.\n");
        return false;
    }

    // otherwise ok
    return true;
}

// check index file version, return true if OK
bool BamToolsIndex::BamToolsIndexPrivate::CheckVersion(void) {

    // read version from file
    size_t elementsRead = fread(&m_inputVersion, sizeof(m_inputVersion), 1, m_parent->m_indexStream);
    if ( elementsRead != 1 ) return false;
    if ( m_isBigEndian ) SwapEndian_32(m_inputVersion);

    // if version is negative, or zero
    if ( m_inputVersion <= 0 ) {
        fprintf(stderr, "Problem with index file - invalid version.\n");
        return false;
    }

    // if version is newer than can be supported by this version of bamtools
    else if ( m_inputVersion > m_outputVersion ) {
        fprintf(stderr, "Problem with index file - attempting to use an outdated version of BamTools with a newer index file.\n");
        fprintf(stderr, "Please update BamTools to a more recent version to support this index file.\n");
        return false;
    }

    // ------------------------------------------------------------------
    // check for deprecated, unsupported versions
    // (typically whose format did not accomodate a particular bug fix)

    else if ( (Version)m_inputVersion == BTI_1_0 ) {
        fprintf(stderr, "\nProblem with index file - this version of the index contains a bug related to accessing data near reference ends.\n");
        fprintf(stderr, "\nPlease run \'bamtools index -bti -in yourData.bam\' to generate an up-to-date BamToolsIndex.\n\n");
        return false;
    }

    else if ( (Version)m_inputVersion == BTI_1_1 ) {
        fprintf(stderr, "\nProblem with index file - this version of the index contains a bug related to handling empty references.\n");
        fprintf(stderr, "\nPlease run \'bamtools index -bti -in yourData.bam\' to generate an up-to-date BamToolsIndex.\n\n");
        return false;
    }

    // otherwise ok
    else return true;
}

// clear all current index offset data in memory
void BamToolsIndex::BamToolsIndexPrivate::ClearAllData(void) {
    BamToolsIndexData::const_iterator indexIter = m_indexData.begin();
    BamToolsIndexData::const_iterator indexEnd  = m_indexData.end();
    for ( ; indexIter != indexEnd; ++indexIter ) {
        const int& refId = (*indexIter).first;
        ClearReferenceOffsets(refId);
    }
}

// clear all index offset data for desired reference
void BamToolsIndex::BamToolsIndexPrivate::ClearReferenceOffsets(const int& refId) {
    if ( m_indexData.find(refId) == m_indexData.end() ) return;
    vector<BamToolsIndexEntry>& offsets = m_indexData[refId].Offsets;
    offsets.clear();
    m_hasFullDataCache = false;
}

// return file position after header metadata
const off_t BamToolsIndex::BamToolsIndexPrivate::DataBeginOffset(void) const {
    return m_dataBeginOffset;
}

// calculate BAM file offset for desired region
// return true if no error (*NOT* equivalent to "has alignments or valid offset")
//   check @hasAlignmentsInRegion to determine this status
// @region - target region
// @offset - resulting seek target
// @hasAlignmentsInRegion - sometimes a file just lacks data in region, this flag indicates that status
// N.B. - ignores isRightBoundSpecified
bool BamToolsIndex::BamToolsIndexPrivate::GetOffset(const BamRegion& region, int64_t& offset, bool* hasAlignmentsInRegion) { 
  
    // return false if leftBound refID is not found in index data
    BamToolsIndexData::const_iterator indexIter = m_indexData.find(region.LeftRefID);
    if ( indexIter == m_indexData.end()) return false;

    // load index data for region if not already cached
    if ( !IsDataLoaded(region.LeftRefID) ) {
        bool loadedOk = true;
        loadedOk &= SkipToReference(region.LeftRefID);
        loadedOk &= LoadReference(region.LeftRefID);
        if ( !loadedOk ) return false;
    }

    // localize index data for this reference (& sanity check that data actually exists)
    indexIter = m_indexData.find(region.LeftRefID);
    if ( indexIter == m_indexData.end()) return false;
    const vector<BamToolsIndexEntry>& referenceOffsets = (*indexIter).second.Offsets;
    if ( referenceOffsets.empty() ) return false;

    // -------------------------------------------------------
    // calculate nearest index to jump to
    
    // save first offset
    offset = (*referenceOffsets.begin()).StartOffset;
    
    // iterate over offsets entries on this reference
    vector<BamToolsIndexEntry>::const_iterator offsetIter = referenceOffsets.begin();
    vector<BamToolsIndexEntry>::const_iterator offsetEnd  = referenceOffsets.end();
    for ( ; offsetIter != offsetEnd; ++offsetIter ) {
	const BamToolsIndexEntry& entry = (*offsetIter);
        // break if alignment 'entry' overlaps region
        if ( entry.MaxEndPosition >= region.LeftPosition ) break;
	offset = (*offsetIter).StartOffset;
    }
  
    // set flag based on whether an index entry was found for this region
    *hasAlignmentsInRegion = ( offsetIter != offsetEnd );

    // if cache mode set to none, dump the data we just loaded
    if (m_parent->m_cacheMode == BamIndex::NoIndexCaching )
        ClearReferenceOffsets(region.LeftRefID);

    // return success
    return true; 
}

// returns whether reference has alignments or no
bool BamToolsIndex::BamToolsIndexPrivate::HasAlignments(const int& refId) const {

    BamToolsIndexData::const_iterator indexIter = m_indexData.find(refId);
    if ( indexIter == m_indexData.end()) return false;
    const BamToolsReferenceEntry& refEntry = (*indexIter).second;
    return refEntry.HasAlignments;
}

// return true if all index data is cached
bool BamToolsIndex::BamToolsIndexPrivate::HasFullDataCache(void) const {
    return m_hasFullDataCache;
}

// returns true if index cache has data for desired reference
bool BamToolsIndex::BamToolsIndexPrivate::IsDataLoaded(const int& refId) const {

    BamToolsIndexData::const_iterator indexIter = m_indexData.find(refId);
    if ( indexIter == m_indexData.end()) return false;
    const BamToolsReferenceEntry& refEntry = (*indexIter).second;

    if ( !refEntry.HasAlignments ) return true; // no data period

    // return whether offsets list contains data
    return !refEntry.Offsets.empty();
}

// attempts to use index to jump to region; returns success/fail
bool BamToolsIndex::BamToolsIndexPrivate::Jump(const BamRegion& region, bool* hasAlignmentsInRegion) {
  
    // clear flag
    *hasAlignmentsInRegion = false;
  
    // localize parent data
    if ( m_parent == 0 ) return false;
    BamReader* reader     = m_parent->m_reader;
    BgzfData*  mBGZF      = m_parent->m_BGZF;
    RefVector& references = m_parent->m_references;
  
    // check valid BamReader state
    if ( reader == 0 || mBGZF == 0 || !reader->IsOpen() ) {
        fprintf(stderr, "ERROR: Could not jump: invalid BamReader state.\n");
        return false;
    }
  
    // make sure left-bound position is valid
    if ( region.LeftPosition > references.at(region.LeftRefID).RefLength )
        return false;
  
    // calculate nearest offset to jump to
    int64_t offset;
    if ( !GetOffset(region, offset, hasAlignmentsInRegion) ) {
        fprintf(stderr, "ERROR: Could not jump - unable to calculate offset for specified region.\n");
        return false;
    }
    
    // return success/failure of seek
    return mBGZF->Seek(offset);    
}

// clears index data from all references except the first
void BamToolsIndex::BamToolsIndexPrivate::KeepOnlyFirstReferenceOffsets(void) {
    BamToolsIndexData::const_iterator indexBegin = m_indexData.begin();
    KeepOnlyReferenceOffsets( (*indexBegin).first );
}

// clears index data from all references except the one specified
void BamToolsIndex::BamToolsIndexPrivate::KeepOnlyReferenceOffsets(const int& refId) {
    BamToolsIndexData::iterator mapIter = m_indexData.begin();
    BamToolsIndexData::iterator mapEnd  = m_indexData.end();
    for ( ; mapIter != mapEnd; ++mapIter ) {
        const int entryRefId = (*mapIter).first;
        if ( entryRefId != refId )
            ClearReferenceOffsets(entryRefId);
    }
}

// load index data for all references, return true if loaded OK
bool BamToolsIndex::BamToolsIndexPrivate::LoadAllReferences(bool saveData) {

    // skip if data already loaded
    if ( m_hasFullDataCache ) return true;

    // read in number of references
    int32_t numReferences;
    if ( !LoadReferenceCount(numReferences) ) return false;
    //SetReferenceCount(numReferences);

    // iterate over reference entries
    bool loadedOk = true;
    for ( int i = 0; i < numReferences; ++i )
        loadedOk &= LoadReference(i, saveData);

    // set flag
    if ( loadedOk && saveData )
        m_hasFullDataCache = true;

    // return success/failure of load
    return loadedOk;
}

// load header data from index file, return true if loaded OK
bool BamToolsIndex::BamToolsIndexPrivate::LoadHeader(void) {

    // check magic number
    if ( !CheckMagicNumber() ) return false;

    // check BTI version
    if ( !CheckVersion() ) return false;

    // read in block size
    size_t elementsRead = fread(&m_blockSize, sizeof(m_blockSize), 1, m_parent->m_indexStream);
    if ( elementsRead != 1 ) return false;
    if ( m_isBigEndian ) SwapEndian_32(m_blockSize);

    // store offset of beginning of data
    m_dataBeginOffset = ftell64(m_parent->m_indexStream);

    // return success/failure of load
    return (elementsRead == 1);
}

// load a single index entry from file, return true if loaded OK
// @saveData - save data in memory if true, just read & discard if false
bool BamToolsIndex::BamToolsIndexPrivate::LoadIndexEntry(const int& refId, bool saveData) {
    
    // read in index entry data members
    size_t elementsRead = 0;
    BamToolsIndexEntry entry;
    elementsRead += fread(&entry.MaxEndPosition, sizeof(entry.MaxEndPosition), 1, m_parent->m_indexStream);
    elementsRead += fread(&entry.StartOffset,    sizeof(entry.StartOffset),    1, m_parent->m_indexStream);
    elementsRead += fread(&entry.StartPosition,  sizeof(entry.StartPosition),  1, m_parent->m_indexStream);
    if ( elementsRead != 3 ) {
        cerr << "Error reading index entry. Expected 3 elements, read in: " << elementsRead << endl;
        return false;
    }

    // swap endian-ness if necessary
    if ( m_isBigEndian ) {
        SwapEndian_32(entry.MaxEndPosition);
        SwapEndian_64(entry.StartOffset);
        SwapEndian_32(entry.StartPosition);
    }

    // save data
    if ( saveData )
        SaveOffsetEntry(refId, entry);

    // return success/failure of load
    return true;
}

// load a single reference from file, return true if loaded OK
// @saveData - save data in memory if true, just read & discard if false
bool BamToolsIndex::BamToolsIndexPrivate::LoadFirstReference(bool saveData) {
    BamToolsIndexData::const_iterator indexBegin = m_indexData.begin();
    return LoadReference( (*indexBegin).first, saveData );
}

// load a single reference from file, return true if loaded OK
// @saveData - save data in memory if true, just read & discard if false
bool BamToolsIndex::BamToolsIndexPrivate::LoadReference(const int& refId, bool saveData) {
  
    // read in number of offsets for this reference
    uint32_t numOffsets;
    size_t elementsRead = fread(&numOffsets, sizeof(numOffsets), 1, m_parent->m_indexStream);
    if ( elementsRead != 1 ) return false;
    if ( m_isBigEndian ) SwapEndian_32(numOffsets);
    
    // initialize offsets container for this reference
    SetOffsetCount(refId, (int)numOffsets);
    
    // iterate over offset entries
    for ( unsigned int j = 0; j < numOffsets; ++j )
        LoadIndexEntry(refId, saveData);

    // return success/failure of load
    return true;
}

// loads number of references, return true if loaded OK
bool BamToolsIndex::BamToolsIndexPrivate::LoadReferenceCount(int& numReferences) {

    size_t elementsRead = 0;

    // read reference count
    elementsRead += fread(&numReferences, sizeof(numReferences), 1, m_parent->m_indexStream);
    if ( m_isBigEndian ) SwapEndian_32(numReferences);

    // return success/failure of load
    return ( elementsRead == 1 );
}

// saves an index offset entry in memory
void BamToolsIndex::BamToolsIndexPrivate::SaveOffsetEntry(const int& refId, const BamToolsIndexEntry& entry) {
    BamToolsReferenceEntry& refEntry = m_indexData[refId];
    refEntry.HasAlignments = true;
    refEntry.Offsets.push_back(entry);
}

// pre-allocates size for offset vector
void BamToolsIndex::BamToolsIndexPrivate::SetOffsetCount(const int& refId, const int& offsetCount) {
    BamToolsReferenceEntry& refEntry = m_indexData[refId];
    refEntry.Offsets.reserve(offsetCount);
    refEntry.HasAlignments = ( offsetCount > 0);
}

// initializes index data structure to hold @count references
void BamToolsIndex::BamToolsIndexPrivate::SetReferenceCount(const int& count) {
    for ( int i = 0; i < count; ++i )
        m_indexData[i].HasAlignments = false;
}

// position file pointer to first reference begin, return true if skipped OK
bool BamToolsIndex::BamToolsIndexPrivate::SkipToFirstReference(void) {
    BamToolsIndexData::const_iterator indexBegin = m_indexData.begin();
    return SkipToReference( (*indexBegin).first );
}

// position file pointer to desired reference begin, return true if skipped OK
bool BamToolsIndex::BamToolsIndexPrivate::SkipToReference(const int& refId) {

    // attempt rewind
    if ( !m_parent->Rewind() ) return false;

    // read in number of references
    int32_t numReferences;
    size_t elementsRead = fread(&numReferences, sizeof(numReferences), 1, m_parent->m_indexStream);
    if ( elementsRead != 1 ) return false;
    if ( m_isBigEndian ) SwapEndian_32(numReferences);

    // iterate over reference entries
    bool skippedOk = true;
    int currentRefId = 0;
    while (currentRefId != refId) {
        skippedOk &= LoadReference(currentRefId, false);
        ++currentRefId;
    }

    // return success/failure of skip
    return skippedOk;
}

// write header to new index file
bool BamToolsIndex::BamToolsIndexPrivate::WriteHeader(void) {

    size_t elementsWritten = 0;

    // write BTI index format 'magic number'
    elementsWritten += fwrite("BTI\1", 1, 4, m_parent->m_indexStream);

    // write BTI index format version
    int32_t currentVersion = (int32_t)m_outputVersion;
    if ( m_isBigEndian ) SwapEndian_32(currentVersion);
    elementsWritten += fwrite(&currentVersion, sizeof(currentVersion), 1, m_parent->m_indexStream);

    // write block size
    int32_t blockSize = m_blockSize;
    if ( m_isBigEndian ) SwapEndian_32(blockSize);
    elementsWritten += fwrite(&blockSize, sizeof(blockSize), 1, m_parent->m_indexStream);

    // store offset of beginning of data
    m_dataBeginOffset = ftell64(m_parent->m_indexStream);

    // return success/failure of write
    return ( elementsWritten == 6 );
}

// write index data for all references to new index file
bool BamToolsIndex::BamToolsIndexPrivate::WriteAllReferences(void) {

    size_t elementsWritten = 0;

    // write number of references
    int32_t numReferences = (int32_t)m_indexData.size();
    if ( m_isBigEndian ) SwapEndian_32(numReferences);
    elementsWritten += fwrite(&numReferences, sizeof(numReferences), 1, m_parent->m_indexStream);

    // iterate through references in index
    bool refOk = true;
    BamToolsIndexData::const_iterator refIter = m_indexData.begin();
    BamToolsIndexData::const_iterator refEnd  = m_indexData.end();
    for ( ; refIter != refEnd; ++refIter )
        refOk &= WriteReferenceEntry( (*refIter).second );

    return ( (elementsWritten == 1) && refOk );
}

// write current reference index data to new index file
bool BamToolsIndex::BamToolsIndexPrivate::WriteReferenceEntry(const BamToolsReferenceEntry& refEntry) {

    size_t elementsWritten = 0;

    // write number of offsets listed for this reference
    uint32_t numOffsets = refEntry.Offsets.size();
    if ( m_isBigEndian ) SwapEndian_32(numOffsets);
    elementsWritten += fwrite(&numOffsets, sizeof(numOffsets), 1, m_parent->m_indexStream);

    // iterate over offset entries
    bool entriesOk = true;
    vector<BamToolsIndexEntry>::const_iterator offsetIter = refEntry.Offsets.begin();
    vector<BamToolsIndexEntry>::const_iterator offsetEnd  = refEntry.Offsets.end();
    for ( ; offsetIter != offsetEnd; ++offsetIter )
        entriesOk &= WriteIndexEntry( (*offsetIter) );

    return ( (elementsWritten == 1) && entriesOk );
}

// write current index offset entry to new index file
bool BamToolsIndex::BamToolsIndexPrivate::WriteIndexEntry(const BamToolsIndexEntry& entry) {

    // copy entry data
    int32_t maxEndPosition = entry.MaxEndPosition;
    int64_t startOffset    = entry.StartOffset;
    int32_t startPosition  = entry.StartPosition;

    // swap endian-ness if necessary
    if ( m_isBigEndian ) {
        SwapEndian_32(maxEndPosition);
        SwapEndian_64(startOffset);
        SwapEndian_32(startPosition);
    }

    // write the reference index entry
    size_t elementsWritten = 0;
    elementsWritten += fwrite(&maxEndPosition, sizeof(maxEndPosition), 1, m_parent->m_indexStream);
    elementsWritten += fwrite(&startOffset,    sizeof(startOffset),    1, m_parent->m_indexStream);
    elementsWritten += fwrite(&startPosition,  sizeof(startPosition),  1, m_parent->m_indexStream);
    return ( elementsWritten == 3 );
}

// ---------------------------------------------------
// BamToolsIndex implementation

BamToolsIndex::BamToolsIndex(BgzfData* bgzf, BamReader* reader)
    : BamIndex(bgzf, reader)
{ 
    d = new BamToolsIndexPrivate(this);
}    

BamToolsIndex::~BamToolsIndex(void) { 
    delete d;
    d = 0;
}

// BamIndex interface implementation
bool BamToolsIndex::Build(void) { return d->Build(); }
void BamToolsIndex::ClearAllData(void) { d->ClearAllData(); }
const off_t BamToolsIndex::DataBeginOffset(void) const { return d->DataBeginOffset(); }
bool BamToolsIndex::HasAlignments(const int& referenceID) const { return d->HasAlignments(referenceID); }
bool BamToolsIndex::HasFullDataCache(void) const { return d->HasFullDataCache(); }
bool BamToolsIndex::Jump(const BamRegion& region, bool* hasAlignmentsInRegion) { return d->Jump(region, hasAlignmentsInRegion); }
void BamToolsIndex::KeepOnlyFirstReferenceOffsets(void) { d->KeepOnlyFirstReferenceOffsets(); }
bool BamToolsIndex::LoadAllReferences(bool saveData) { return d->LoadAllReferences(saveData); }
bool BamToolsIndex::LoadFirstReference(bool saveData) { return d->LoadFirstReference(saveData); }
bool BamToolsIndex::LoadHeader(void) { return d->LoadHeader(); }
bool BamToolsIndex::SkipToFirstReference(void) { return d->SkipToFirstReference(); }
bool BamToolsIndex::WriteAllReferences(void) { return d->WriteAllReferences(); }
bool BamToolsIndex::WriteHeader(void) { return d->WriteHeader(); }

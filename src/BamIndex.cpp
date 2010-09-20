// ***************************************************************************
// BamIndex.cpp (c) 2009 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ---------------------------------------------------------------------------
// Last modified: 18 September 2010 (DB)
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

BamIndex::BamIndex(BamTools::BgzfData* bgzf, BamTools::BamReader* reader, bool isBigEndian) 
    : m_BGZF(bgzf)
    , m_reader(reader)
    , m_isBigEndian(isBigEndian)
{ 
    if ( m_reader && m_reader->IsOpen() ) 
        m_references = m_reader->GetReferenceData();
}

bool BamIndex::HasAlignments(const int& referenceID) {
    
    // return false if invalid ID
    if ( (referenceID < 0) || (referenceID >= (int)m_references.size()) ) 
        return false;
    
    // else return status of reference (has alignments?)
    else
        return m_references.at(referenceID).RefHasAlignments;
}

// #########################################################################################
// #########################################################################################

// -------------------------------
// BamStandardIndex structs & typedefs 
 
namespace BamTools { 

// BAM index constants
const int MAX_BIN           = 37450;    // =(8^6-1)/7+1
const int BAM_LIDX_SHIFT    = 14;  
  
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
    
    // constructor
    ReferenceIndex(const BamBinMap& binMap = BamBinMap(),
                   const LinearOffsetVector& offsets = LinearOffsetVector())
        : Bins(binMap)
        , Offsets(offsets)
    { }
};

typedef vector<ReferenceIndex> BamStandardIndexData;

} // namespace BamTools
 
// -------------------------------
// BamStandardIndexPrivate implementation
  
struct BamStandardIndex::BamStandardIndexPrivate { 
  
    // -------------------------
    // data members
    
    BamStandardIndexData m_indexData;
    BamStandardIndex*    m_parent;
    
    // -------------------------
    // ctor & dtor
    
    BamStandardIndexPrivate(BamStandardIndex* parent) : m_parent(parent) { }
    ~BamStandardIndexPrivate(void) { m_indexData.clear(); }
    
    // -------------------------
    // 'public' methods
    
    // creates index data (in-memory) from current reader data
    bool Build(void);
    // attempts to use index to jump to region; returns success/fail
    bool Jump(const BamTools::BamRegion& region, bool* hasAlignmentsInRegion);
      // loads existing data from file into memory
    bool Load(const std::string& filename);
    // writes in-memory index data out to file 
    // N.B. - (this is the original BAM filename, method will modify it to use applicable extension)
    bool Write(const std::string& bamFilename);
    
    // -------------------------
    // internal methods
    
    // calculate bins that overlap region
    int BinsFromRegion(const BamRegion& region, const bool isRightBoundSpecified, uint16_t bins[BamTools::MAX_BIN]);
    // calculates offset(s) for a given region
    bool GetOffsets(const BamRegion& region, const bool isRightBoundSpecified, vector<int64_t>& offsets);
    // saves BAM bin entry for index
    void InsertBinEntry(BamBinMap& binMap, const uint32_t& saveBin, const uint64_t& saveOffset, const uint64_t& lastOffset);
    // saves linear offset entry for index
    void InsertLinearOffset(LinearOffsetVector& offsets, const BamAlignment& bAlignment, const uint64_t& lastOffset);
    // simplifies index by merging 'chunks'
    void MergeChunks(void);
};

// calculate bins that overlap region
int BamStandardIndex::BamStandardIndexPrivate::BinsFromRegion(const BamRegion& region, const bool isRightBoundSpecified, uint16_t bins[MAX_BIN]) {
  
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

bool BamStandardIndex::BamStandardIndexPrivate::Build(void) { 
  
    // localize parent data
    if ( m_parent == 0 ) return false;
    BamReader* reader = m_parent->m_reader;
    BgzfData*  mBGZF  = m_parent->m_BGZF;
    RefVector& references = m_parent->m_references;
  
    // be sure reader & BGZF file are valid & open for reading
    if ( reader == 0 || mBGZF == 0 || !mBGZF->IsOpen ) 
        return false;

    // move file pointer to beginning of alignments
    reader->Rewind();

    // get reference count, reserve index space
    int numReferences = (int)references.size();
    for ( int i = 0; i < numReferences; ++i ) 
        m_indexData.push_back(ReferenceIndex());
    
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
            fprintf(stderr, "Alignment %s : %d > %d on reference (id = %d)", bAlignment.Name.c_str(), lastCoordinate, bAlignment.Position, bAlignment.RefID);
            exit(1);
        }

        // if valid reference && BAM bin spans some minimum cutoff (smaller bin ids span larger regions)
        if ( (bAlignment.RefID >= 0) && (bAlignment.Bin < 4681) ) {

            // save linear offset entry (matched to BAM entry refID)
            ReferenceIndex& refIndex = m_indexData.at(bAlignment.RefID);
            LinearOffsetVector& offsets = refIndex.Offsets;
            InsertLinearOffset(offsets, bAlignment, lastOffset);
        }

        // if current BamAlignment bin != lastBin, "then possibly write the binning index"
        if ( bAlignment.Bin != lastBin ) {

            // if not first time through
            if ( saveBin != defaultValue ) {

                // save Bam bin entry
                ReferenceIndex& refIndex = m_indexData.at(saveRefID);
                BamBinMap& binMap = refIndex.Bins;
                InsertBinEntry(binMap, saveBin, saveOffset, lastOffset);
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
        ReferenceIndex& refIndex = m_indexData.at(saveRefID);
        BamBinMap& binMap = refIndex.Bins;
        InsertBinEntry(binMap, saveBin, saveOffset, lastOffset);
    }

    // simplify index by merging chunks
    MergeChunks();
    
    // iterate through references in index
    // store whether reference has data &
    // sort offsets in linear offset vector
    BamStandardIndexData::iterator indexIter = m_indexData.begin();
    BamStandardIndexData::iterator indexEnd  = m_indexData.end();
    for ( int i = 0; indexIter != indexEnd; ++indexIter, ++i ) {

        // get reference index data
        ReferenceIndex& refIndex = (*indexIter);
        BamBinMap& binMap = refIndex.Bins;
        LinearOffsetVector& offsets = refIndex.Offsets;

        // store whether reference has alignments or no
        references[i].RefHasAlignments = ( binMap.size() > 0 );

        // sort linear offsets
        sort(offsets.begin(), offsets.end());
    }

    // rewind file pointer to beginning of alignments, return success/fail
    return reader->Rewind();
}

bool BamStandardIndex::BamStandardIndexPrivate::GetOffsets(const BamRegion& region, const bool isRightBoundSpecified, vector<int64_t>& offsets) { 
  
    // calculate which bins overlap this region
    uint16_t* bins = (uint16_t*)calloc(MAX_BIN, 2);
    int numBins = BinsFromRegion(region, isRightBoundSpecified, bins);

    // get bins for this reference
    const ReferenceIndex& refIndex = m_indexData.at(region.LeftRefID);
    const BamBinMap& binMap        = refIndex.Bins;

    // get minimum offset to consider
    const LinearOffsetVector& linearOffsets = refIndex.Offsets;
    uint64_t minOffset = ( (unsigned int)(region.LeftPosition>>BAM_LIDX_SHIFT) >= linearOffsets.size() ) ? 0 : linearOffsets.at(region.LeftPosition>>BAM_LIDX_SHIFT);

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
    
    // return whether any offsets were found
    return ( offsets.size() != 0 );
}

// saves BAM bin entry for index
void BamStandardIndex::BamStandardIndexPrivate::InsertBinEntry(BamBinMap&      binMap,
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
void BamStandardIndex::BamStandardIndexPrivate::InsertLinearOffset(LinearOffsetVector& offsets,
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

bool BamStandardIndex::BamStandardIndexPrivate::Jump(const BamRegion& region, bool* hasAlignmentsInRegion) {
  
    // localize parent data
    if ( m_parent == 0 ) return false;
    BamReader* reader = m_parent->m_reader;
    BgzfData*  mBGZF  = m_parent->m_BGZF;
    RefVector& references = m_parent->m_references;
  
    // be sure reader & BGZF file are valid & open for reading
    if ( reader == 0 || mBGZF == 0 || !mBGZF->IsOpen ) 
        return false;
    
    // see if left-bound reference of region has alignments
    if ( !m_parent->HasAlignments(region.LeftRefID) ) return false; 
    
    // make sure left-bound position is valid
    if ( region.LeftPosition > references.at(region.LeftRefID).RefLength ) return false;
        
    // calculate offsets for this region
    // if failed, print message, set flag, and return failure
    vector<int64_t> offsets;
    if ( !GetOffsets(region, region.isRightBoundSpecified(), offsets) ) {
        fprintf(stderr, "ERROR: Could not jump: unable to calculate offset for specified region.\n");
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
        if ( (bAlignment.RefID == region.LeftRefID && bAlignment.Position + bAlignment.Length > region.LeftPosition) || (bAlignment.RefID > region.LeftRefID) ) {
            if ( o != offsets.begin() ) --o;
            return mBGZF->Seek(*o);
        }
    }
    
    // if error in jumping, print message & set flag
    if ( !result ) {
        fprintf(stderr, "ERROR: Could not jump: unable to calculate offset for specified region.\n");
        *hasAlignmentsInRegion = false;
    }
    
    // return success/failure
    return result;
}

bool BamStandardIndex::BamStandardIndexPrivate::Load(const string& filename)  { 
    
    // localize parent data
    if ( m_parent == 0 ) return false;
    const bool isBigEndian = m_parent->m_isBigEndian;
    RefVector& references  = m_parent->m_references;
  
    // open index file, abort on error
    FILE* indexStream = fopen(filename.c_str(), "rb");
    if( !indexStream ) {
        fprintf(stderr, "ERROR: Unable to open the BAM index file %s for reading.\n", filename.c_str());
        return false;
    }

    // set placeholder to receive input byte count (suppresses compiler warnings)
    size_t elementsRead = 0;
        
    // see if index is valid BAM index
    char magic[4];
    elementsRead = fread(magic, 1, 4, indexStream);
    if ( strncmp(magic, "BAI\1", 4) ) {
        fprintf(stderr, "Problem with index file - invalid format.\n");
        fclose(indexStream);
        return false;
    }

    // get number of reference sequences
    uint32_t numRefSeqs;
    elementsRead = fread(&numRefSeqs, 4, 1, indexStream);
    if ( isBigEndian ) SwapEndian_32(numRefSeqs);
    
    // intialize space for BamStandardIndexData data structure
    m_indexData.reserve(numRefSeqs);

    // iterate over reference sequences
    for ( unsigned int i = 0; i < numRefSeqs; ++i ) {

        // get number of bins for this reference sequence
        int32_t numBins;
        elementsRead = fread(&numBins, 4, 1, indexStream);
        if ( isBigEndian ) SwapEndian_32(numBins);

        if ( numBins > 0 ) {
            RefData& refEntry = references[i];
            refEntry.RefHasAlignments = true;
        }

        // intialize BinVector
        BamBinMap binMap;

        // iterate over bins for that reference sequence
        for ( int j = 0; j < numBins; ++j ) {

            // get binID
            uint32_t binID;
            elementsRead = fread(&binID, 4, 1, indexStream);

            // get number of regionChunks in this bin
            uint32_t numChunks;
            elementsRead = fread(&numChunks, 4, 1, indexStream);

            if ( isBigEndian ) { 
              SwapEndian_32(binID);
              SwapEndian_32(numChunks);
            }
            
            // intialize ChunkVector
            ChunkVector regionChunks;
            regionChunks.reserve(numChunks);

            // iterate over regionChunks in this bin
            for ( unsigned int k = 0; k < numChunks; ++k ) {

                // get chunk boundaries (left, right)
                uint64_t left;
                uint64_t right;
                elementsRead = fread(&left,  8, 1, indexStream);
                elementsRead = fread(&right, 8, 1, indexStream);

                if ( isBigEndian ) {
                    SwapEndian_64(left);
                    SwapEndian_64(right);
                }
                
                // save ChunkPair
                regionChunks.push_back( Chunk(left, right) );
            }

            // sort chunks for this bin
            sort( regionChunks.begin(), regionChunks.end(), ChunkLessThan );

            // save binID, chunkVector for this bin
            binMap.insert( pair<uint32_t, ChunkVector>(binID, regionChunks) );
        }

        // -----------------------------------------------------
        // load linear index for this reference sequence

        // get number of linear offsets
        int32_t numLinearOffsets;
        elementsRead = fread(&numLinearOffsets, 4, 1, indexStream);
        if ( isBigEndian ) SwapEndian_32(numLinearOffsets);

        // intialize LinearOffsetVector
        LinearOffsetVector offsets;
        offsets.reserve(numLinearOffsets);

        // iterate over linear offsets for this reference sequeence
        uint64_t linearOffset;
        for ( int j = 0; j < numLinearOffsets; ++j ) {
            // read a linear offset & store
            elementsRead = fread(&linearOffset, 8, 1, indexStream);
            if ( isBigEndian ) SwapEndian_64(linearOffset);
            offsets.push_back(linearOffset);
        }

        // sort linear offsets
        sort( offsets.begin(), offsets.end() );

        // store index data for that reference sequence
        m_indexData.push_back( ReferenceIndex(binMap, offsets) );
    }

    // close index file (.bai) and return
    fclose(indexStream);
    return true;
}

// merges 'alignment chunks' in BAM bin (used for index building)
void BamStandardIndex::BamStandardIndexPrivate::MergeChunks(void) {

    // iterate over reference enties
    BamStandardIndexData::iterator indexIter = m_indexData.begin();
    BamStandardIndexData::iterator indexEnd  = m_indexData.end();
    for ( ; indexIter != indexEnd; ++indexIter ) {

        // get BAM bin map for this reference
        ReferenceIndex& refIndex = (*indexIter);
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

                // if currentChunk.Stop(shifted) == iterator Chunk.Start(shifted)
                if ( currentChunk.Stop>>16 == iteratorChunk.Start>>16 ) {

                    // set currentChunk.Stop to iteratorChunk.Stop
                    currentChunk.Stop = iteratorChunk.Stop;
                }

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

// writes in-memory index data out to file 
// N.B. - (this is the original BAM filename, method will modify it to use applicable extension)
bool BamStandardIndex::BamStandardIndexPrivate::Write(const std::string& bamFilename) { 

    // localize parent data
    if ( m_parent == 0 ) return false;
    const bool isBigEndian = m_parent->m_isBigEndian;
  
    string indexFilename = bamFilename + ".bai";
    FILE* indexStream = fopen(indexFilename.c_str(), "wb");
    if ( indexStream == 0 ) {
        fprintf(stderr, "ERROR: Could not open file to save index.\n");
        return false;
    }

    // write BAM index header
    fwrite("BAI\1", 1, 4, indexStream);

    // write number of reference sequences
    int32_t numReferenceSeqs = m_indexData.size();
    if ( isBigEndian ) SwapEndian_32(numReferenceSeqs);
    fwrite(&numReferenceSeqs, 4, 1, indexStream);

    // iterate over reference sequences
    BamStandardIndexData::const_iterator indexIter = m_indexData.begin();
    BamStandardIndexData::const_iterator indexEnd  = m_indexData.end();
    for ( ; indexIter != indexEnd; ++ indexIter ) {

        // get reference index data
        const ReferenceIndex& refIndex = (*indexIter);
        const BamBinMap& binMap = refIndex.Bins;
        const LinearOffsetVector& offsets = refIndex.Offsets;

        // write number of bins
        int32_t binCount = binMap.size();
        if ( isBigEndian ) SwapEndian_32(binCount);
        fwrite(&binCount, 4, 1, indexStream);

        // iterate over bins
        BamBinMap::const_iterator binIter = binMap.begin();
        BamBinMap::const_iterator binEnd  = binMap.end();
        for ( ; binIter != binEnd; ++binIter ) {

            // get bin data (key and chunk vector)
            uint32_t binKey = (*binIter).first;
            const ChunkVector& binChunks = (*binIter).second;

            // save BAM bin key
            if ( isBigEndian ) SwapEndian_32(binKey);
            fwrite(&binKey, 4, 1, indexStream);

            // save chunk count
            int32_t chunkCount = binChunks.size();
            if ( isBigEndian ) SwapEndian_32(chunkCount); 
            fwrite(&chunkCount, 4, 1, indexStream);

            // iterate over chunks
            ChunkVector::const_iterator chunkIter = binChunks.begin();
            ChunkVector::const_iterator chunkEnd  = binChunks.end();
            for ( ; chunkIter != chunkEnd; ++chunkIter ) {

                // get current chunk data
                const Chunk& chunk    = (*chunkIter);
                uint64_t start = chunk.Start;
                uint64_t stop  = chunk.Stop;

                if ( isBigEndian ) {
                    SwapEndian_64(start);
                    SwapEndian_64(stop);
                }
                
                // save chunk offsets
                fwrite(&start, 8, 1, indexStream);
                fwrite(&stop,  8, 1, indexStream);
            }
        }

        // write linear offsets size
        int32_t offsetSize = offsets.size();
        if ( isBigEndian ) SwapEndian_32(offsetSize);
        fwrite(&offsetSize, 4, 1, indexStream);

        // iterate over linear offsets
        LinearOffsetVector::const_iterator offsetIter = offsets.begin();
        LinearOffsetVector::const_iterator offsetEnd  = offsets.end();
        for ( ; offsetIter != offsetEnd; ++offsetIter ) {

            // write linear offset value
            uint64_t linearOffset = (*offsetIter);
            if ( isBigEndian ) SwapEndian_64(linearOffset);
            fwrite(&linearOffset, 8, 1, indexStream);
        }
    }

    // flush buffer, close file, and return success
    fflush(indexStream);
    fclose(indexStream);
    return true;
}
 
// ---------------------------------------------------------------
// BamStandardIndex implementation
 
BamStandardIndex::BamStandardIndex(BgzfData* bgzf, BamReader* reader, bool isBigEndian)
    : BamIndex(bgzf, reader, isBigEndian)
{
    d = new BamStandardIndexPrivate(this);
}    

BamStandardIndex::~BamStandardIndex(void) {
    delete d;
    d = 0;
}

// creates index data (in-memory) from current reader data
bool BamStandardIndex::Build(void) { return d->Build(); }

// attempts to use index to jump to region; returns success/fail
bool BamStandardIndex::Jump(const BamRegion& region, bool* hasAlignmentsInRegion) { return d->Jump(region, hasAlignmentsInRegion); }

// loads existing data from file into memory
bool BamStandardIndex::Load(const string& filename) { return d->Load(filename); }

// writes in-memory index data out to file 
// N.B. - (this is the original BAM filename, method will modify it to use applicable extension)
bool BamStandardIndex::Write(const string& bamFilename) { return d->Write(bamFilename); }

// #########################################################################################
// #########################################################################################

// ---------------------------------------------------
// BamToolsIndex structs & typedefs

namespace BamTools {
  
// individual index entry
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

// the actual index data structure
typedef map<int, vector<BamToolsIndexEntry> > BamToolsIndexData;
  
} // namespace BamTools

// ---------------------------------------------------
// BamToolsIndexPrivate implementation

struct BamToolsIndex::BamToolsIndexPrivate {
    
    // keep a list of any supported versions here 
    // (might be useful later to handle any 'legacy' versions if the format changes)
    // listed for example like: BTI_1_1 = 2, BTI_1_2 = 3, BTI_2_0 = 4, and so on
    //
    // so a change introduced in (hypothetical) BTI_1_2 would be handled from then on by: 
    //
    // if ( indexVersion >= BTI_1_2 ) 
    //   do something new 
    // else 
    //   do the old thing
    enum Version { BTI_1_0 = 1 
                 , BTI_1_1
                 };  
  
    // -------------------------
    // data members
    BamToolsIndexData m_indexData;
    BamToolsIndex*    m_parent;
    int32_t           m_blockSize;
    Version           m_version;
    
    // -------------------------
    // ctor & dtor
    
    BamToolsIndexPrivate(BamToolsIndex* parent) 
        : m_parent(parent)
        , m_blockSize(1000)
        , m_version(BTI_1_1) // latest version - used for writing new index files
    { }
    
    ~BamToolsIndexPrivate(void) { }
    
    // -------------------------
    // 'public' methods
    
    // creates index data (in-memory) from current reader data
    bool Build(void);
    // returns supported file extension
    const std::string Extension(void) const { return std::string(".bti"); }
    // attempts to use index to jump to region; returns success/fail
    bool Jump(const BamTools::BamRegion& region, bool* hasAlignmentsInRegion);
      // loads existing data from file into memory
    bool Load(const std::string& filename);
    // writes in-memory index data out to file 
    // N.B. - (this is the original BAM filename, method will modify it to use applicable extension)
    bool Write(const std::string& bamFilename);
    
    // -------------------------
    // internal methods
    
    // calculates offset for a given region
    bool GetOffset(const BamRegion& region, int64_t& offset, bool* hasAlignmentsInRegion);
};

bool BamToolsIndex::BamToolsIndexPrivate::Build(void) { 
  
    // localize parent data
    if ( m_parent == 0 ) return false;
    BamReader* reader = m_parent->m_reader;
    BgzfData*  mBGZF  = m_parent->m_BGZF;
    RefVector& references = m_parent->m_references;
  
    // be sure reader & BGZF file are valid & open for reading
    if ( reader == 0 || mBGZF == 0 || !mBGZF->IsOpen ) 
        return false;

    // move file pointer to beginning of alignments
    // quit if failed
    if ( !reader->Rewind() ) 
        return false;
    
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
          
            // make sure reference flag is set
            references[blockRefId].RefHasAlignments = true;
          
            // store previous data
            m_indexData[blockRefId].push_back( BamToolsIndexEntry(blockMaxEndPosition, blockStartOffset, blockStartPosition) );
            
            // intialize new block
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
          
            // make sure reference flag is set
            references[blockRefId].RefHasAlignments = true;
          
            m_indexData[blockRefId].push_back( BamToolsIndexEntry(blockMaxEndPosition, blockStartOffset, blockStartPosition) );
            blockStartOffset  = mBGZF->Tell();
            currentBlockCount = 0;
        }
        
        // not the best name, but for the next iteration, this value will be the offset of the *current* alignment
        // necessary because we won't know if this next alignment is on a new reference until we actually read it
        currentAlignmentOffset = mBGZF->Tell();  
    }
    
    // store final block
    m_indexData[blockRefId].push_back( BamToolsIndexEntry(blockMaxEndPosition, blockStartOffset, blockStartPosition) );
    
    // attempt to rewind back to begnning of alignments
    // return success/fail
    return reader->Rewind();
}

// N.B. - ignores isRightBoundSpecified
bool BamToolsIndex::BamToolsIndexPrivate::GetOffset(const BamRegion& region, int64_t& offset, bool* hasAlignmentsInRegion) { 
  
    // return false if leftBound refID is not found in index data
    if ( m_indexData.find(region.LeftRefID) == m_indexData.end() ) return false;
    
    const vector<BamToolsIndexEntry> referenceOffsets = m_indexData[region.LeftRefID];
    if ( referenceOffsets.empty() ) return false;
    
    // -------------------------------------------------------
    // calculate nearest index to jump to
    
    // save first offset
    offset = (*referenceOffsets.begin()).StartOffset;
    
    // iterate over offsets entries on this reference
    vector<BamToolsIndexEntry>::const_iterator indexIter = referenceOffsets.begin();
    vector<BamToolsIndexEntry>::const_iterator indexEnd  = referenceOffsets.end();
    for ( ; indexIter != indexEnd; ++indexIter ) {
        const BamToolsIndexEntry& entry = (*indexIter);
        // break if alignment 'entry' overlaps region
        if ( entry.MaxEndPosition >= region.LeftPosition ) break;
        offset = (*indexIter).StartOffset;
    }
  
    // set flag based on whether an index entry was found for this region
    *hasAlignmentsInRegion = ( indexIter != indexEnd );
    return true; 
}

bool BamToolsIndex::BamToolsIndexPrivate::Jump(const BamRegion& region, bool* hasAlignmentsInRegion) {
  
    // clear flag
    *hasAlignmentsInRegion = false;
  
    // localize parent data
    if ( m_parent == 0 ) return false;
    BamReader* reader = m_parent->m_reader;
    BgzfData*  mBGZF  = m_parent->m_BGZF;
    RefVector& references = m_parent->m_references;
  
    // check valid BamReader state
    if ( reader == 0 || mBGZF == 0 || !reader->IsOpen() ) {
        fprintf(stderr, "ERROR: Could not jump: invalid BamReader state.\n");
        return false;
    }
  
    // see if left-bound reference of region has alignments
    if ( !m_parent->HasAlignments(region.LeftRefID) ) return false; 
    
    // make sure left-bound position is valid
    if ( region.LeftPosition > references.at(region.LeftRefID).RefLength ) return false;
  
    // calculate nearest offset to jump to
    int64_t offset;
    if ( !GetOffset(region, offset, hasAlignmentsInRegion) ) {
        fprintf(stderr, "ERROR: Could not jump - unable to calculate offset for specified region.\n");
        return false;
    }
  
    // attempt seek in file, return success/fail
    return mBGZF->Seek(offset);    
}

bool BamToolsIndex::BamToolsIndexPrivate::Load(const string& filename) { 
  
    // localize parent data
    if ( m_parent == 0 ) return false;
    const bool isBigEndian = m_parent->m_isBigEndian;
    RefVector& references  = m_parent->m_references;
  
    // open index file, abort on error
    FILE* indexStream = fopen(filename.c_str(), "rb");
    if( !indexStream ) {
        fprintf(stderr, "ERROR: Unable to open the BAM index file %s for reading.\n", filename.c_str());
        return false;
    }

    // set placeholder to receive input byte count (suppresses compiler warnings)
    size_t elementsRead = 0;
        
    // see if index is valid BAM index
    char magic[4];
    elementsRead = fread(magic, 1, 4, indexStream);
    if ( strncmp(magic, "BTI\1", 4) ) {
        fprintf(stderr, "Problem with index file - invalid format.\n");
        fclose(indexStream);
        return false;
    }

    // read in version
    int32_t indexVersion;
    elementsRead = fread(&indexVersion, sizeof(indexVersion), 1, indexStream);
    if ( isBigEndian ) SwapEndian_32(indexVersion);
    if ( indexVersion <= 0 || indexVersion > m_version ) {
        fprintf(stderr, "Problem with index file - unsupported version.\n");
        fclose(indexStream);
        return false;
    }

    if ( (Version)indexVersion < BTI_1_1 ) {
        fprintf(stderr, "\nProblem with index file - this version of the index contains a bug related to accessing data near reference ends.\n");
        fprintf(stderr, "\nPlease run \'bamtools index -bti -in yourData.bam\' to generate an up-to-date BamToolsIndex.\n\n");
        fclose(indexStream);
        exit(1);
    }

    // read in block size
    elementsRead = fread(&m_blockSize, sizeof(m_blockSize), 1, indexStream);
    if ( isBigEndian ) SwapEndian_32(m_blockSize);
    
    // read in number of references
    int32_t numReferences;
    elementsRead = fread(&numReferences, sizeof(numReferences), 1, indexStream);
    if ( isBigEndian ) SwapEndian_32(numReferences);
    
    // iterate over reference entries
    for ( int i = 0; i < numReferences; ++i ) {
      
        // read in number of offsets for this reference
        uint32_t numOffsets;
        elementsRead = fread(&numOffsets, sizeof(numOffsets), 1, indexStream);
        if ( isBigEndian ) SwapEndian_32(numOffsets);
        
        // initialize offsets container for this reference
        vector<BamToolsIndexEntry> offsets;
        offsets.reserve(numOffsets);
        
        // iterate over offset entries
        for ( unsigned int j = 0; j < numOffsets; ++j ) {
          
            // copy entry data
            int32_t maxEndPosition;
            int64_t startOffset;
            int32_t startPosition;
          
            // read in data
            elementsRead = fread(&maxEndPosition, sizeof(maxEndPosition), 1, indexStream);
            elementsRead = fread(&startOffset,    sizeof(startOffset),    1, indexStream);
            elementsRead = fread(&startPosition,  sizeof(startPosition),  1, indexStream);
            
            // swap endian-ness if necessary
            if ( isBigEndian ) {
                SwapEndian_32(maxEndPosition);
                SwapEndian_64(startOffset);
                SwapEndian_32(startPosition);
            }
          
            // save current index entry
            offsets.push_back( BamToolsIndexEntry(maxEndPosition, startOffset, startPosition) );
        }
        
        // save refID, offsetVector entry into data structure
        m_indexData.insert( make_pair(i, offsets) );
        
        // set ref.HasAlignments flag
        references[i].RefHasAlignments = (numOffsets != 0);
    }

    // close index file and return
    fclose(indexStream);
    return true;
}

// writes in-memory index data out to file 
// N.B. - (this is the original BAM filename, method will modify it to use applicable extension)
bool BamToolsIndex::BamToolsIndexPrivate::Write(const std::string& bamFilename) { 
    
    // localize parent data
    if ( m_parent == 0 ) return false;
    const bool isBigEndian = m_parent->m_isBigEndian;
  
    // open index file for writing
    string indexFilename = bamFilename + ".bti";
    FILE* indexStream = fopen(indexFilename.c_str(), "wb");
    if ( indexStream == 0 ) {
        fprintf(stderr, "ERROR: Could not open file to save index.\n");
        return false;
    }

    // write BTI index format 'magic number'
    fwrite("BTI\1", 1, 4, indexStream);

    // write BTI index format version
    int32_t currentVersion = (int32_t)m_version;
    if ( isBigEndian ) SwapEndian_32(currentVersion);
    fwrite(&currentVersion, sizeof(int32_t), 1, indexStream);
    
    // write block size
    int32_t blockSize = m_blockSize;
    if ( isBigEndian ) SwapEndian_32(blockSize);
    fwrite(&blockSize, sizeof(blockSize), 1, indexStream);
    
    // write number of references
    int32_t numReferences = (int32_t)m_indexData.size();
    if ( isBigEndian ) SwapEndian_32(numReferences);
    fwrite(&numReferences, sizeof(numReferences), 1, indexStream);
    
    // iterate through references in index 
    BamToolsIndexData::const_iterator refIter = m_indexData.begin();
    BamToolsIndexData::const_iterator refEnd  = m_indexData.end();
    for ( ; refIter != refEnd; ++refIter ) {
      
        const vector<BamToolsIndexEntry> offsets = (*refIter).second;
        
        // write number of offsets listed for this reference
        uint32_t numOffsets = offsets.size();
        if ( isBigEndian ) SwapEndian_32(numOffsets);
        fwrite(&numOffsets, sizeof(numOffsets), 1, indexStream);
       
        // iterate over offset entries
        vector<BamToolsIndexEntry>::const_iterator offsetIter = offsets.begin();
        vector<BamToolsIndexEntry>::const_iterator offsetEnd  = offsets.end();
        for ( ; offsetIter != offsetEnd; ++offsetIter ) {
          
            // get reference index data
            const BamToolsIndexEntry& entry = (*offsetIter);
            
            // copy entry data
            int32_t maxEndPosition = entry.MaxEndPosition;
            int64_t startOffset    = entry.StartOffset;
            int32_t startPosition  = entry.StartPosition;
            
            // swap endian-ness if necessary
            if ( isBigEndian ) {
                SwapEndian_32(maxEndPosition);
                SwapEndian_64(startOffset);
                SwapEndian_32(startPosition);
            }
            
            // write the reference index entry
            fwrite(&maxEndPosition, sizeof(maxEndPosition), 1, indexStream);
            fwrite(&startOffset,    sizeof(startOffset),    1, indexStream);
            fwrite(&startPosition,  sizeof(startPosition),  1, indexStream);
        }
    }

    // flush any remaining output, close file, and return success
    fflush(indexStream);
    fclose(indexStream);
    return true;
}

// ---------------------------------------------------
// BamToolsIndex implementation

BamToolsIndex::BamToolsIndex(BgzfData* bgzf, BamReader* reader, bool isBigEndian)
    : BamIndex(bgzf, reader, isBigEndian)
{ 
    d = new BamToolsIndexPrivate(this);
}    

BamToolsIndex::~BamToolsIndex(void) { 
    delete d;
    d = 0;
}

// creates index data (in-memory) from current reader data
bool BamToolsIndex::Build(void) { return d->Build(); }

// attempts to use index to jump to region; returns success/fail
// a "successful" jump indicates no error, but not whether this region has data
//   * thus, the method sets a flag to indicate whether there are alignments 
//     available after the jump position
bool BamToolsIndex::Jump(const BamRegion& region, bool* hasAlignmentsInRegion) { return d->Jump(region, hasAlignmentsInRegion); }

// loads existing data from file into memory
bool BamToolsIndex::Load(const string& filename) { return d->Load(filename); }

// writes in-memory index data out to file 
// N.B. - (this is the original BAM filename, method will modify it to use applicable extension)
bool BamToolsIndex::Write(const string& bamFilename) { return d->Write(bamFilename); }

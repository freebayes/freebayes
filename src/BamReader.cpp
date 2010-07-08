// ***************************************************************************
// BamReader.cpp (c) 2009 Derek Barnett, Michael Str�mberg
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ---------------------------------------------------------------------------
// Last modified: 22 June 2010 (DB)
// ---------------------------------------------------------------------------
// Uses BGZF routines were adapted from the bgzf.c code developed at the Broad
// Institute.
// ---------------------------------------------------------------------------
// Provides the basic functionality for reading BAM files
// ***************************************************************************

// C++ includes
#include <algorithm>
#include <iterator>
#include <string>
#include <vector>
#include <iostream>

// BamTools includes
#include "BGZF.h"
#include "BamReader.h"
using namespace BamTools;
using namespace std;

struct BamReader::BamReaderPrivate {

    // -------------------------------
    // structs, enums, typedefs
    // -------------------------------
    enum RegionState { BEFORE_REGION = 0
                      , WITHIN_REGION
                      , AFTER_REGION
          };

    // -------------------------------
    // data members
    // -------------------------------

    // general file data
    BgzfData  mBGZF;
    string    HeaderText;
    BamIndex  Index;
    RefVector References;
    bool      IsIndexLoaded;
    int64_t   AlignmentsBeginOffset;
    string    Filename;
    string    IndexFilename;
    
    // system data
    bool IsBigEndian;

    // user-specified region values
    BamRegion Region;
    bool IsLeftBoundSpecified;
    bool IsRightBoundSpecified;
    
    bool IsRegionSpecified;
    int  CurrentRefID;
    int  CurrentLeft;

    // BAM character constants
    const char* DNA_LOOKUP;
    const char* CIGAR_LOOKUP;

    // -------------------------------
    // constructor & destructor
    // -------------------------------
    BamReaderPrivate(void);
    ~BamReaderPrivate(void);

    // -------------------------------
    // "public" interface
    // -------------------------------

    // file operations
    void Close(void);
    bool Jump(int refID, int position = 0);
    bool Open(const string& filename, const string& indexFilename = "");
    bool Rewind(void);
    bool SetRegion(const BamRegion& region);

    // access alignment data
    bool GetNextAlignment(BamAlignment& bAlignment);
    bool GetNextAlignmentCore(BamAlignment& bAlignment);

    // access auxiliary data
    int GetReferenceID(const string& refName) const;

    // index operations
    bool CreateIndex(void);

    // -------------------------------
    // internal methods
    // -------------------------------

    // *** reading alignments and auxiliary data *** //

    // calculate bins that overlap region
    int BinsFromRegion(uint16_t bins[MAX_BIN]);
    // fills out character data for BamAlignment data
    bool BuildCharData(BamAlignment& bAlignment);
    // calculate file offset for first alignment chunk overlapping specified region
    int64_t GetOffset(std::vector<int64_t>& chunkStarts);
    // checks to see if alignment overlaps current region
    RegionState IsOverlap(BamAlignment& bAlignment);
    // retrieves header text from BAM file
    void LoadHeaderData(void);
    // retrieves BAM alignment under file pointer
    bool LoadNextAlignment(BamAlignment& bAlignment);
    // builds reference data structure from BAM file
    void LoadReferenceData(void);

    // *** index file handling *** //

    // calculates index for BAM file
    bool BuildIndex(void);
    // clear out inernal index data structure
    void ClearIndex(void);
    // saves BAM bin entry for index
    void InsertBinEntry(BamBinMap& binMap, const uint32_t& saveBin, const uint64_t& saveOffset, const uint64_t& lastOffset);
    // saves linear offset entry for index
    void InsertLinearOffset(LinearOffsetVector& offsets, const BamAlignment& bAlignment, const uint64_t& lastOffset);
    // loads index from BAM index file
    bool LoadIndex(void);
    // simplifies index by merging 'chunks'
    void MergeChunks(void);
    // saves index to BAM index file
    bool WriteIndex(void);
};

// -----------------------------------------------------
// BamReader implementation (wrapper around BRPrivate)
// -----------------------------------------------------

// constructor
BamReader::BamReader(void) {
    d = new BamReaderPrivate;
}

// destructor
BamReader::~BamReader(void) {
    delete d;
    d = 0;
}

// file operations
void BamReader::Close(void) { d->Close(); }
bool BamReader::Jump(int refID, int position) { 
    d->Region.LeftRefID = refID;
    d->Region.LeftPosition = position;
    d->IsLeftBoundSpecified = true;
    d->IsRightBoundSpecified = false;
    return d->Jump(refID, position); 
}
bool BamReader::Open(const string& filename, const string& indexFilename) { return d->Open(filename, indexFilename); }
bool BamReader::Rewind(void) { return d->Rewind(); }
bool BamReader::SetRegion(const BamRegion& region) { return d->SetRegion(region); }
bool BamReader::SetRegion(const int& leftRefID, const int& leftBound, const int& rightRefID, const int& rightBound) {
    return d->SetRegion( BamRegion(leftRefID, leftBound, rightRefID, rightBound) );
}

// access alignment data
bool BamReader::GetNextAlignment(BamAlignment& bAlignment) { return d->GetNextAlignment(bAlignment); }
bool BamReader::GetNextAlignmentCore(BamAlignment& bAlignment) { return d->GetNextAlignmentCore(bAlignment); }

// access auxiliary data
const string BamReader::GetHeaderText(void) const { return d->HeaderText; }
int BamReader::GetReferenceCount(void) const { return d->References.size(); }
const RefVector BamReader::GetReferenceData(void) const { return d->References; }
int BamReader::GetReferenceID(const string& refName) const { return d->GetReferenceID(refName); }
const std::string BamReader::GetFilename(void) const { return d->Filename; }

// index operations
bool BamReader::CreateIndex(void) { return d->CreateIndex(); }

// -----------------------------------------------------
// BamReaderPrivate implementation
// -----------------------------------------------------

// constructor
BamReader::BamReaderPrivate::BamReaderPrivate(void)
    : IsIndexLoaded(false)
    , AlignmentsBeginOffset(0)
    , IsLeftBoundSpecified(false)
    , IsRightBoundSpecified(false)
    , IsRegionSpecified(false)
    , CurrentRefID(0)
    , CurrentLeft(0)
    , DNA_LOOKUP("=ACMGRSVTWYHKDBN")
    , CIGAR_LOOKUP("MIDNSHP")
{ 
    IsBigEndian = SystemIsBigEndian();
}

// destructor
BamReader::BamReaderPrivate::~BamReaderPrivate(void) {
    Close();
}

// calculate bins that overlap region
int BamReader::BamReaderPrivate::BinsFromRegion(uint16_t list[MAX_BIN]) {

    // get region boundaries
    uint32_t begin = (unsigned int)Region.LeftPosition;
    uint32_t end;
    
    // if right bound specified AND left&right bounds are on same reference
    // OK to use right bound position
    if ( IsRightBoundSpecified && ( Region.LeftRefID == Region.RightRefID ) )
        end = (unsigned int)Region.RightPosition; // -1 ??
    
    // otherwise, use end of left bound reference as cutoff
    else
        end = (unsigned int)References.at(Region.LeftRefID).RefLength - 1;
    
    // initialize list, bin '0' always a valid bin
    int i = 0;
    list[i++] = 0;

    // get rest of bins that contain this region
    unsigned int k;
    for (k =    1 + (begin>>26); k <=    1 + (end>>26); ++k) { list[i++] = k; }
    for (k =    9 + (begin>>23); k <=    9 + (end>>23); ++k) { list[i++] = k; }
    for (k =   73 + (begin>>20); k <=   73 + (end>>20); ++k) { list[i++] = k; }
    for (k =  585 + (begin>>17); k <=  585 + (end>>17); ++k) { list[i++] = k; }
    for (k = 4681 + (begin>>14); k <= 4681 + (end>>14); ++k) { list[i++] = k; }

    // return number of bins stored
    return i;
}

bool BamReader::BamReaderPrivate::BuildCharData(BamAlignment& bAlignment) {
  
    // calculate character lengths/offsets
    const unsigned int dataLength     = bAlignment.SupportData.BlockLength - BAM_CORE_SIZE;
    const unsigned int seqDataOffset  = bAlignment.SupportData.QueryNameLength + (bAlignment.SupportData.NumCigarOperations * 4);
    const unsigned int qualDataOffset = seqDataOffset + (bAlignment.SupportData.QuerySequenceLength+1)/2;
    const unsigned int tagDataOffset  = qualDataOffset + bAlignment.SupportData.QuerySequenceLength;
    const unsigned int tagDataLength  = dataLength - tagDataOffset;
      
    // set up char buffers
    const char* allCharData = bAlignment.SupportData.AllCharData.data();
    const char* seqData     = ((const char*)allCharData) + seqDataOffset;
    const char* qualData    = ((const char*)allCharData) + qualDataOffset;
          char* tagData     = ((char*)allCharData) + tagDataOffset;
  
    // save query sequence
    bAlignment.QueryBases.clear();
    bAlignment.QueryBases.reserve(bAlignment.SupportData.QuerySequenceLength);
    for (unsigned int i = 0; i < bAlignment.SupportData.QuerySequenceLength; ++i) {
        char singleBase = DNA_LOOKUP[ ( ( seqData[(i/2)] >> (4*(1-(i%2)))) & 0xf ) ];
        bAlignment.QueryBases.append(1, singleBase);
    }
  
    // save qualities, converting from numeric QV to 'FASTQ-style' ASCII character
    bAlignment.Qualities.clear();
    bAlignment.Qualities.reserve(bAlignment.SupportData.QuerySequenceLength);
    for (unsigned int i = 0; i < bAlignment.SupportData.QuerySequenceLength; ++i) {
        char singleQuality = (char)(qualData[i]+33);
        bAlignment.Qualities.append(1, singleQuality);
    }
    
    // parse CIGAR to build 'AlignedBases'
    bAlignment.AlignedBases.clear();
    bAlignment.AlignedBases.reserve(bAlignment.SupportData.QuerySequenceLength);
    
    int k = 0;
    vector<CigarOp>::const_iterator cigarIter = bAlignment.CigarData.begin();
    vector<CigarOp>::const_iterator cigarEnd  = bAlignment.CigarData.end();
    for ( ; cigarIter != cigarEnd; ++cigarIter ) {
        
        const CigarOp& op = (*cigarIter);
        switch(op.Type) {
          
            case ('M') :
            case ('I') :
                bAlignment.AlignedBases.append(bAlignment.QueryBases.substr(k, op.Length)); // for 'M', 'I' - write bases
                // fall through
            
            case ('S') :
                k += op.Length;                                     // for 'S' - soft clip, skip over query bases
                break;
                
            case ('D') :
                bAlignment.AlignedBases.append(op.Length, '-');     // for 'D' - write gap character
                break;
                
            case ('P') :
                bAlignment.AlignedBases.append( op.Length, '*' );   // for 'P' - write padding character
                break;
                
            case ('N') :
                bAlignment.AlignedBases.append( op.Length, 'N' );  // for 'N' - write N's, skip bases in original query sequence
                // k+=op.Length; 
                break;
                
            case ('H') :
                break;  // for 'H' - hard clip, do nothing to AlignedBases, move to next op
                
            default:
                printf("ERROR: Invalid Cigar op type\n"); // shouldn't get here
                exit(1);
        }
    }
 
    // -----------------------
    // Added: 3-25-2010 DWB
    // Fixed: endian-correctness for tag data
    // -----------------------
    if ( IsBigEndian ) {
        int i = 0;
        while ( (unsigned int)i < tagDataLength ) {
          
            i += 2; // skip tag type (e.g. "RG", "NM", etc)
            uint8_t type = toupper(tagData[i]);     // lower & upper case letters have same meaning 
            ++i;                                    // skip value type
    
            switch (type) {
                
                case('A') :
                case('C') : 
                    ++i;
                    break;

                case('S') : 
                    SwapEndian_16p(&tagData[i]); 
                    i+=2; // sizeof(uint16_t)
                    break;
                    
                case('F') :
                case('I') : 
                    SwapEndian_32p(&tagData[i]);
                    i+=4; // sizeof(uint32_t)
                    break;
                
                case('D') : 
                    SwapEndian_64p(&tagData[i]);
                    i+=8; // sizeof(uint64_t) 
                    break;
                
                case('H') :
                case('Z') : 
                    while (tagData[i]) { ++i; }
                    ++i; // increment one more for null terminator
                    break;
                
                default : 
                    printf("ERROR: Invalid tag value type\n"); // shouldn't get here
                    exit(1);
            }
        }
    }
    
    // store TagData
    bAlignment.TagData.clear();
    bAlignment.TagData.resize(tagDataLength);
    memcpy((char*)bAlignment.TagData.data(), tagData, tagDataLength);
    
    // set support data parsed flag
    bAlignment.SupportData.IsParsed = true;
    
    // return success
    return true;
}

// populates BAM index data structure from BAM file data
bool BamReader::BamReaderPrivate::BuildIndex(void) {

    // check to be sure file is open
    if (!mBGZF.IsOpen) { return false; }

    // move file pointer to beginning of alignments
    Rewind();

    // get reference count, reserve index space
    int numReferences = References.size();
    for ( int i = 0; i < numReferences; ++i ) {
        Index.push_back(ReferenceIndex());
    }

    // sets default constant for bin, ID, offset, coordinate variables
    const uint32_t defaultValue = 0xffffffffu;

    // bin data
    uint32_t saveBin(defaultValue);
    uint32_t lastBin(defaultValue);

    // reference ID data
    int32_t saveRefID(defaultValue);
    int32_t lastRefID(defaultValue);

    // offset data
    uint64_t saveOffset = mBGZF.Tell();
    uint64_t lastOffset = saveOffset;

    // coordinate data
    int32_t lastCoordinate = defaultValue;

    BamAlignment bAlignment;
    while( GetNextAlignment(bAlignment) ) {

        // change of chromosome, save ID, reset bin
        if ( lastRefID != bAlignment.RefID ) {
            lastRefID = bAlignment.RefID;
            lastBin   = defaultValue;
        }

        // if lastCoordinate greater than BAM position - file not sorted properly
        else if ( lastCoordinate > bAlignment.Position ) {
            printf("BAM file not properly sorted:\n");
            printf("Alignment %s : %d > %d on reference (id = %d)", bAlignment.Name.c_str(), lastCoordinate, bAlignment.Position, bAlignment.RefID);
            exit(1);
        }

        // if valid reference && BAM bin spans some minimum cutoff (smaller bin ids span larger regions)
        if ( (bAlignment.RefID >= 0) && (bAlignment.Bin < 4681) ) {

            // save linear offset entry (matched to BAM entry refID)
            ReferenceIndex& refIndex = Index.at(bAlignment.RefID);
            LinearOffsetVector& offsets = refIndex.Offsets;
            InsertLinearOffset(offsets, bAlignment, lastOffset);
        }

        // if current BamAlignment bin != lastBin, "then possibly write the binning index"
        if ( bAlignment.Bin != lastBin ) {

            // if not first time through
            if ( saveBin != defaultValue ) {

                // save Bam bin entry
                ReferenceIndex& refIndex = Index.at(saveRefID);
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

            // if invalid RefID, break out (why?)
            if ( saveRefID < 0 ) { break; }
        }

        // make sure that current file pointer is beyond lastOffset
        if ( mBGZF.Tell() <= (int64_t)lastOffset  ) {
            printf("Error in BGZF offsets.\n");
            exit(1);
        }

        // update lastOffset
        lastOffset = mBGZF.Tell();

        // update lastCoordinate
        lastCoordinate = bAlignment.Position;
    }

    // save any leftover BAM data (as long as refID is valid)
    if ( saveRefID >= 0 ) {
        // save Bam bin entry
        ReferenceIndex& refIndex = Index.at(saveRefID);
        BamBinMap& binMap = refIndex.Bins;
        InsertBinEntry(binMap, saveBin, saveOffset, lastOffset);
    }

    // simplify index by merging chunks
    MergeChunks();

    // iterate over references
    BamIndex::iterator indexIter = Index.begin();
    BamIndex::iterator indexEnd  = Index.end();
    for ( int i = 0; indexIter != indexEnd; ++indexIter, ++i ) {

        // get reference index data
        ReferenceIndex& refIndex = (*indexIter);
        BamBinMap& binMap = refIndex.Bins;
        LinearOffsetVector& offsets = refIndex.Offsets;

        // store whether reference has alignments or no
        References[i].RefHasAlignments = ( binMap.size() > 0 );

        // sort linear offsets
        sort(offsets.begin(), offsets.end());
    }


    // rewind file pointer to beginning of alignments, return success/fail
    return Rewind();
}


// clear index data structure
void BamReader::BamReaderPrivate::ClearIndex(void) {
    Index.clear(); // sufficient ??
}

// closes the BAM file
void BamReader::BamReaderPrivate::Close(void) {
    mBGZF.Close();
    ClearIndex();
    HeaderText.clear();
    IsLeftBoundSpecified = false;
    IsRightBoundSpecified = false;
    IsRegionSpecified = false;
}

// create BAM index from BAM file (keep structure in memory) and write to default index output file
bool BamReader::BamReaderPrivate::CreateIndex(void) {

    // clear out index
    ClearIndex();

    // build (& save) index from BAM file
    bool ok = true;
    ok &= BuildIndex();
    ok &= WriteIndex();

    // return success/fail
    return ok;
}

// get next alignment (from specified region, if given)
bool BamReader::BamReaderPrivate::GetNextAlignment(BamAlignment& bAlignment) {

    // if valid alignment found, attempt to parse char data, and return success/failure
    if ( GetNextAlignmentCore(bAlignment) )
        return BuildCharData(bAlignment);
    
    // no valid alignment found
    else
        return false;
}

// retrieves next available alignment core data (returns success/fail)
// ** DOES NOT parse any character data (bases, qualities, tag data)
//    these can be accessed, if necessary, from the supportData 
// useful for operations requiring ONLY positional or other alignment-related information
bool BamReader::BamReaderPrivate::GetNextAlignmentCore(BamAlignment& bAlignment) {

    // if valid alignment available
    if ( LoadNextAlignment(bAlignment) ) {

        // if region not specified, return success
        if ( !IsLeftBoundSpecified ) return true;

        // determine region state (before, within, after)
        BamReader::BamReaderPrivate::RegionState state = IsOverlap(bAlignment);
      
        // if alignment lies after region, return false
        if ( state == AFTER_REGION ) 
            return false;

        while ( state != WITHIN_REGION ) {
            // if no valid alignment available (likely EOF) return failure
            if ( !LoadNextAlignment(bAlignment) ) return false;
            // if alignment lies after region, return false (no available read within region)
            state = IsOverlap(bAlignment);
            if ( state == AFTER_REGION) return false;
            
        }

        // return success (alignment found that overlaps region)
        return true;
    }

    // no valid alignment
    else
        return false;
}

// calculate file offset for first alignment chunk overlapping specified region
int64_t BamReader::BamReaderPrivate::GetOffset(std::vector<int64_t>& chunkStarts) {

    // calculate which bins overlap this region
    uint16_t* bins = (uint16_t*)calloc(MAX_BIN, 2);
    int numBins = BinsFromRegion(bins);

    // get bins for this reference
    const ReferenceIndex& refIndex = Index.at(Region.LeftRefID);
    const BamBinMap& binMap        = refIndex.Bins;

    // get minimum offset to consider
    const LinearOffsetVector& offsets = refIndex.Offsets;
    uint64_t minOffset = ( (unsigned int)(Region.LeftPosition>>BAM_LIDX_SHIFT) >= offsets.size() ) ? 0 : offsets.at(Region.LeftPosition>>BAM_LIDX_SHIFT);

    // store offsets to beginning of alignment 'chunks'
    //std::vector<int64_t> chunkStarts;

    // store all alignment 'chunk' starts for bins in this region
    for (int i = 0; i < numBins; ++i ) {
      
        const uint16_t binKey = bins[i];
        map<uint32_t, ChunkVector>::const_iterator binIter = binMap.find(binKey);
        if ( (binIter != binMap.end()) && ((*binIter).first == binKey) ) {

            const ChunkVector& chunks = (*binIter).second;
            std::vector<Chunk>::const_iterator chunksIter = chunks.begin();
            std::vector<Chunk>::const_iterator chunksEnd  = chunks.end();
            for ( ; chunksIter != chunksEnd; ++chunksIter) {
                const Chunk& chunk = (*chunksIter);
                if ( chunk.Stop > minOffset ) {
                    chunkStarts.push_back( chunk.Start );
                }
            }
        }
    }

    // clean up memory
    free(bins);

    // if no alignments found, else return smallest offset for alignment starts
    if ( chunkStarts.size() == 0 ) { return -1; }
    else { return *min_element(chunkStarts.begin(), chunkStarts.end()); }
}

// returns RefID for given RefName (returns References.size() if not found)
int BamReader::BamReaderPrivate::GetReferenceID(const string& refName) const {

    // retrieve names from reference data
    vector<string> refNames;
    RefVector::const_iterator refIter = References.begin();
    RefVector::const_iterator refEnd  = References.end();
    for ( ; refIter != refEnd; ++refIter) {
        refNames.push_back( (*refIter).RefName );
    }

    // return 'index-of' refName ( if not found, returns refNames.size() )
    return distance(refNames.begin(), find(refNames.begin(), refNames.end(), refName));
}

// saves BAM bin entry for index
void BamReader::BamReaderPrivate::InsertBinEntry(BamBinMap&      binMap,
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
void BamReader::BamReaderPrivate::InsertLinearOffset(LinearOffsetVector& offsets,
                                                     const BamAlignment& bAlignment,
                                                     const uint64_t&     lastOffset)
{
    // get converted offsets
    int beginOffset = bAlignment.Position >> BAM_LIDX_SHIFT;
    int endOffset   = (bAlignment.GetEndPosition() - 1) >> BAM_LIDX_SHIFT;

    // resize vector if necessary
    int oldSize = offsets.size();
    int newSize = endOffset + 1;
    if ( oldSize < newSize ) { offsets.resize(newSize, 0); }

    // store offset
    for(int i = beginOffset + 1; i <= endOffset ; ++i) {
        if ( offsets[i] == 0) {
            offsets[i] = lastOffset;
        }
    }
}

// returns region state - whether alignment ends before, overlaps, or starts after currently specified region
// this *internal* method should ONLY called when (at least) IsLeftBoundSpecified == true
BamReader::BamReaderPrivate::RegionState BamReader::BamReaderPrivate::IsOverlap(BamAlignment& bAlignment) {
    
    // --------------------------------------------------
    // check alignment start against right bound cutoff
    
    // if full region of interest was given
    if ( IsRightBoundSpecified ) {
      
        // read starts on right bound reference, but AFTER right bound position
        if ( bAlignment.RefID == Region.RightRefID && bAlignment.Position > Region.RightPosition )
            return AFTER_REGION;
      
        // if read starts on reference AFTER right bound, return false
        if ( bAlignment.RefID > Region.RightRefID ) 
            return AFTER_REGION;
    }
  
    // --------------------------------------------------------
    // no right bound given OR read starts before right bound
    // so, check if it overlaps left bound 
  
    // if read starts on left bound reference AND after left boundary, return success
    if ( bAlignment.RefID == Region.LeftRefID && bAlignment.Position >= Region.LeftPosition)
        return WITHIN_REGION;
  
    // if read is on any reference sequence before left bound, return false
    if ( bAlignment.RefID < Region.LeftRefID )
        return BEFORE_REGION;

    // --------------------------------------------------------
    // read is on left bound reference, but starts before left bound position

    // if it overlaps, return WITHIN_REGION
    if ( bAlignment.GetEndPosition() >= Region.LeftPosition )
        return WITHIN_REGION;
    // else begins before left bound position
    else
        return BEFORE_REGION;
}

// jumps to specified region(refID, leftBound) in BAM file, returns success/fail
bool BamReader::BamReaderPrivate::Jump(int refID, int position) {

    // if data exists for this reference and position is valid    
    if ( References.at(refID).RefHasAlignments && (position <= References.at(refID).RefLength) ) {

        // calculate offset
        std::vector<int64_t> chunkStarts;
        int64_t offset = GetOffset(chunkStarts);
	sort(chunkStarts.begin(), chunkStarts.end());

        // if in valid offset, return failure
        // otherwise return success of seek operation
        if ( offset == -1 ) {
            return false;
        } else {
            //return mBGZF.Seek(offset);
            BamAlignment bAlignment;
            bool result = true;
            for (std::vector<int64_t>::const_iterator o = chunkStarts.begin(); o != chunkStarts.end(); ++o) {
            //    std::cerr << *o << endl;
	    //	std::cerr << "Seeking to offset: " << *o << endl;
                result &= mBGZF.Seek(*o);
                LoadNextAlignment(bAlignment);
            //    std::cerr << "alignment: " << bAlignment.RefID 
            //        << ":" << bAlignment.Position << ".." << bAlignment.Position + bAlignment.Length << endl;
                if ((bAlignment.RefID == refID && bAlignment.Position + bAlignment.Length > position) || bAlignment.RefID > refID) {
             //       std::cerr << "here i am" << endl;
	     //	    std::cerr << "seeking to offset: " << (*(o-1)) << endl;
	     //	    std::cerr << "*** Finished jumping ***" << endl;
                    return mBGZF.Seek(*o);
                }
            }

            //std::cerr << "*** Finished jumping ***" << endl;
            return result;
        }
    }

    // invalid jump request parameters, return failure
    return false;
}

// load BAM header data
void BamReader::BamReaderPrivate::LoadHeaderData(void) {

    // check to see if proper BAM header
    char buffer[4];
    if (mBGZF.Read(buffer, 4) != 4) {
        printf("Could not read header type\n");
        exit(1);
    }

    if (strncmp(buffer, "BAM\001", 4)) {
        printf("wrong header type!\n");
        exit(1);
    }

    // get BAM header text length
    mBGZF.Read(buffer, 4);
    unsigned int headerTextLength = BgzfData::UnpackUnsignedInt(buffer);
    if ( IsBigEndian ) { SwapEndian_32(headerTextLength); }
    
    // get BAM header text
    char* headerText = (char*)calloc(headerTextLength + 1, 1);
    mBGZF.Read(headerText, headerTextLength);
    HeaderText = (string)((const char*)headerText);

    // clean up calloc-ed temp variable
    free(headerText);
}

// load existing index data from BAM index file (".bai"), return success/fail
bool BamReader::BamReaderPrivate::LoadIndex(void) {

    // clear out index data
    ClearIndex();

    // skip if index file empty
    if ( IndexFilename.empty() ) { return false; }

    // open index file, abort on error
    FILE* indexStream = fopen(IndexFilename.c_str(), "rb");
    if(!indexStream) {
        printf("ERROR: Unable to open the BAM index file %s for reading.\n", IndexFilename.c_str() );
        return false;
    }

    size_t elementsRead = 0;
        
    // see if index is valid BAM index
    char magic[4];
    elementsRead = fread(magic, 1, 4, indexStream);
    if (strncmp(magic, "BAI\1", 4)) {
        printf("Problem with index file - invalid format.\n");
        fclose(indexStream);
        return false;
    }

    // get number of reference sequences
    uint32_t numRefSeqs;
    elementsRead = fread(&numRefSeqs, 4, 1, indexStream);
    if ( IsBigEndian ) { SwapEndian_32(numRefSeqs); }
    
    // intialize space for BamIndex data structure
    Index.reserve(numRefSeqs);

    // iterate over reference sequences
    for (unsigned int i = 0; i < numRefSeqs; ++i) {

        // get number of bins for this reference sequence
        int32_t numBins;
        elementsRead = fread(&numBins, 4, 1, indexStream);
        if ( IsBigEndian ) { SwapEndian_32(numBins); }

        if (numBins > 0) {
            RefData& refEntry = References[i];
            refEntry.RefHasAlignments = true;
        }

        // intialize BinVector
        BamBinMap binMap;

        // iterate over bins for that reference sequence
        for (int j = 0; j < numBins; ++j) {

            // get binID
            uint32_t binID;
            elementsRead = fread(&binID, 4, 1, indexStream);

            // get number of regionChunks in this bin
            uint32_t numChunks;
            elementsRead = fread(&numChunks, 4, 1, indexStream);

            if ( IsBigEndian ) { 
              SwapEndian_32(binID);
              SwapEndian_32(numChunks);
            }
            
            // intialize ChunkVector
            ChunkVector regionChunks;
            regionChunks.reserve(numChunks);

            // iterate over regionChunks in this bin
            for (unsigned int k = 0; k < numChunks; ++k) {

                // get chunk boundaries (left, right)
                uint64_t left;
                uint64_t right;
                elementsRead = fread(&left, 8, 1, indexStream);
                elementsRead = fread(&right, 8, 1, indexStream);

                if ( IsBigEndian ) {
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

        // load linear index for this reference sequence

        // get number of linear offsets
        int32_t numLinearOffsets;
        elementsRead = fread(&numLinearOffsets, 4, 1, indexStream);
        if ( IsBigEndian ) { SwapEndian_32(numLinearOffsets); }

        // intialize LinearOffsetVector
        LinearOffsetVector offsets;
        offsets.reserve(numLinearOffsets);

        // iterate over linear offsets for this reference sequeence
        uint64_t linearOffset;
        for (int j = 0; j < numLinearOffsets; ++j) {
            // read a linear offset & store
            elementsRead = fread(&linearOffset, 8, 1, indexStream);
            if ( IsBigEndian ) { SwapEndian_64(linearOffset); }
            offsets.push_back(linearOffset);
        }

        // sort linear offsets
        sort( offsets.begin(), offsets.end() );

        // store index data for that reference sequence
        Index.push_back( ReferenceIndex(binMap, offsets) );
    }

    // close index file (.bai) and return
    fclose(indexStream);
    return true;
}

// populates BamAlignment with alignment data under file pointer, returns success/fail
bool BamReader::BamReaderPrivate::LoadNextAlignment(BamAlignment& bAlignment) {

    // read in the 'block length' value, make sure it's not zero
    char buffer[4];
    mBGZF.Read(buffer, 4);
    bAlignment.SupportData.BlockLength = BgzfData::UnpackUnsignedInt(buffer);
    if ( IsBigEndian ) { SwapEndian_32(bAlignment.SupportData.BlockLength); }
    if ( bAlignment.SupportData.BlockLength == 0 ) { return false; }

    // read in core alignment data, make sure the right size of data was read
    char x[BAM_CORE_SIZE];
    if ( mBGZF.Read(x, BAM_CORE_SIZE) != BAM_CORE_SIZE ) { return false; }

    if ( IsBigEndian ) {
        for ( int i = 0; i < BAM_CORE_SIZE; i+=sizeof(uint32_t) ) { 
          SwapEndian_32p(&x[i]); 
        }
    }
    
    // set BamAlignment 'core' and 'support' data
    bAlignment.RefID    = BgzfData::UnpackSignedInt(&x[0]);  
    bAlignment.Position = BgzfData::UnpackSignedInt(&x[4]);
    
    unsigned int tempValue = BgzfData::UnpackUnsignedInt(&x[8]);
    bAlignment.Bin        = tempValue >> 16;
    bAlignment.MapQuality = tempValue >> 8 & 0xff;
    bAlignment.SupportData.QueryNameLength = tempValue & 0xff;

    tempValue = BgzfData::UnpackUnsignedInt(&x[12]);
    bAlignment.AlignmentFlag = tempValue >> 16;
    bAlignment.SupportData.NumCigarOperations = tempValue & 0xffff;

    bAlignment.SupportData.QuerySequenceLength = BgzfData::UnpackUnsignedInt(&x[16]);
    bAlignment.MateRefID    = BgzfData::UnpackSignedInt(&x[20]);
    bAlignment.MatePosition = BgzfData::UnpackSignedInt(&x[24]);
    bAlignment.InsertSize   = BgzfData::UnpackSignedInt(&x[28]);
    
    // store 'all char data' and cigar ops
    const unsigned int dataLength      = bAlignment.SupportData.BlockLength - BAM_CORE_SIZE;
    const unsigned int cigarDataOffset = bAlignment.SupportData.QueryNameLength;
    
    char*     allCharData = (char*)calloc(sizeof(char), dataLength);
    uint32_t* cigarData   = (uint32_t*)(allCharData + cigarDataOffset);
    
    // read in character data - make sure proper data size was read
    if ( mBGZF.Read(allCharData, dataLength) != (signed int)dataLength) { return false; }
    else {
     
        // store alignment name and length
        bAlignment.Name.assign((const char*)(allCharData));
        bAlignment.Length = bAlignment.SupportData.QuerySequenceLength;
      
        // store remaining 'allCharData' in supportData structure
        bAlignment.SupportData.AllCharData.assign((const char*)allCharData, dataLength);
        
        // save CigarOps for BamAlignment
        CigarOp op;
        bAlignment.CigarData.clear();
        bAlignment.CigarData.reserve(bAlignment.SupportData.NumCigarOperations);
        for (unsigned int i = 0; i < bAlignment.SupportData.NumCigarOperations; ++i) {

            // swap if necessary
            if ( IsBigEndian ) { SwapEndian_32(cigarData[i]); }
          
            // build CigarOp structure
            op.Length = (cigarData[i] >> BAM_CIGAR_SHIFT);
            op.Type   = CIGAR_LOOKUP[ (cigarData[i] & BAM_CIGAR_MASK) ];

            // save CigarOp
            bAlignment.CigarData.push_back(op);
        }
    }

    free(allCharData);
    return true;
}

// loads reference data from BAM file
void BamReader::BamReaderPrivate::LoadReferenceData(void) {

    // get number of reference sequences
    char buffer[4];
    mBGZF.Read(buffer, 4);
    unsigned int numberRefSeqs = BgzfData::UnpackUnsignedInt(buffer);
    if ( IsBigEndian ) { SwapEndian_32(numberRefSeqs); }
    if (numberRefSeqs == 0) { return; }
    References.reserve((int)numberRefSeqs);

    // iterate over all references in header
    for (unsigned int i = 0; i != numberRefSeqs; ++i) {

        // get length of reference name
        mBGZF.Read(buffer, 4);
        unsigned int refNameLength = BgzfData::UnpackUnsignedInt(buffer);
        if ( IsBigEndian ) { SwapEndian_32(refNameLength); }
        char* refName = (char*)calloc(refNameLength, 1);

        // get reference name and reference sequence length
        mBGZF.Read(refName, refNameLength);
        mBGZF.Read(buffer, 4);
        int refLength = BgzfData::UnpackSignedInt(buffer);
        if ( IsBigEndian ) { SwapEndian_32(refLength); }

        // store data for reference
        RefData aReference;
        aReference.RefName   = (string)((const char*)refName);
        aReference.RefLength = refLength;
        References.push_back(aReference);

        // clean up calloc-ed temp variable
        free(refName);
    }
}

// merges 'alignment chunks' in BAM bin (used for index building)
void BamReader::BamReaderPrivate::MergeChunks(void) {

    // iterate over reference enties
    BamIndex::iterator indexIter = Index.begin();
    BamIndex::iterator indexEnd  = Index.end();
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
            if ( binChunks.size() == 0 ) { continue; }

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

// opens BAM file (and index)
bool BamReader::BamReaderPrivate::Open(const string& filename, const string& indexFilename) {

    Filename = filename;
    IndexFilename = indexFilename;

    // open the BGZF file for reading, return false on failure
    if ( !mBGZF.Open(filename, "rb") ) 
        return false;
    
    // retrieve header text & reference data
    LoadHeaderData();
    LoadReferenceData();

    // store file offset of first alignment
    AlignmentsBeginOffset = mBGZF.Tell();

    // open index file & load index data (if exists)
    if ( !IndexFilename.empty() ) {
        LoadIndex();
    }
    
    // return success
    return true;
}

// returns BAM file pointer to beginning of alignment data
bool BamReader::BamReaderPrivate::Rewind(void) {

    // find first reference that has alignments in the BAM file
    int refID = 0;
    int refCount = References.size();
    for ( ; refID < refCount; ++refID ) {
        if ( References.at(refID).RefHasAlignments ) { break; }
    }

    // reset default region info
    Region.LeftRefID = refID;
    Region.LeftPosition = 0;
    Region.RightRefID = -1;
    Region.RightPosition = -1;
    IsLeftBoundSpecified = false;
    IsRightBoundSpecified = false;

    // return success/failure of seek
    return mBGZF.Seek(AlignmentsBeginOffset);
}

// sets a region of interest (with left & right bound reference/position)
// attempts a Jump() to left bound as well
// returns success/failure of Jump()
bool BamReader::BamReaderPrivate::SetRegion(const BamRegion& region) {
    
    // save region of interest
    Region = region;
    
    // set flags
    if ( region.LeftRefID >= 0 && region.LeftPosition >= 0 ) 
        IsLeftBoundSpecified = true;
    if ( region.RightRefID >= 0 && region.RightPosition >= 0 ) 
        IsRightBoundSpecified = true;
    
    // attempt jump to beginning of region, return success/fail of Jump()
    return Jump( Region.LeftRefID, Region.LeftPosition );
}

// saves index data to BAM index file (".bai"), returns success/fail
bool BamReader::BamReaderPrivate::WriteIndex(void) {

    IndexFilename = Filename + ".bai";
    FILE* indexStream = fopen(IndexFilename.c_str(), "wb");
    if ( indexStream == 0 ) {
        printf("ERROR: Could not open file to save index\n");
        return false;
    }

    // write BAM index header
    fwrite("BAI\1", 1, 4, indexStream);

    // write number of reference sequences
    int32_t numReferenceSeqs = Index.size();
    if ( IsBigEndian ) { SwapEndian_32(numReferenceSeqs); }
    fwrite(&numReferenceSeqs, 4, 1, indexStream);

    // iterate over reference sequences
    BamIndex::const_iterator indexIter = Index.begin();
    BamIndex::const_iterator indexEnd  = Index.end();
    for ( ; indexIter != indexEnd; ++ indexIter ) {

        // get reference index data
        const ReferenceIndex& refIndex = (*indexIter);
        const BamBinMap& binMap = refIndex.Bins;
        const LinearOffsetVector& offsets = refIndex.Offsets;

        // write number of bins
        int32_t binCount = binMap.size();
        if ( IsBigEndian ) { SwapEndian_32(binCount); }
        fwrite(&binCount, 4, 1, indexStream);

        // iterate over bins
        BamBinMap::const_iterator binIter = binMap.begin();
        BamBinMap::const_iterator binEnd  = binMap.end();
        for ( ; binIter != binEnd; ++binIter ) {

            // get bin data (key and chunk vector)
            uint32_t binKey = (*binIter).first;
            const ChunkVector& binChunks = (*binIter).second;

            // save BAM bin key
            if ( IsBigEndian ) { SwapEndian_32(binKey); }
            fwrite(&binKey, 4, 1, indexStream);

            // save chunk count
            int32_t chunkCount = binChunks.size();
            if ( IsBigEndian ) { SwapEndian_32(chunkCount); }
            fwrite(&chunkCount, 4, 1, indexStream);

            // iterate over chunks
            ChunkVector::const_iterator chunkIter = binChunks.begin();
            ChunkVector::const_iterator chunkEnd  = binChunks.end();
            for ( ; chunkIter != chunkEnd; ++chunkIter ) {

                // get current chunk data
                const Chunk& chunk    = (*chunkIter);
                uint64_t start = chunk.Start;
                uint64_t stop  = chunk.Stop;

                if ( IsBigEndian ) {
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
        if ( IsBigEndian ) { SwapEndian_32(offsetSize); }
        fwrite(&offsetSize, 4, 1, indexStream);

        // iterate over linear offsets
        LinearOffsetVector::const_iterator offsetIter = offsets.begin();
        LinearOffsetVector::const_iterator offsetEnd  = offsets.end();
        for ( ; offsetIter != offsetEnd; ++offsetIter ) {

            // write linear offset value
            uint64_t linearOffset = (*offsetIter);
            if ( IsBigEndian ) { SwapEndian_64(linearOffset); }
            fwrite(&linearOffset, 8, 1, indexStream);
        }
    }

    // flush buffer, close file, and return success
    fflush(indexStream);
    fclose(indexStream);
    return true;
}

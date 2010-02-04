// ***************************************************************************
// BamReader.cpp (c) 2009 Derek Barnett, Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ---------------------------------------------------------------------------
// Last modified: 15 July 2009 (DB)
// ---------------------------------------------------------------------------
// The BGZF routines were adapted from the bgzf.c code developed at the Broad
// Institute.
// ---------------------------------------------------------------------------
// Provides the basic functionality for reading BAM files
// ***************************************************************************

// BamTools includes
#include "BamReader.h"
using namespace BamTools;
using namespace std;

// static character constants
const char* BamReader::DNA_LOOKUP   = "=ACMGRSVTWYHKDBN";
const char* BamReader::CIGAR_LOOKUP = "MIDNSHP";

// constructor
BamReader::BamReader(void)
	: m_index(NULL)
	, m_isIndexLoaded(false)
	, m_alignmentsBeginOffset(0)
	, m_isRegionSpecified(false)
	, m_currentRefID(0)
	, m_currentLeft(0)
{ }

// destructor
BamReader::~BamReader(void) {
	Close();
}

// checks BGZF block header
bool BamReader::BgzfCheckBlockHeader(char* header) {

	return (header[0] == GZIP_ID1 &&
            header[1] == (char)GZIP_ID2 &&
            header[2] == Z_DEFLATED &&
            (header[3] & FLG_FEXTRA) != 0 &&
            BgzfUnpackUnsignedShort(&header[10]) == BGZF_XLEN &&
            header[12] == BGZF_ID1 &&
            header[13] == BGZF_ID2 &&
            BgzfUnpackUnsignedShort(&header[14]) == BGZF_LEN
		   );
}

// closes the BAM file
void BamReader::BgzfClose(void) {
	m_BGZF.IsOpen = false;
	fclose(m_BGZF.Stream);
}

// de-compresses the current block
int BamReader::BgzfInflateBlock(int blockLength) {
	
	// Inflate the block in m_BGZF.CompressedBlock into m_BGZF.UncompressedBlock
    z_stream zs;
    zs.zalloc    = NULL;
    zs.zfree     = NULL;
    zs.next_in   = (Bytef*)m_BGZF.CompressedBlock + 18;
    zs.avail_in  = blockLength - 16;
    zs.next_out  = (Bytef*)m_BGZF.UncompressedBlock;
    zs.avail_out = m_BGZF.UncompressedBlockSize;

    int status = inflateInit2(&zs, GZIP_WINDOW_BITS);
    if (status != Z_OK) {
        printf("inflateInit failed\n");
        exit(1);
    }

    status = inflate(&zs, Z_FINISH);
    if (status != Z_STREAM_END) {
        inflateEnd(&zs);
        printf("inflate failed\n");
        exit(1);
    }

    status = inflateEnd(&zs);
    if (status != Z_OK) {
        printf("inflateEnd failed\n");
        exit(1);
    }

    return zs.total_out;
}

// opens the BAM file for reading
void BamReader::BgzfOpen(const string& filename) {

	m_BGZF.Stream = fopen(filename.c_str(), "rb");
	if(!m_BGZF.Stream) {
		printf("ERROR: Unable to open the BAM file %s for reading.\n", filename.c_str() );
		exit(1);
	}

	m_BGZF.IsOpen = true;
}

// reads BGZF data into buffer
unsigned int BamReader::BgzfRead(char* data, const unsigned int dataLength) {

    if (dataLength == 0) { return 0; }

	char* output = data;
    unsigned int numBytesRead = 0;
    while (numBytesRead < dataLength) {

        int bytesAvailable = m_BGZF.BlockLength - m_BGZF.BlockOffset;
        if (bytesAvailable <= 0) {
            if ( BgzfReadBlock() != 0 ) { return -1; }
            bytesAvailable = m_BGZF.BlockLength - m_BGZF.BlockOffset;
            if ( bytesAvailable <= 0 ) { break; }
        }

		char* buffer   = m_BGZF.UncompressedBlock;
        int copyLength = min( (int)(dataLength-numBytesRead), bytesAvailable );
        memcpy(output, buffer + m_BGZF.BlockOffset, copyLength);

        m_BGZF.BlockOffset += copyLength;
        output             += copyLength;
        numBytesRead       += copyLength;
    }

    if ( m_BGZF.BlockOffset == m_BGZF.BlockLength ) {
		m_BGZF.BlockAddress = ftello64(m_BGZF.Stream);						
        m_BGZF.BlockOffset  = 0;
        m_BGZF.BlockLength  = 0;
    }

	return numBytesRead;
}

int BamReader::BgzfReadBlock(void) {

    char    header[BLOCK_HEADER_LENGTH];
    int64_t blockAddress = ftello(m_BGZF.Stream);

    int count = fread(header, 1, sizeof(header), m_BGZF.Stream);
	if (count == 0) {
        m_BGZF.BlockLength = 0;
        return 0;
    }

    if (count != sizeof(header)) {
        printf("read block failed - count != sizeof(header)\n");
        return -1;
    }

    if (!BgzfCheckBlockHeader(header)) {
        printf("read block failed - CheckBgzfBlockHeader() returned false\n");
        return -1;
    }

    int blockLength = BgzfUnpackUnsignedShort(&header[16]) + 1;
    char* compressedBlock = m_BGZF.CompressedBlock;
    memcpy(compressedBlock, header, BLOCK_HEADER_LENGTH);
    int remaining = blockLength - BLOCK_HEADER_LENGTH;

    count = fread(&compressedBlock[BLOCK_HEADER_LENGTH], 1, remaining, m_BGZF.Stream);
    if (count != remaining) {
        printf("read block failed - count != remaining\n");
        return -1;
    }

    count = BgzfInflateBlock(blockLength);
    if (count < 0) { return -1; }

    if (m_BGZF.BlockLength != 0) {
        m_BGZF.BlockOffset = 0;
    }

    m_BGZF.BlockAddress = blockAddress;
    m_BGZF.BlockLength  = count;
    return 0;
}

// move file pointer to specified offset
bool BamReader::BgzfSeek(int64_t position) {

	int     blockOffset  = (position & 0xFFFF);
    int64_t blockAddress = (position >> 16) & 0xFFFFFFFFFFFFLL;
    if (fseeko(m_BGZF.Stream, blockAddress, SEEK_SET) != 0) {
        printf("ERROR: Unable to seek in BAM file\n");
		exit(1);
    }

    m_BGZF.BlockLength  = 0;
    m_BGZF.BlockAddress = blockAddress;
    m_BGZF.BlockOffset  = blockOffset;
	return true;
}

// get file position in BAM file
int64_t BamReader::BgzfTell(void) {
	return ( (m_BGZF.BlockAddress << 16) | (m_BGZF.BlockOffset & 0xFFFF) );
}

int BamReader::BinsFromRegion(int refID, unsigned int left, uint16_t list[MAX_BIN]) {

	// get region boundaries
	uint32_t begin = left;
	uint32_t end   = m_references.at(refID).RefLength - 1;

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

unsigned int BamReader::CalculateAlignmentEnd(const unsigned int& position, const vector<CigarOp>& cigarData) {

	// initialize alignment end to starting position
	unsigned int alignEnd = position;

	// iterate over cigar operations
	vector<CigarOp>::const_iterator cigarIter = cigarData.begin();
	vector<CigarOp>::const_iterator cigarEnd  = cigarData.end();
	for ( ; cigarIter != cigarEnd; ++cigarIter) {
		char cigarType = (*cigarIter).Type;
		if ( cigarType == 'M' || cigarType == 'D' || cigarType == 'N' ) {
			alignEnd += (*cigarIter).Length;
		}
	}
	return alignEnd;
}

void BamReader::ClearIndex(void) {

	if ( m_index ) {
		// iterate over references
		vector<RefIndex*>::iterator refIter = m_index->begin();
		vector<RefIndex*>::iterator refEnd  = m_index->end();
		for ( ; refIter != refEnd; ++refIter) {
			RefIndex* aRef = (*refIter);
			if ( aRef ) {
				// clear out BAM bins
				if ( aRef->first ) {
					BinVector::iterator binIter = (aRef->first)->begin();
					BinVector::iterator binEnd  = (aRef->first)->end();
					for ( ; binIter != binEnd; ++binIter ) {
						ChunkVector* chunks = (*binIter).second;
						if ( chunks ) { delete chunks; chunks = NULL;}
					}
					delete aRef->first;
					aRef->first = NULL;
				}
				// clear BAM linear offsets
				if ( aRef->second ) { delete aRef->second; aRef->second = NULL; }
				delete aRef;
				aRef = NULL;
			}
		}
		delete m_index;
		m_index = NULL;
	}
}

// closes the BAM file
void BamReader::Close(void) {
	if(m_BGZF.IsOpen) { BgzfClose(); }	
	ClearIndex();
	m_headerText.clear();
	m_isRegionSpecified = false;
}

const string BamReader::GetHeaderText(void) const {
	return m_headerText;
}

const int BamReader::GetReferenceCount(void) const {
	return m_references.size();
}

const RefVector BamReader::GetReferenceData(void) const {
	return m_references;
}

const int BamReader::GetReferenceID(const string& refName) const {

	// retrieve names from reference data
	vector<string> refNames;
	RefVector::const_iterator refIter = m_references.begin();
    RefVector::const_iterator refEnd  = m_references.end();
    for ( ; refIter != refEnd; ++refIter) {
		refNames.push_back( (*refIter).RefName );
    }

	// return 'index-of' refName ( if not found, returns refNames.size() )
	return Index( refNames.begin(), refNames.end(), refName );
}

// get next alignment (from specified region, if given)
bool BamReader::GetNextAlignment(BamAlignment& bAlignment) {

	// if valid alignment available
	if ( LoadNextAlignment(bAlignment) ) {
		
		// if region not specified, return success
		if ( !m_isRegionSpecified ) { return true; }
		
		// load next alignment until region overlap is found
		while ( !IsOverlap(bAlignment) ) {
			// if no valid alignment available (likely EOF) return failure
			if ( !LoadNextAlignment(bAlignment) ) { return false; }
		}
		
		// return success (alignment found that overlaps region)
		return true;
	} 
	
	// no valid alignment
	else { return false; }
}

int64_t BamReader::GetOffset(int refID, unsigned int left) {

	// calculate which bins overlap this region
	uint16_t* bins = (uint16_t*)calloc(MAX_BIN, 2);
	int numBins = BinsFromRegion(refID, left, bins);

	// get bins for this reference
	RefIndex* refIndex = m_index->at(refID);
	BinVector* refBins = refIndex->first;

	// get minimum offset to consider
	LinearOffsetVector* linearOffsets = refIndex->second;
	uint64_t minOffset = ((left>>BAM_LIDX_SHIFT) >= linearOffsets->size()) ? 0 : linearOffsets->at(left>>BAM_LIDX_SHIFT);

	// store offsets to beginning of alignment 'chunks' 
	std::vector<int64_t> chunkStarts;
	
	// reference bin iterators
	BinVector::const_iterator binIter;
	BinVector::const_iterator binBegin = refBins->begin();
	BinVector::const_iterator binEnd   = refBins->end();
	
	// store all alignment 'chunk' starts for bins in this region
	for (int i = 0; i < numBins; ++i ) {
		binIter = lower_bound(binBegin, binEnd, bins[i], LookupKeyCompare<uint32_t, ChunkVector*>() );
		if ( (binIter != binEnd) && ( (*binIter).first == bins[i]) ) {
			ChunkVector* binChunks = (*binIter).second;
			ChunkVector::const_iterator chunkIter = binChunks->begin();
			ChunkVector::const_iterator chunkEnd  = binChunks->end();
			for ( ; chunkIter != chunkEnd; ++chunkIter) {
				if ( (*chunkIter).second > minOffset ) {
					chunkStarts.push_back( (*chunkIter).first );
				}	
			}
		}
	}
	
	// clean up memory
	free(bins);
	
	// if no alignments found
	if ( chunkStarts.empty() ) { return -1; }
	
	// else return smallest offset for alignment starts 
	else { return *min_element(chunkStarts.begin(), chunkStarts.end()); }
}

bool BamReader::IsOverlap(BamAlignment& bAlignment) {

	// if on different reference sequence, quit
	if ( bAlignment.RefID != (unsigned int)m_currentRefID ) { return false; }

	// read starts after left boundary
	if ( bAlignment.Position >= m_currentLeft) { return true; }

	// return whether alignment end overlaps left boundary
	return ( CalculateAlignmentEnd(bAlignment.Position, bAlignment.CigarData) >= m_currentLeft );
}

bool BamReader::Jump(int refID, unsigned int position) {

	// if index available, and region is valid
	if ( (m_index->size() != 0) && m_references.at(refID).RefHasAlignments && (position <= m_references.at(refID).RefLength) ) { 
		m_currentRefID = refID;
		m_currentLeft  = position;
		m_isRegionSpecified = true;
		
		int64_t offset = GetOffset(m_currentRefID, m_currentLeft);
		if ( offset == -1 ) { return false; }
		else { return BgzfSeek(offset); }
	}
	return false;
}

void BamReader::LoadHeaderData(void) {
	
	// check to see if proper BAM header
	char buffer[4];
	if (BgzfRead(buffer, 4) != 4) { 
		printf("Could not read header type\n");
		exit(1); 
	}
	if (strncmp(buffer, "BAM\001", 4)) {
		printf("wrong header type!\n");
		exit(1);
	}
	
	// get BAM header text length
	BgzfRead(buffer, 4);
	const unsigned int headerTextLength = BgzfUnpackUnsignedInt(buffer);

	// get BAM header text
	char* headerText = (char*)calloc(headerTextLength + 1, 1);
	BgzfRead(headerText, headerTextLength);
	m_headerText = (string)((const char*)headerText);
	
	// clean up calloc-ed temp variable
	free(headerText);
}

void BamReader::LoadIndexData(FILE* indexStream) {

	// see if index is valid BAM index
	char magic[4];
	fread(magic, 1, 4, indexStream);
	if (strncmp(magic, "BAI\1", 4)) {
		printf("Problem with index file - invalid format.\n");
		fclose(indexStream);
		exit(1);
	}

	// get number of reference sequences
	uint32_t numRefSeqs;
	fread(&numRefSeqs, 4, 1, indexStream);
	
	// intialize BamIndex data structure
	m_index = new BamIndex;
	m_index->reserve(numRefSeqs);

	// iterate over reference sequences
	for (unsigned int i = 0; i < numRefSeqs; ++i) {
		
		// get number of bins for this reference sequence
		int32_t numBins;
		fread(&numBins, 4, 1, indexStream);
		
		if (numBins > 0) { m_references.at(i).RefHasAlignments = true; }

		// intialize BinVector
		BinVector* bins = new BinVector;
		bins->reserve(numBins);
		
		// iterate over bins for that reference sequence
		for (int j = 0; j < numBins; ++j) {
			
			// get binID 
			uint32_t binID;
			fread(&binID, 4, 1, indexStream);
			
			// get number of regionChunks in this bin
			uint32_t numChunks;
			fread(&numChunks, 4, 1, indexStream);
			
			// intialize ChunkVector
			ChunkVector* regionChunks = new ChunkVector;
			regionChunks->reserve(numChunks);
			
			// iterate over regionChunks in this bin
			for (unsigned int k = 0; k < numChunks; ++k) {
				
				// get chunk boundaries (left, right) 
				uint64_t left;
				uint64_t right;
				fread(&left, 8, 1, indexStream);
				fread(&right, 8, 1, indexStream);
				
				// save ChunkPair
				regionChunks->push_back( ChunkPair(left, right) );
			}
			
			// sort chunks for this bin
			sort( regionChunks->begin(), regionChunks->end(), LookupKeyCompare<uint64_t, uint64_t>() );

			// save binID, chunkVector for this bin
			bins->push_back( BamBin(binID, regionChunks) );
		}
		
		// sort bins by binID
		sort(bins->begin(), bins->end(), LookupKeyCompare<uint32_t, ChunkVector*>() );

		// load linear index for this reference sequence
		
		// get number of linear offsets
		int32_t numLinearOffsets;
		fread(&numLinearOffsets, 4, 1, indexStream);
		
		// intialize LinearOffsetVector
		LinearOffsetVector* linearOffsets = new LinearOffsetVector;
		linearOffsets->reserve(numLinearOffsets);
		
		// iterate over linear offsets for this reference sequeence
		for (int j = 0; j < numLinearOffsets; ++j) {
			// get a linear offset
			uint64_t linearOffset;
			fread(&linearOffset, 8, 1, indexStream);
			// store linear offset
			linearOffsets->push_back(linearOffset);
		}
		
		// sort linear offsets
		sort( linearOffsets->begin(), linearOffsets->end() );

		// store index data for that reference sequence
		m_index->push_back( new RefIndex(bins, linearOffsets) );
	}
	
	// close index file (.bai) and return
	fclose(indexStream);
}

bool BamReader::LoadNextAlignment(BamAlignment& bAlignment) {

	// read in the 'block length' value, make sure it's not zero
	char buffer[4];
	BgzfRead(buffer, 4);
	const unsigned int blockLength = BgzfUnpackUnsignedInt(buffer);
	if ( blockLength == 0 ) { return false; }

	// keep track of bytes read as method progresses
	int bytesRead = 4;

	// read in core alignment data, make sure the right size of data was read 
	char x[BAM_CORE_SIZE];
	if ( BgzfRead(x, BAM_CORE_SIZE) != BAM_CORE_SIZE ) { return false; }
	bytesRead += BAM_CORE_SIZE;

	// set BamAlignment 'core' data and character data lengths
	unsigned int tempValue;
	unsigned int queryNameLength;
	unsigned int numCigarOperations;
	unsigned int querySequenceLength;

	bAlignment.RefID    = BgzfUnpackUnsignedInt(&x[0]);
	bAlignment.Position = BgzfUnpackUnsignedInt(&x[4]);

	tempValue             = BgzfUnpackUnsignedInt(&x[8]);		
	bAlignment.Bin        = tempValue >> 16;
	bAlignment.MapQuality = tempValue >> 8 & 0xff;
	queryNameLength       = tempValue & 0xff;

	tempValue                = BgzfUnpackUnsignedInt(&x[12]);	
	bAlignment.AlignmentFlag = tempValue >> 16;
	numCigarOperations       = tempValue & 0xffff;

	querySequenceLength     = BgzfUnpackUnsignedInt(&x[16]);
	bAlignment.MateRefID    = BgzfUnpackUnsignedInt(&x[20]);
	bAlignment.MatePosition = BgzfUnpackUnsignedInt(&x[24]);
	bAlignment.InsertSize   = BgzfUnpackUnsignedInt(&x[28]);

	// calculate lengths/offsets
	const unsigned int dataLength      = blockLength - BAM_CORE_SIZE;
	const unsigned int cigarDataOffset = queryNameLength;
	const unsigned int seqDataOffset   = cigarDataOffset + (numCigarOperations * 4);
	const unsigned int qualDataOffset  = seqDataOffset + (querySequenceLength+1)/2;
	const unsigned int tagDataOffset   = qualDataOffset + querySequenceLength;
	const unsigned int tagDataLen      = dataLength - tagDataOffset;
	
	// set up destination buffers for character data
	char* allCharData   = (char*)calloc(sizeof(char), dataLength);
	uint32_t* cigarData = (uint32_t*)(allCharData + cigarDataOffset);
	char* seqData       = ((char*)allCharData) + seqDataOffset;
	char* qualData      = ((char*)allCharData) + qualDataOffset;
	char* tagData       = ((char*)allCharData) + tagDataOffset;
	
	// get character data - make sure proper data size was read
	if ( BgzfRead(allCharData, dataLength) != dataLength) { return false; }
	else {
		
		bytesRead += dataLength;
		
		// clear out any previous string data
		bAlignment.Name.clear();
		bAlignment.QueryBases.clear();
		bAlignment.Qualities.clear();
		bAlignment.AlignedBases.clear();
		bAlignment.CigarData.clear();
		bAlignment.TagData.clear();
		
		// save name
		bAlignment.Name = (string)((const char*)(allCharData));
		
		// save query sequence
		for (unsigned int i = 0; i < querySequenceLength; ++i) { 	
			char singleBase = DNA_LOOKUP[ ( ( seqData[(i/2)] >> (4*(1-(i%2)))) & 0xf ) ];
			bAlignment.QueryBases.append( 1, singleBase );
		}
		
		// save sequence length
		bAlignment.Length = bAlignment.QueryBases.length();
		
		// save qualities
		for (unsigned int i = 0; i < querySequenceLength; ++i) { 
			char singleQuality = (char)(qualData[i]+33);			// conversion from QV to FASTQ character
			bAlignment.Qualities.append( 1, singleQuality );
		}
		
		// save CIGAR-related data;
		int k = 0;
		for (unsigned int i = 0; i < numCigarOperations; ++i) {
			
			// build CigarOp struct
			CigarOp op;
			op.Length = (cigarData[i] >> BAM_CIGAR_SHIFT);
			op.Type   = CIGAR_LOOKUP[ (cigarData[i] & BAM_CIGAR_MASK) ];
			
			// save CigarOp
			bAlignment.CigarData.push_back(op);
			
			// build AlignedBases string
			switch (op.Type) {
				
				case ('M') : 
				case ('I') : bAlignment.AlignedBases.append( bAlignment.QueryBases.substr(k, op.Length) ); 	// for 'M', 'I' - write bases
				case ('S') : k += op.Length;                                                     			// for 'S' - skip over query bases
							 break;
							 
				case ('D') : bAlignment.AlignedBases.append( op.Length, '-' );	// for 'D' - write gap character
							 break;
							
				case ('P') : bAlignment.AlignedBases.append( op.Length, '*' );	// for 'P' - write padding character;
							 break;
							 
				case ('N') : bAlignment.AlignedBases.append( op.Length, 'N' );  // for 'N' - write N's, skip bases in query sequence
							 k += op.Length;
							 break;
							 
				case ('H') : break; 											// for 'H' - do nothing, move to next op
							 
				default    : printf("ERROR: Invalid Cigar op type\n");			// shouldn't get here
							 exit(1);
			}
		}
		
		// read in the tag data
		bAlignment.TagData.resize(tagDataLen);
		memcpy((char*)bAlignment.TagData.data(), tagData, tagDataLen);
	}

	free(allCharData);
	return true;
}

void BamReader::LoadReferenceData(void) {

	// get number of reference sequences
	char buffer[4];
	BgzfRead(buffer, 4);
	const unsigned int numberRefSeqs = BgzfUnpackUnsignedInt(buffer);
	if (numberRefSeqs == 0) { return; }
	m_references.reserve((int)numberRefSeqs);
	
	// iterate over all references in header
	for (unsigned int i = 0; i != numberRefSeqs; ++i) {

		// get length of reference name
		BgzfRead(buffer, 4);
		const unsigned int refNameLength = BgzfUnpackUnsignedInt(buffer);
		char* refName = (char*)calloc(refNameLength, 1);
		
		// get reference name and reference sequence length
		BgzfRead(refName, refNameLength);
		BgzfRead(buffer, 4);
		const unsigned int refLength = BgzfUnpackUnsignedInt(buffer);
		
		// store data for reference
		RefData aReference;
		aReference.RefName   = (string)((const char*)refName);
		aReference.RefLength = refLength;
		m_references.push_back(aReference);
		
		// clean up calloc-ed temp variable
		free(refName);
	}
}

// opens BAM file (and index)
void BamReader::Open(const string& filename, const string& indexFilename) {

	// open the BGZF file for reading, retrieve header text & reference data
	BgzfOpen(filename);
	LoadHeaderData();	
	LoadReferenceData();

	// store file offset of first alignment
	m_alignmentsBeginOffset = BgzfTell();

	// open index file & load index data (if exists)
	OpenIndex(indexFilename);
}

void BamReader::OpenIndex(const string& indexFilename) {

	// if index file exists
	if (!indexFilename.empty()) {

		// open index
		FILE* indexStream = fopen(indexFilename.c_str(), "rb");
		
		// abort on error
		if(!indexStream) {
			printf("ERROR: Unable to open the BAM index file %s for reading.\n", indexFilename.c_str() );
			exit(1);
		}
	
		// build up index data structure
		LoadIndexData(indexStream);
	}
}

bool BamReader::Rewind(void) {

	// find first reference that has alignments in the BAM file
	int refID = 0;
	int refCount = m_references.size();
	for ( ; refID < refCount; ++refID ) {
		if ( m_references.at(refID).RefHasAlignments ) { break; } 
	}

	// store default bounds for first alignment
	m_currentRefID = refID;
	m_currentLeft = 0;
	m_isRegionSpecified = false;

	// return success/failure of seek
	return BgzfSeek(m_alignmentsBeginOffset);
}	

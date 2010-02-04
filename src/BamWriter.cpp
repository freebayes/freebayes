// ***************************************************************************
// BamWriter.cpp (c) 2009 Michael Strömberg, Derek Barnett
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ---------------------------------------------------------------------------
// Last modified: 24 June 2009 (DB)
// ---------------------------------------------------------------------------
// The BGZF routines were adapted from the bgzf.c code developed at the Broad
// Institute.
// ---------------------------------------------------------------------------
// Provides the basic functionality for producing BAM files
// ***************************************************************************

#include "BamWriter.h"
using namespace BamTools;
using namespace std;

// constructor
BamWriter::BamWriter(void)
{}

// destructor
BamWriter::~BamWriter(void) {
	if(mBGZF.IsOpen) BgzfClose();
}

// closes the BAM file
void BamWriter::BgzfClose(void) {

	mBGZF.IsOpen = false;

	// flush the BGZF block
	BgzfFlushBlock();

	// flush and close
	fflush(mBGZF.Stream);
	fclose(mBGZF.Stream);
}

// compresses the current block
int BamWriter::BgzfDeflateBlock(void) {

	// initialize the gzip header
	char* buffer = mBGZF.CompressedBlock;
	unsigned int bufferSize = mBGZF.CompressedBlockSize;

	memset(buffer, 0, 18);
	buffer[0]  = GZIP_ID1;
	buffer[1]  = (char)GZIP_ID2;
	buffer[2]  = CM_DEFLATE;
	buffer[3]  = FLG_FEXTRA;
	buffer[9]  = (char)OS_UNKNOWN;
	buffer[10] = BGZF_XLEN;
	buffer[12] = BGZF_ID1;
	buffer[13] = BGZF_ID2;
	buffer[14] = BGZF_LEN;

	// loop to retry for blocks that do not compress enough
	int inputLength = mBGZF.BlockOffset;
	int compressedLength = 0;

	while(true) {

		z_stream zs;
		zs.zalloc    = NULL;
		zs.zfree     = NULL;
		zs.next_in   = (Bytef*)mBGZF.UncompressedBlock;
		zs.avail_in  = inputLength;
		zs.next_out  = (Bytef*)&buffer[BLOCK_HEADER_LENGTH];
		zs.avail_out = bufferSize - BLOCK_HEADER_LENGTH - BLOCK_FOOTER_LENGTH;

		// initialize the zlib compression algorithm
		if(deflateInit2(&zs, Z_DEFAULT_COMPRESSION, Z_DEFLATED, GZIP_WINDOW_BITS, Z_DEFAULT_MEM_LEVEL, Z_DEFAULT_STRATEGY) != Z_OK) {
			printf("ERROR: zlib deflate initialization failed.\n");
			exit(1);
		}

		// compress the data
		int status = deflate(&zs, Z_FINISH);
		if(status != Z_STREAM_END) {
			deflateEnd(&zs);

			// reduce the input length and try again
			if(status == Z_OK) {
				inputLength -= 1024;
				if(inputLength < 0) {
					printf("ERROR: input reduction failed.\n");
					exit(1);
				}
				continue;
			}

			printf("ERROR: zlib deflate failed.\n");
			exit(1);
		}

		// finalize the compression routine
		if(deflateEnd(&zs) != Z_OK) {
			printf("ERROR: deflate end failed.\n");
			exit(1);
		}

		compressedLength = zs.total_out;
		compressedLength += BLOCK_HEADER_LENGTH + BLOCK_FOOTER_LENGTH;

		if(compressedLength > MAX_BLOCK_SIZE) {
			printf("ERROR: deflate overflow.\n");
			exit(1);
		}

		break;
	}

	// store the compressed length
	BgzfPackUnsignedShort(&buffer[16], (unsigned short)(compressedLength - 1));

	// store the CRC32 checksum
	unsigned int crc = crc32(0, NULL, 0);
	crc = crc32(crc, (Bytef*)mBGZF.UncompressedBlock, inputLength);
	BgzfPackUnsignedInt(&buffer[compressedLength - 8], crc);
	BgzfPackUnsignedInt(&buffer[compressedLength - 4], inputLength);

	// ensure that we have less than a block of data left
	int remaining = mBGZF.BlockOffset - inputLength;
	if(remaining > 0) {
		if(remaining > inputLength) {
			printf("ERROR: remainder too large.\n");
			exit(1);
		}

		memcpy(mBGZF.UncompressedBlock, mBGZF.UncompressedBlock + inputLength, remaining);
	}

	mBGZF.BlockOffset = remaining;
	return compressedLength;
}

// flushes the data in the BGZF block
void BamWriter::BgzfFlushBlock(void) {

	// flush all of the remaining blocks
	while(mBGZF.BlockOffset > 0) {

		// compress the data block
		int blockLength = BgzfDeflateBlock();

		// flush the data to our output stream
		int numBytesWritten = fwrite(mBGZF.CompressedBlock, 1, blockLength, mBGZF.Stream);

		if(numBytesWritten != blockLength) {
			printf("ERROR: Expected to write %u bytes during flushing, but wrote %u bytes.\n", blockLength, numBytesWritten);
			exit(1);
		}

		mBGZF.BlockAddress += blockLength;
	}
}

// opens the BAM file for writing
void BamWriter::BgzfOpen(const string& filename) {

	mBGZF.Stream = fopen(filename.c_str(), "wb");

	if(!mBGZF.Stream) {
		printf("ERROR: Unable to open the BAM file (%s) for writing.\n", filename.c_str());
		exit(1);
	}

	mBGZF.IsOpen = true;
}

// writes the supplied data into the BGZF buffer
unsigned int BamWriter::BgzfWrite(const char* data, const unsigned int dataLen) {

	// initialize
	unsigned int numBytesWritten = 0;
	const char* input = data;
	unsigned int blockLength = mBGZF.UncompressedBlockSize;

	// copy the data to the buffer
	while(numBytesWritten < dataLen) {
		unsigned int copyLength = min(blockLength - mBGZF.BlockOffset, dataLen - numBytesWritten);
		char* buffer = mBGZF.UncompressedBlock;
		memcpy(buffer + mBGZF.BlockOffset, input, copyLength);

		mBGZF.BlockOffset += copyLength;
		input             += copyLength;
		numBytesWritten   += copyLength;

		if(mBGZF.BlockOffset == blockLength) BgzfFlushBlock();
	}

	return numBytesWritten;
}

// closes the alignment archive
void BamWriter::Close(void) {
	if(mBGZF.IsOpen) BgzfClose();
}

// creates a cigar string from the supplied alignment
void BamWriter::CreatePackedCigar(const vector<CigarOp>& cigarOperations, string& packedCigar) {

	// initialize
	const unsigned int numCigarOperations = cigarOperations.size();
	packedCigar.resize(numCigarOperations * SIZEOF_INT);

	// pack the cigar data into the string
	unsigned int* pPackedCigar = (unsigned int*)packedCigar.data();

	unsigned int cigarOp;
	vector<CigarOp>::const_iterator coIter;
	for(coIter = cigarOperations.begin(); coIter != cigarOperations.end(); coIter++) {

		switch(coIter->Type) {
	  case 'M':
		  cigarOp = BAM_CMATCH;
		  break;
	  case 'I':
		  cigarOp = BAM_CINS;
		  break;
	  case 'D':
		  cigarOp = BAM_CDEL;
		  break;
	  case 'N':
		  cigarOp = BAM_CREF_SKIP;
		  break;
	  case 'S':
		  cigarOp = BAM_CSOFT_CLIP;
		  break;
	  case 'H':
		  cigarOp = BAM_CHARD_CLIP;
		  break;
	  case 'P':
		  cigarOp = BAM_CPAD;
		  break;
	  default:
		  printf("ERROR: Unknown cigar operation found: %c\n", coIter->Type);
		  exit(1);
		}

		*pPackedCigar = coIter->Length << BAM_CIGAR_SHIFT | cigarOp;
		pPackedCigar++;
	}
}

// encodes the supplied query sequence into 4-bit notation
void BamWriter::EncodeQuerySequence(const string& query, string& encodedQuery) {

	// prepare the encoded query string
	const unsigned int queryLen = query.size();
	const unsigned int encodedQueryLen = (unsigned int)((queryLen / 2.0) + 0.5);
	encodedQuery.resize(encodedQueryLen);
	char* pEncodedQuery = (char*)encodedQuery.data();
	const char* pQuery = (const char*)query.data();

	unsigned char nucleotideCode;
	bool useHighWord = true;

	while(*pQuery) {

		switch(*pQuery) {
			case '=':
				nucleotideCode = 0;
				break;
			case 'A':
				nucleotideCode = 1;
				break;
			case 'C':
				nucleotideCode = 2;
				break;
			case 'G':
				nucleotideCode = 4;
				break;
			case 'T':
				nucleotideCode = 8;
				break;
			case 'N':
				nucleotideCode = 15;
				break;
			default:
				printf("ERROR: Only the following bases are supported in the BAM format: {=, A, C, G, T, N}. Found [%c]\n", *pQuery);
				exit(1);
		}

		// pack the nucleotide code
		if(useHighWord) {
			*pEncodedQuery = nucleotideCode << 4;
			useHighWord = false;
		} else {
			*pEncodedQuery |= nucleotideCode;
			pEncodedQuery++;
			useHighWord = true;
		}

		// increment the query position
		pQuery++;
	}
}

// opens the alignment archive
void BamWriter::Open(const string& filename, const string& samHeader, const RefVector& referenceSequences) {

	// open the BGZF file for writing
	BgzfOpen(filename);

	// ================
	// write the header
	// ================

	// write the BAM signature
	const unsigned char SIGNATURE_LENGTH = 4;
	const char* BAM_SIGNATURE = "BAM\1";
	BgzfWrite(BAM_SIGNATURE, SIGNATURE_LENGTH);

	// write the SAM header text length
	const unsigned int samHeaderLen = samHeader.size();
	BgzfWrite((char*)&samHeaderLen, SIZEOF_INT);

	// write the SAM header text
	if(samHeaderLen > 0) BgzfWrite(samHeader.data(), samHeaderLen);

	// write the number of reference sequences
	const unsigned int numReferenceSequences = referenceSequences.size();
	BgzfWrite((char*)&numReferenceSequences, SIZEOF_INT);

	// =============================
	// write the sequence dictionary
	// =============================

	RefVector::const_iterator rsIter;
	for(rsIter = referenceSequences.begin(); rsIter != referenceSequences.end(); rsIter++) {

		// write the reference sequence name length
		const unsigned int referenceSequenceNameLen = rsIter->RefName.size() + 1;
		BgzfWrite((char*)&referenceSequenceNameLen, SIZEOF_INT);

		// write the reference sequence name
		BgzfWrite(rsIter->RefName.c_str(), referenceSequenceNameLen);

		// write the reference sequence length
		BgzfWrite((char*)&rsIter->RefLength, SIZEOF_INT);
	}
}

// saves the alignment to the alignment archive
void BamWriter::SaveAlignment(const BamAlignment& al) {

	// initialize
	const unsigned int nameLen            = al.Name.size() + 1;
	const unsigned int queryLen           = al.QueryBases.size();
	const unsigned int numCigarOperations = al.CigarData.size();

	// create our packed cigar string
	string packedCigar;
	CreatePackedCigar(al.CigarData, packedCigar);
	const unsigned int packedCigarLen = packedCigar.size();

	// encode the query
	string encodedQuery;
	EncodeQuerySequence(al.QueryBases, encodedQuery);
	const unsigned int encodedQueryLen = encodedQuery.size();

	// store the tag data length
	const unsigned int tagDataLength = al.TagData.size() + 1;

	// assign the BAM core data
	unsigned int buffer[8];
	buffer[0] = al.RefID;
	buffer[1] = al.Position;
	buffer[2] = (al.Bin << 16) | (al.MapQuality << 8) | nameLen;
	buffer[3] = (al.AlignmentFlag << 16) | numCigarOperations;
	buffer[4] = queryLen;
	buffer[5] = al.MateRefID;
	buffer[6] = al.MatePosition;
	buffer[7] = al.InsertSize;

	// write the block size
	const unsigned int dataBlockSize = nameLen + packedCigarLen + encodedQueryLen + queryLen + tagDataLength;
	const unsigned int blockSize = BAM_CORE_SIZE + dataBlockSize;

	BgzfWrite((char*)&blockSize, SIZEOF_INT);

	// write the BAM core
	BgzfWrite((char*)&buffer, BAM_CORE_SIZE);

	// write the query name
	BgzfWrite(al.Name.c_str(), nameLen);

	// write the packed cigar
	BgzfWrite(packedCigar.data(), packedCigarLen);

	// write the encoded query sequence
	BgzfWrite(encodedQuery.data(), encodedQueryLen);

	// write the base qualities
	string baseQualities = al.Qualities;
	char* pBaseQualities = (char*)al.Qualities.data();
	for(unsigned int i = 0; i < queryLen; i++) pBaseQualities[i] -= 33;
	BgzfWrite(pBaseQualities, queryLen);

	// write the read group tag
	BgzfWrite(al.TagData.data(), tagDataLength);
}

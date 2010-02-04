// ***************************************************************************
// BamAux.h (c) 2009 Derek Barnett, Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ---------------------------------------------------------------------------
// Last modified: 24 June 2009 (DB)
// ---------------------------------------------------------------------------
// The BGZF routines were adapted from the bgzf.c code developed at the Broad
// Institute.
// ---------------------------------------------------------------------------
// Defines common constants, typedefs, & data structures for BamTools.
// ***************************************************************************

/*! 
	\file BamAux.h
	\brief BamTools constants, typedefs, & data structures
*/

#pragma once

// C++ includes
#include <exception>
#include <string>
#include <utility>
#include <vector>

// C includes
#include <cstdio>
#include <cstdlib>
#include <cstring>

// Platform-specific type definitions
#ifdef WIN32
	typedef char                 int8_t;
	typedef unsigned char       uint8_t;
	typedef short               int16_t;
	typedef unsigned short     uint16_t;
	typedef int                 int32_t;
	typedef unsigned int       uint32_t;
	typedef long long           int64_t;
	typedef unsigned long long uint64_t;
#else
	#include <stdint.h>
#endif

//! \namespace BamTools
namespace BamTools {

	//! \cond
	// --------------------------------------------------------------------------------------
	// This section is purely internal and can be excluded from main generated documentation.
	
	// zlib constants
	const int GZIP_ID1   = 31;
	const int GZIP_ID2   = 139;
	const int CM_DEFLATE = 8;
	const int FLG_FEXTRA = 4;
	const int OS_UNKNOWN = 255;
	const int BGZF_XLEN  = 6;
	const int BGZF_ID1   = 66;
	const int BGZF_ID2   = 67;
	const int BGZF_LEN   = 2;
	const int GZIP_WINDOW_BITS = -15;
	const int Z_DEFAULT_MEM_LEVEL = 8;

	// BZGF constants
	const int BLOCK_HEADER_LENGTH = 18;
	const int BLOCK_FOOTER_LENGTH = 8;
	const int MAX_BLOCK_SIZE      = 65536;
	const int DEFAULT_BLOCK_SIZE  = 65536;

	// BAM constants
	const unsigned int BAM_CORE_SIZE = 32;
	const int BAM_CMATCH      = 0;
	const int BAM_CINS        = 1;
	const int BAM_CDEL        = 2;
	const int BAM_CREF_SKIP   = 3;
	const int BAM_CSOFT_CLIP  = 4;
	const int BAM_CHARD_CLIP  = 5;
	const int BAM_CPAD        = 6;
	const int BAM_CIGAR_SHIFT = 4;
	const int BAM_CIGAR_MASK  = ((1 << BAM_CIGAR_SHIFT) - 1);

	// BAM index constants
	const int MAX_BIN           = 37450;	// =(8^6-1)/7+1
	const int BAM_MIN_CHUNK_GAP = 32768;
	const int BAM_LIDX_SHIFT    = 14;

	// Explicit variable sizes
	const int SIZEOF_INT = 4;
	
	struct BgzfData {
		unsigned int UncompressedBlockSize;
		unsigned int CompressedBlockSize;
		unsigned int BlockLength;
		unsigned int BlockOffset;
		uint64_t BlockAddress;
		bool     IsOpen;
		FILE*    Stream;
		char*    UncompressedBlock;
		char*    CompressedBlock;
		
		// constructor
		BgzfData(void)
			: UncompressedBlockSize(DEFAULT_BLOCK_SIZE)
			, CompressedBlockSize(MAX_BLOCK_SIZE)
			, BlockLength(0)
			, BlockOffset(0)
			, BlockAddress(0)
			, IsOpen(false)
			, Stream(NULL)
			, UncompressedBlock(NULL)
			, CompressedBlock(NULL)
		{
			try {
				CompressedBlock   = new char[CompressedBlockSize];
				UncompressedBlock = new char[UncompressedBlockSize];
			} catch( std::bad_alloc& ba ) {
				printf("ERROR: Unable to allocate memory for our BGZF object.\n");
				exit(1);
			}
		}
		
		// destructor
		~BgzfData(void) {
			if(CompressedBlock)   delete [] CompressedBlock;
			if(UncompressedBlock) delete [] UncompressedBlock;
		}
	};
	//! \endcond
	
	// --------------------------------------------------------------------------------------
	// Data structures
	
	//! \brief Cigar operation data structure
	struct CigarOp {
		char     Type;   //!< Operation type (MIDNSHP)
		uint32_t Length; //!< Operation length (number of bases)
		   
	};

	//! Reference sequence data structure
	struct RefData {
		std::string  RefName;          //!< Name of reference sequence
		unsigned int RefLength;        //!< Length of reference sequence
		bool         RefHasAlignments; //!< True if BAM file contains alignments mapped to reference sequence
		
		// constructor
		RefData(void)
			: RefLength(0)
			, RefHasAlignments(false)
		{ }
	};

	//! BAM alignment data structure
	struct BamAlignment {
		
		// Queries against alignment flag
		public:
			//! Returns true if this read is a PCR duplicate (determined by external app)
			bool IsDuplicate(void) const { return ( (AlignmentFlag & DUPLICATE) != 0 ); }
			//! Returns true if this read failed quality control (determined by external app)
			bool IsFailedQC(void) const { return ( (AlignmentFlag & QC_FAILED) != 0 ); }	
			//! Returns true if alignment is first mate on read
			bool IsFirstMate(void) const { return ( (AlignmentFlag & READ_1) != 0 ); }
			//! Returns true if alignment is mapped			
			bool IsMapped(void) const { return ( (AlignmentFlag & UNMAPPED) == 0 ); }
			//! Returns true if alignment's mate is mapped
			bool IsMateMapped(void) const { return ( (AlignmentFlag & MATE_UNMAPPED) == 0 ); }	
			//! Returns true if alignment's mate mapped to reverse strand
			bool IsMateReverseStrand(void) const { return ( (AlignmentFlag & MATE_REVERSE)  != 0 ); }
			//! Returns true if alignment part of paired-end read
			bool IsPaired(void) const { return ( (AlignmentFlag & PAIRED) != 0 ); }	
			//! Returns true if this position is primary alignment (determined by external app)
			bool IsPrimaryAlignment(void) const  { return ( (AlignmentFlag & SECONDARY) == 0 ); }	
			//! Returns true if alignment is part of read that satisfied paired-end resolution (determined by external app)
			bool IsProperPair(void) const { return ( (AlignmentFlag & PROPER_PAIR) != 0 ); }	
			//! Returns true if alignment mapped to reverse strand
			bool IsReverseStrand(void) const { return ( (AlignmentFlag & REVERSE) != 0 ); }	
			//! Returns true if alignment is second mate on read
			bool IsSecondMate(void) const { return ( (AlignmentFlag & READ_2) != 0 ); }	
		
		public:
			/*! 
				\brief Get alignment's read group text.
				
				Assigns read group text to readGroup.
				
				\return True if read group data successfully retrieved.
			*/
			bool GetReadGroup(std::string& readGroup) const {
				
				if ( TagData.empty() ) { return false; }
				
				// localize the tag data
				char* pTagData = (char*)TagData.data();
				const unsigned int tagDataLen = TagData.size();
				unsigned int numBytesParsed = 0;
				
				bool foundReadGroupTag = false;
				while( numBytesParsed < tagDataLen ) {
					
					const char* pTagType = pTagData;
					const char* pTagStorageType = pTagData + 2;
					pTagData       += 3;
					numBytesParsed += 3;
					
					// check the current tag
					if ( strncmp(pTagType, "RG", 2) == 0 ) {
						foundReadGroupTag = true;
						break;
					}
					
					// get the storage class and find the next tag
					SkipToNextTag( *pTagStorageType, pTagData, numBytesParsed );
				}
				
				// return if the read group tag was not present
				if ( !foundReadGroupTag ) { return false; }
				
				// assign the read group
				const unsigned int readGroupLen = strlen(pTagData);
				readGroup.resize(readGroupLen);
				memcpy( (char*)readGroup.data(), pTagData, readGroupLen );
				return true;
			}
		
		private:
			// skips to the next tag
			static void SkipToNextTag(const char storageType, char* &pTagData, unsigned int& numBytesParsed) {
				switch(storageType) {
					
					case 'A':
					case 'c':
					case 'C':
						++numBytesParsed;
						++pTagData;
						break;
						
					case 's':
					case 'S':
					case 'f':
						numBytesParsed += 2;
						pTagData       += 2;
						break;
						
					case 'i':
					case 'I':
						numBytesParsed += 4;
						pTagData       += 4;
						break;
						
					case 'Z':
					case 'H':
						while(*pTagData) {
							++numBytesParsed;
							++pTagData;
						}
						break;
						
					default:
						printf("ERROR: Unknown tag storage class encountered: [%c]\n", *pTagData);
						exit(1);
				}
			}
		
		// Data members
		public:
			std::string  Name;              //!< Read name
			unsigned int Length;            //!< Query length
			std::string  QueryBases;        //!< 'Original' sequence (as reported from sequencing machine)
			std::string  AlignedBases;      //!< 'Aligned' sequence (includes any indels, padding, clipping) 
			std::string  Qualities;         //!< FASTQ qualities (ASCII characters, not numeric values)
			std::string  TagData;           //!< Tag data (accessor methods will pull the requested information out)
			unsigned int RefID;             //!< ID number for reference sequence
			unsigned int Position;          //!< Position (0-based) where alignment starts 
			unsigned int Bin;               //!< Bin in BAM file where this alignment resides
			unsigned int MapQuality;        //!< Mapping quality score 
			unsigned int AlignmentFlag;     //!< Alignment bit-flag - see Is<something>() methods for available queries
			std::vector<CigarOp> CigarData; //!< CIGAR operations for this alignment
			unsigned int MateRefID;         //!< ID number for reference sequence where alignment's mate was aligned
			unsigned int MatePosition;      //!< Position (0-based) where alignment's mate starts
			unsigned int InsertSize;        //!< Mate-pair insert size
		
		// Alignment flag query constants
		private:
			enum { PAIRED        = 1,
				   PROPER_PAIR   = 2,
				   UNMAPPED      = 4,
				   MATE_UNMAPPED = 8,
				   REVERSE       = 16,
				   MATE_REVERSE  = 32,
				   READ_1        = 64,
				   READ_2        = 128,
				   SECONDARY     = 256,
				   QC_FAILED     = 512,
				   DUPLICATE     = 1024
			     };
	};

	// ----------------------------------------------------------------
	// Typedefs
	
	/*!
		\typedef RefVector
		\brief Vector of RefData objects
	*/
	typedef std::vector<RefData> RefVector;
	
	/*! 
		\typedef BamAlignmentVector
		\brief Vector of BamAlignments
	*/
	typedef std::vector< BamAlignment > BamAlignmentVector;
	
	//! \cond
	// ----------------------------------------------------------------
	// Typedefs (internal - can exclude from main documentation)
	
	//Offsets for linear indexing
	typedef std::vector<uint64_t> LinearOffsetVector;

	// Alignment 'chunk' boundaries
	typedef std::pair<uint64_t, uint64_t> ChunkPair;
	// Vector of alignment 'chunks'
	typedef std::vector<ChunkPair> ChunkVector;

	// BAM bin contains a bin ID & a vector of alignment 'chunks'
	typedef std::pair<uint32_t, ChunkVector*> BamBin;
	// Vector of BAM bins
	typedef std::vector<BamBin> BinVector;

	// Reference sequence index data
	typedef std::pair<BinVector*, LinearOffsetVector*> RefIndex;

	// Full BAM file index data structure 
	typedef std::vector<RefIndex*> BamIndex;
	// ----------------------------------------------------------------
	//! \endcond
}

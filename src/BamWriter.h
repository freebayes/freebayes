// ***************************************************************************
// BamWriter.h (c) 2009 Michael Strömberg, Derek Barnett
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

/*! 
	\file BamWriter.h
	\brief API for writing BAM files.
*/

#pragma once

// C++ includes
#include <string>
#include <vector>

// zlib includes
#include <zlib.h>

// BamTools includes
#include "BamAux.h"

namespace BamTools { 

	//! API for writing BAM files.
	class BamWriter {
		
		public:
			
			//! Constructor.
			BamWriter(void);
			
			//! Destructor.
			~BamWriter(void);
		
		public:
			
			/*! 
				\brief Closes the alignment archive
				\sa Open()
			*/
			void Close(void);
			
			/*! 
				\brief Opens the alignment archive
				\param filename output BAM file
				\param samHeader SAM-format header text
				\param referenceSequences Reference sequence data 
				\sa Close()
			*/
			void Open(const std::string& filename, const std::string& samHeader, const RefVector& referenceSequences);
			
			/*! 
				\brief Saves an alignment to the archive
				\param al The BamAlignment to be saved
			*/
			void SaveAlignment(const BamAlignment& al);
		
		// --------------------------------------------------------------------------------------
		// internal methods
		private:
			// closes the BAM file
			void BgzfClose(void);
			// compresses the current block
			int BgzfDeflateBlock(void);
			// flushes the data in the BGZF block
			void BgzfFlushBlock(void);
			// opens the BAM file for writing
			void BgzfOpen(const std::string& filename);
			// packs an unsigned integer into the specified buffer
			static inline void BgzfPackUnsignedInt(char* buffer, unsigned int value);
			// packs an unsigned short into the specified buffer
			static inline void BgzfPackUnsignedShort(char* buffer, unsigned short value);
			// writes the supplied data into the BGZF buffer
			unsigned int BgzfWrite(const char* data, const unsigned int dataLen);
			// calculates the minimum bin that contains a region [begin, end)
			static inline unsigned int CalculateMinimumBin(unsigned int begin, unsigned int end);
			// creates a packed cigar string from the supplied alignment
			static void CreatePackedCigar(const std::vector<CigarOp>& cigarOperations, std::string& packedCigar);
			// encodes the supplied query sequence into 4-bit notation
			static void EncodeQuerySequence(const std::string& query, std::string& encodedQuery);
			// our BGZF output object
			BgzfData mBGZF;
	};

	//! \cond
	// --------------------------------------------------------------------------------------
	// static inline methods (internal - can exclude from main documentation)
	
	// packs an unsigned integer into the specified buffer
	inline void BamWriter::BgzfPackUnsignedInt(char* buffer, unsigned int value) {
		buffer[0] = (char)value;
		buffer[1] = (char)(value >> 8);
		buffer[2] = (char)(value >> 16);
		buffer[3] = (char)(value >> 24);
	}

	// packs an unsigned short into the specified buffer
	inline void BamWriter::BgzfPackUnsignedShort(char* buffer, unsigned short value) {
		buffer[0] = (char)value;
		buffer[1] = (char)(value >> 8);
	}
	// --------------------------------------------------------------------------------------
	//! \endcond
	
} // end BamTools namespace
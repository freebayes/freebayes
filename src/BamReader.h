// ***************************************************************************
// BamReader.h (c) 2009 Derek Barnett, Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ---------------------------------------------------------------------------
// Last modified: 24 June 2009 (DB)
// ---------------------------------------------------------------------------
// The BGZF routines were adapted from the bgzf.c code developed at the Broad
// Institute.
// ---------------------------------------------------------------------------
// Provides the basic functionality for reading BAM files
// ***************************************************************************

/*! 
	\file BamReader.h
	\brief API for reading BAM files.
*/

#pragma once

// C++ includes
#include <algorithm>
#include <iterator>
#include <string>
#include <utility>
#include <vector>

// zlib includes
#include <zlib.h>

// BamTools includes
#include "BamAux.h"

namespace BamTools { 

	//! API for reading BAM files.
	class BamReader {
		
		public:
			
			//! Constructor
			BamReader(void);
			
			//! Destructor
			~BamReader(void);
		
		public:
			
			/*! 
				\brief Closes the BAM file.
				
				Also closes index file and clears index data, if provided.
				
				\sa Open()
			*/
			void Close(void);
			
			/*! 
				\brief Access SAM format header data.
				
				See SAM format documentation for detailed header description.
				
				\return Full header text (no parsing of tags)
			*/
			const std::string GetHeaderText(void) const;
			
			/*! 
				\brief Retrieve next alignment.
				
				Stores result in bAlignment.  
				
				If reference and position are specified by a prior call to Jump(), this method stores the next aligmment that either: 
				a) overlaps, or 
				b) begins on/after that specified position.
				
				Otherwise this simply stores next alignment, if one exists.
				
				Note that this method does not specifiy a 'right bound' position.  
				If a range is desired, you should supply some stopping criteria. For example:
				
				\code
				BamReader bReader;
				bReader.Open(bamFile, bamIndexfile);
				if ( bReader.Jump( someID, somePosition ) {
					BamAlignment bAlignment;
					while ( bReader.GetNextAlignment(bAlignment) && (bAlignment.Position <= someUpperBound) ) {
						// do something
					}
				}
				\endcode
				
				\param bAlignment destination for alignment data
				\return success/failure
				\sa Jump(), Rewind()
			*/
			bool GetNextAlignment(BamAlignment& bAlignment);
			
			/*! 
				\brief Get number of reference sequences in BAM file.
				\return Number of references
				\sa GetReferenceData(), GetReferenceID()
			*/
			const int GetReferenceCount(void) const;
			
			/*! 
				\brief Access reference data.
				\return Vector of RefData entry
				\sa GetReferenceCount(), GetReferenceID()
			*/
			const RefVector GetReferenceData(void) const;
			
			/*! 
				\brief Get reference ID from name.
				\param refName name of reference sequence
				\return reference ID number
				\sa GetReferenceCount(), GetReferenceData()
			*/
			const int GetReferenceID(const std::string& refName) const;
			
			/*! 
				\brief Random access in BAM file.
				
				Jump to a specified position on reference sequence. 
				Position is optional - defaults to beginning of specified reference.
				
				Reference and position are stored for use by subsequent calls to GetNextAlignment().
				
				\param refID ID number of desired reference
				\param position left-bound position  
				\return success/failure
				\sa GetNextAlignment(), Rewind()
			*/
			bool Jump(int refID, unsigned int position = 0);
			
			/*! 
				\brief Opens a BAM file.
				
				Index file is optional - sequential reading through a BAM file does not require an index.
				
				However, the index is required to perform random-access alignment retrival (using the Jump() method ).  
				
				See SAMtools documentation for BAM index generation.
				
				\param filename BAM file
				\param indexFilename BAM index file
				\sa Jump(), Close()
			*/
			void Open(const std::string& filename, const std::string& indexFilename = "");
			
			/*!
				\brief Moves file pointer to beginning of alignment data.
				
				A subsequent call to GetNextAlignment() would retrieve the first alignment in the BAM file.
				Clears any reference and position set by a prior call to Jump()
				
				\return success/failure
				\sa GetNextAlignment(), Jump()
			*/
			bool Rewind(void);
		
		// --------------------------------------------------------------------------------------
		// internal methods
		private:
			// checks BGZF block header
			bool BgzfCheckBlockHeader(char* header);
			// closes the BAM file
			void BgzfClose(void);
			// de-compresses the current block
			int BgzfInflateBlock(int blockLength);
			// opens the BAM file for reading
			void BgzfOpen(const std::string& filename);
			// reads BGZF data into a byte buffer
			unsigned int BgzfRead(char* data, const unsigned int dataLen);
			// reads BGZF block
			int BgzfReadBlock(void);
			// seek to position in BAM file
			bool BgzfSeek(int64_t position);
			// get file position in BAM file
			int64_t BgzfTell(void);
			// unpacks a buffer into an unsigned int
			static inline unsigned int BgzfUnpackUnsignedInt(char* buffer);
			// unpacks a buffer into an unsigned short
			static inline unsigned short BgzfUnpackUnsignedShort(char* buffer);
			// calculate bins that overlap region ( left to reference end for now )
			int BinsFromRegion(int, unsigned int, uint16_t[MAX_BIN]);
			// calculates alignment end position based on starting position and provided CIGAR operations
			unsigned int CalculateAlignmentEnd(const unsigned int& position, const std::vector<CigarOp>& cigarData);
			// clear out (delete pointers in) index data structure
			void ClearIndex(void);
			// calculate file offset for first alignment chunk overlapping 'left'
			int64_t GetOffset(int refID, unsigned int left);
			// checks to see if alignment overlaps current region
			bool IsOverlap(BamAlignment& bAlignment);
			// retrieves header text from BAM file
			void LoadHeaderData(void);
			// builds BamIndex data structure from BAM index file
			void LoadIndexData(FILE* indexStream);
			// retrieves BAM alignment under file pointer
			bool LoadNextAlignment(BamAlignment& bAlignment);
			// builds reference data structure from BAM file
			void LoadReferenceData(void);
			// open BAM index file (if successful, loads index)
			void OpenIndex(const std::string& indexFilename);
		
		// aligment file & index file data
		private:
			BgzfData    m_BGZF;
			std::string m_headerText;
			BamIndex*   m_index;
			RefVector   m_references;
			bool        m_isIndexLoaded;
			int64_t     m_alignmentsBeginOffset;
		
		// user-specified region values
		private:
			bool         m_isRegionSpecified;
			int          m_currentRefID;
			unsigned int m_currentLeft;
		
		// BAM character constants
		private:
			static const char* DNA_LOOKUP;
			static const char* CIGAR_LOOKUP;
	};

	//! \cond
	// --------------------------------------------------------------------------------------
	// static inline methods (internal - can exclude from main documentation)
	
	// unpacks a buffer into an unsigned int
	inline unsigned int BamReader::BgzfUnpackUnsignedInt(char* buffer) {
		union { unsigned int value; unsigned char valueBuffer[sizeof(unsigned int)]; } un;
		un.valueBuffer[0] = buffer[0];
		un.valueBuffer[1] = buffer[1];
		un.valueBuffer[2] = buffer[2];
		un.valueBuffer[3] = buffer[3];
		return un.value;
	}

	// unpacks a buffer into an unsigned short
	inline unsigned short BamReader::BgzfUnpackUnsignedShort(char* buffer) {
		union { unsigned short value; unsigned char valueBuffer[sizeof(unsigned short)];} un;
		un.valueBuffer[0] = buffer[0];
		un.valueBuffer[1] = buffer[1];
		return un.value;
	}

	// --------------------------------------------------------------------------------------
	// template classes/methods (internal - can exclude from main documentation)
	
	// allows sorting/searching of a vector of pairs (instead of using maps)
	template <typename Key, typename Value>
	class LookupKeyCompare {

		typedef std::pair< Key, Value > LookupData;
		typedef typename LookupData::first_type Key_t;
		
		public:
			bool operator() (const LookupData& lhs, const LookupData& rhs) const { return keyLess(lhs.first, rhs.first); }
			bool operator() (const LookupData& lhs, const Key_t& k) const        { return keyLess(lhs.first, k); }
			bool operator() (const Key_t& k, const LookupData& rhs) const        { return keyLess(k, rhs.first); }
		private:
			bool keyLess(const Key_t& k1, const Key_t& k2) const { return k1 < k2; }
	};

	// return index of item if found, else return container.size()
	template < typename InputIterator, typename EqualityComparable >
	typename std::iterator_traits<InputIterator>::difference_type
	Index(const InputIterator& begin, const InputIterator& end, const EqualityComparable& item) {
		return std::distance(begin, std::find(begin, end, item));
	}
	//! \endcond
	
}

#pragma once

#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#ifndef WIN32
#include <tr1/unordered_map>
#else
#include <unordered_map>
#endif
#include <vector>
#include <math.h>
#include <string.h>
#include "fastlz.h"

using namespace std;
using namespace std::tr1;

// add some large file support
#ifdef WIN32
#define ftell64(a)     _ftelli64(a)
#define fseek64(a,b,c) _fseeki64(a,b,c)
typedef __int64 off_type;
#elif SPARC
#define ftell64(a)     ftello(a)
#define fseek64(a,b,c) fseeko(a,b,c)
typedef off64_t off_type;
#elif __APPLE__
#define ftell64(a)     ftello(a)
#define fseek64(a,b,c) fseeko(a,b,c)
#define fstat64(a,b)   fstat(a,b)
#define stat64 stat
typedef off_t off_type;
#else
#define ftell64(a)     ftello(a)
#define fseek64(a,b,c) fseeko(a,b,c)
typedef off_t off_type;
#endif

// define some data type sizes
//#define SIZEOF_CHAR          1
//#define SIZEOF_SHORT         2
#define SIZEOF_INT           4
#define SIZEOF_UINT64        8
#define SIZEOF_OFF_TYPE      8

typedef unsigned char ReferenceSequenceStatus;

#define REF_UNKNOWN    0 
#define REF_COLORSPACE 1

// the buffer is currently set to 10 MB
#define FASTLZ_IO_BUFFER_LEN        10485760

// the buffer must be at least 5 % larger
#define FASTLZ_IO_OUTPUT_BUFFER_LEN  9986438

#define FASTLZ_BETTER_SPEED       1
#define FASTLZ_BETTER_COMPRESSION 2

namespace Mosaik {

	// define our 64-bit data type
	typedef unsigned long long uint64_t;

	// define our reference sequence
	struct ReferenceSequence {
		off_type BasesOffset;
		uint64_t NumAligned;
		unsigned int Begin;
		unsigned int End;
		unsigned int NumBases;
		string Name;
		string Bases;
		string GenomeAssemblyID;
		string Species;
		string MD5;
		string URI;

		// constructor
		ReferenceSequence()
			: BasesOffset(0)
			, NumAligned(0)
			, Begin(0)
			, End(0)
			, NumBases(0)
		{}
	};

	class CFastLZIO {
	public:
		// constructor
		CFastLZIO(void);
		// destructor
		~CFastLZIO(void);
		// clears the buffer
		void Clear(void);
		// our input method
		void Read(char* &buffer, unsigned int& bufferLen, FILE* stm);
		// our input method (STL string)
		void Read(string& s, FILE* stm);
		// our output method
		void Write(const char* buffer, const unsigned int bufferLen, FILE* stm);
	private:
		// our buffer
		char* mBuffer;
		unsigned int mBufferLen;
	};

	class CReferenceSequenceReader {
	public:
		// constructor
		CReferenceSequenceReader(void);
		// destructor
		~CReferenceSequenceReader(void);
		// checks to see if this is truly a MOSAIK reference sequence archive
		static bool CheckFile(const string& filename, const bool showError);
		// closes the reference sequence archive
		void Close(void);
		// returns the reference sequence length
		unsigned int GetReferenceSequenceLength(void) const;
		// retrieves the desired reference sequence and places it in the specified string
		void GetReferenceSequence(const string& name, string& bases);
		// adds the reference sequences to the supplied vector
		void GetReferenceSequences(vector<ReferenceSequence>& referenceSequences);
		// initializes the supplied pointer with the concatenated reference sequence
		void LoadConcatenatedSequence(char* &referenceSequence);
		//initializes the supplied pointer with the 2-bit concatenated reference sequence
		void Load2BitConcatenatedSequence(char* &referenceSequence, char* &maskSequence, unsigned int& numMaskedPositions);
		// opens the reference sequence archive
		void Open(const string& filename);

	private:
		// define a comparison function for sorting our alignment positions (ascending)
		struct SortReferenceSequencesByBeginAsc {
			bool operator()(const ReferenceSequence& ar1, const ReferenceSequence& ar2) {
				return ar1.Begin < ar2.Begin;
			}
		};
		// stores the file state
		bool mIsOpen;
		// our input file stream
		FILE* mInStream;
		// our offsets
		off_type mConcatenatedOffset;
		off_type mConcatenated2bOffset;
		off_type mIndexOffset;
		off_type mReferenceBasesOffset;
		off_type mMaskedRegionsOffset;
		// the number of reference sequences contained in the anchor archive
		unsigned int mNumReferenceSequences;
		// the concatenated sequence length
		unsigned int mConcatenatedLen;
		// the concatenated 2-bit sequence length
		unsigned int mConcatenated2bLen;
		// our file status
		ReferenceSequenceStatus mStatus;
		// stores the index for our reference sequences
		unordered_map<string, ReferenceSequence> mIndex;
		// our FastLZ IO reader and writer
		CFastLZIO mFIO;
	};
}

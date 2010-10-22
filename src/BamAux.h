// ***************************************************************************
// BamAux.h (c) 2009 Derek Barnett, Michael Str�mberg
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ---------------------------------------------------------------------------
// Last modified: 9 October 2010 (DB)
// ---------------------------------------------------------------------------
// Provides the basic constants, data structures, utilities etc. 
// used throughout the API for handling BAM files
// ***************************************************************************

#ifndef BAMAUX_H
#define BAMAUX_H

#include <fstream> 
#include <iostream>
#include <string>
#include <vector>

// Platform-specific large-file support
#ifndef BAMTOOLS_LFS
#define BAMTOOLS_LFS
    #ifdef WIN32
        #define ftell64(a)     _ftelli64(a)
        #define fseek64(a,b,c) _fseeki64(a,b,c)
    #else
        #define ftell64(a)     ftello(a)
        #define fseek64(a,b,c) fseeko(a,b,c)
    #endif
#endif // BAMTOOLS_LFS

// Platform-specific type definitions
#ifndef BAMTOOLS_TYPES
#define BAMTOOLS_TYPES
    #ifdef _MSC_VER
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
#endif // BAMTOOLS_TYPES

namespace BamTools {

// ----------------------------------------------------------------
// ----------------------------------------------------------------
// BAM constants

const int BAM_CMATCH      = 0;
const int BAM_CINS        = 1;
const int BAM_CDEL        = 2;
const int BAM_CREF_SKIP   = 3;
const int BAM_CSOFT_CLIP  = 4;
const int BAM_CHARD_CLIP  = 5;
const int BAM_CPAD        = 6;
const int BAM_CIGAR_SHIFT = 4;
const int BAM_CIGAR_MASK  = ((1 << BAM_CIGAR_SHIFT) - 1);
const int BAM_CORE_SIZE   = 32;
const int BT_SIZEOF_INT   = 4;

// ----------------------------------------------------------------
// ----------------------------------------------------------------
// Data structs & typedefs

// CIGAR operation data structure
struct CigarOp {
  
    // data members
    char     Type;   // Operation type (MIDNSHP)
    uint32_t Length; // Operation length (number of bases)
    
    // constructor
    CigarOp(const char type = '\0', 
            const uint32_t length = 0) 
        : Type(type)
        , Length(length) 
    { }
};

// Reference data entry
struct RefData {
   
    // data members
    std::string RefName;          // Name of reference sequence
    int32_t     RefLength;        // Length of reference sequence
    bool        RefHasAlignments; // True if BAM file contains alignments mapped to reference sequence
    
    // constructor
    RefData(const int32_t& length = 0, 
            bool ok = false)
        : RefLength(length)
        , RefHasAlignments(ok)
    { }
};
typedef std::vector<RefData> RefVector;

// General (sequential) genome region
struct BamRegion {
  
    // data members
    int LeftRefID;
    int LeftPosition;
    int RightRefID;
    int RightPosition;
    
    // constructor
    BamRegion(const int& leftID   = -1, 
              const int& leftPos  = -1,
              const int& rightID  = -1,
              const int& rightPos = -1)
        : LeftRefID(leftID)
        , LeftPosition(leftPos)
        , RightRefID(rightID)
        , RightPosition(rightPos)
    { }
    
    // copy constructor
    BamRegion(const BamRegion& other)
	: LeftRefID(other.LeftRefID)
	, LeftPosition(other.LeftPosition)
	, RightRefID(other.RightRefID)
	, RightPosition(other.RightPosition)
    { }
    
    // member functions
    void clear(void) { LeftRefID = -1; LeftPosition = -1; RightRefID = -1; RightPosition = -1; }
    bool isLeftBoundSpecified(void) const { return ( LeftRefID != -1 && LeftPosition != -1 ); }
    bool isNull(void) const { return ( !isLeftBoundSpecified() && !isRightBoundSpecified() ); }
    bool isRightBoundSpecified(void) const { return ( RightRefID != -1 && RightPosition != -1 ); }
};

// ----------------------------------------------------------------
// ----------------------------------------------------------------
// General utilities 

// returns true if system is big endian
inline bool SystemIsBigEndian(void) {
   const uint16_t one = 0x0001;
   return ((*(char*) &one) == 0 );
}

// swaps endianness of 16-bit value 'in place'
inline void SwapEndian_16(int16_t& x) {
    x = ((x >> 8) | (x << 8));
}

inline void SwapEndian_16(uint16_t& x) {
    x = ((x >> 8) | (x << 8));
}

// swaps endianness of 32-bit value 'in-place'
inline void SwapEndian_32(int32_t& x) {
    x = ( (x >> 24) | 
         ((x << 8) & 0x00FF0000) | 
         ((x >> 8) & 0x0000FF00) | 
          (x << 24)
        );
}

inline void SwapEndian_32(uint32_t& x) {
    x = ( (x >> 24) | 
         ((x << 8) & 0x00FF0000) | 
         ((x >> 8) & 0x0000FF00) | 
          (x << 24)
        );
}

// swaps endianness of 64-bit value 'in-place'
inline void SwapEndian_64(int64_t& x) {
    x = ( (x >> 56) | 
         ((x << 40) & 0x00FF000000000000ll) |
         ((x << 24) & 0x0000FF0000000000ll) |
         ((x << 8)  & 0x000000FF00000000ll) |
         ((x >> 8)  & 0x00000000FF000000ll) |
         ((x >> 24) & 0x0000000000FF0000ll) |
         ((x >> 40) & 0x000000000000FF00ll) |
          (x << 56)
        );
}

inline void SwapEndian_64(uint64_t& x) {
    x = ( (x >> 56) | 
         ((x << 40) & 0x00FF000000000000ll) |
         ((x << 24) & 0x0000FF0000000000ll) |
         ((x << 8)  & 0x000000FF00000000ll) |
         ((x >> 8)  & 0x00000000FF000000ll) |
         ((x >> 24) & 0x0000000000FF0000ll) |
         ((x >> 40) & 0x000000000000FF00ll) |
          (x << 56)
        );
}

// swaps endianness of 'next 2 bytes' in a char buffer (in-place)
inline void SwapEndian_16p(char* data) {
    uint16_t& value = (uint16_t&)*data; 
    SwapEndian_16(value);
}

// swaps endianness of 'next 4 bytes' in a char buffer (in-place)
inline void SwapEndian_32p(char* data) {
    uint32_t& value = (uint32_t&)*data; 
    SwapEndian_32(value);
}

// swaps endianness of 'next 8 bytes' in a char buffer (in-place)
inline void SwapEndian_64p(char* data) {
    uint64_t& value = (uint64_t&)*data; 
    SwapEndian_64(value);
}

// returns whether file exists (can be opened OK)
inline bool FileExists(const std::string& filename) {
    std::ifstream f(filename.c_str(), std::ifstream::in);
    return !f.fail();
}

} // namespace BamTools

#endif // BAMAUX_H

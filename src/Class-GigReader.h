//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// Class-GigReader
// Class definition for GIG format binary assembly archive files
// Copyright 2008 Gabor T. Marth, Boston College
// All rights reserved
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

#ifndef CLASS_GIGREADER_H
#define CLASS_GIGREADER_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <deque>
#include <map>
#include <iterator>
#include <algorithm>

// "hash_map" true hashes
#include <ext/hash_map>

// tell g++ how to handle string or long long keys in hashmaps
namespace __gnu_cxx {
  template<> struct hash<std::string> {
    size_t operator()( const std::string& x ) const {
      return hash<const char*>()( x.c_str() );
    }
  };

  template<> struct hash<long long> {
    size_t operator()( const long long& x ) const {
      long long ret = (x >> 32L) ^ (x & 0xffffffff);
      return (size_t)ret;
    }
  };
}

using std::ios;
using std::ifstream;
using std::ofstream;
using std::ostream;
using std::istream;
using std::fstream;
using std::cin;
using std::cout;
using std::clog;
using std::endl;
using std::string;
using std::vector;
using std::deque;
using std::map;
using namespace __gnu_cxx;

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// utility routines
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// type definitions
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// Assembly data structure
//------------------------------------------------------------------------------
struct AssemblyData {
  string name;
  int numReads;
  int numContigs;
  vector<string> contigNames;
  bool valid;
};

//------------------------------------------------------------------------------
// Basesegment data structure
//------------------------------------------------------------------------------
struct BasesegmentData {
  string readName;
  int contigBegin;
  int contigEnd;
};

//------------------------------------------------------------------------------
// Contig data structure
//------------------------------------------------------------------------------
struct ContigData {
  string name;
  string dna;
  vector<short> qual;
  int length;
  int numReads;
  int numBasesegments;
  vector<string> readNames;
  vector<BasesegmentData> basesegments;
  bool valid;
};

//------------------------------------------------------------------------------
// Read data structure
//------------------------------------------------------------------------------
struct ReadData {
  string name;
  string dna;
  vector<short> qual;
  int length;
  int contigFitLeft;
  int contigFitRight;
  int readQualClipLeft;
  int readQualClipRight;
  int readAlignmentClipLeft;
  int readAlignmentClipRight;
  bool complemented;
  bool valid;
};

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// class definition
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------


class GigReader {

public:

  //----------------------------------------------------------------------------
  // public functions
  //----------------------------------------------------------------------------

  // constructor
  GigReader(
	    FILE *
	   );
  
  // getAssBlockOffset
  off64_t getAssBlockOffset();

  // getConBlockOffset
  off64_t getConBlockOffset();

  // getConIndexOffset
  off64_t getConIndexOffset();
  
  // getAssemblyData
  AssemblyData getAssemblyData();

  // getContigData
  ContigData getContigData(string);

  // getReadData
  ReadData getReadData(string, string);

  // prepareReads
  bool prepareReads(string);

  // getNextRead
  bool getNextRead(ReadData &);

  // getReads
  unsigned int getReads(int, deque<ReadData> &);
  
  //----------------------------------------------------------------------------
  // public variables
  //----------------------------------------------------------------------------

private:

  //----------------------------------------------------------------------------
  // private variables
  //----------------------------------------------------------------------------

  // file handle
  FILE * gig;

  // main file offsets
  off64_t conBlockOffset;
  off64_t conIndexOffset;
  off64_t assBlockOffset;

  // archive data characteristics
  bool haveQual;
  bool haveBpos;

  // assembly data
  AssemblyData assemblyData;

  // assembly name
  string assName;

  // number of contigs in archive
  int numCons;

  // list of contigs
  vector<string> conNames;

  // number of assembled reads in archive
  int assNumReads;

  // contig-specific assembled read lists
  hash_map<string, vector<string> > conAssReadNames;

  // hash_map of contig-specific assembled read indexes
  hash_map<string, off64_t> conOffsetHash;
  hash_map<string, hash_map<string, off64_t> > conAssReadOffset;
};

#endif


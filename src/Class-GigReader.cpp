//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// Class-GigReader
// Class code for GigReader
// Copyright 2008 Gabor T. Marth, Boston College
// All rights reserved
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

#ifndef CLASS_GIGREADER_CPP
#define CLASS_GIGREADER_CPP

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <deque>
#include <map>
#include <iterator>
#include <algorithm>
#include <cmath>

#include "Function-Generic.h"
#include "Function-Sequence.h"
#include "Class-GigReader.h"
#include "LargeFileSupport.h"

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
using std::min;
using std::max;
using namespace __gnu_cxx;

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// BasReader object class utility routines
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// constants
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// static class-wide variable initializations
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// class methods
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// constructor
//------------------------------------------------------------------------------
GigReader::GigReader(
		     FILE * gigIn
		     ) {

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // input
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  gig = gigIn;

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // initializations
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------

  // store current position
  off64_t currentOffset = ftell64(gig);

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // read archive id
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------

  // read length of archive id
  int idLength;
  fread((char*)&idLength, sizeof(int), 1, gig);

  char archiveIdCstyle[idLength];
  fread(archiveIdCstyle, sizeof(char), idLength+1, gig);
  string archiveId = archiveIdCstyle;

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // read assembly info
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------

  // assembly block offset
  assBlockOffset = ftell64(gig);

  // position to start of assembly block
  fseek64(gig, assBlockOffset, SEEK_SET);

  // read number of contigs
  int numberAssemblyContigs;
  fread((char*)&numberAssemblyContigs, sizeof(int), 1, gig);

  // read offset of contig index
  fread((char*)&conIndexOffset, sizeof(off64_t), 1, gig);

  // read offset of contig block
  fread((char*)&conBlockOffset, sizeof(off64_t), 1, gig);

  // read number of reads
  int numberAssemblyReads;
  fread((char*)&numberAssemblyReads, sizeof(int), 1, gig);

  // read length of assemblyName
  int assemblyNameCstyleLength; 
  fread((char*)&assemblyNameCstyleLength, sizeof(int), 1, gig);
  
  // read assemblyName
  char assemblyNameCstyle[assemblyNameCstyleLength];
  fread(assemblyNameCstyle, sizeof(char), assemblyNameCstyleLength+1, gig);

  // assign assembly name
  assName = assemblyNameCstyle;

  // assign number of reads
  assNumReads = numberAssemblyReads;

  // assign assembly info
  assemblyData.name = assName;
  assemblyData.numReads = numberAssemblyReads;
  assemblyData.numContigs = numberAssemblyContigs;

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // read contig index
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------

  // position to start of contig index
  fseek64(gig, conIndexOffset, SEEK_SET);

  //----------------------------------------------------------------------------
  // iterate through every contig in index, and register offset
  //----------------------------------------------------------------------------
  int conCount = 0;
  for (int i=0; i<numberAssemblyContigs; i++) {

    // increment conCount
    conCount++;

    // read length of conName
    int conNameCstyleLength; 
    fread((char*)&conNameCstyleLength, sizeof(int), 1, gig);

    // read conName
    char conNameCstyle[conNameCstyleLength];
    fread(conNameCstyle, sizeof(char), conNameCstyleLength+1, gig);
  
    // read conOffset
    off64_t conOffset; 
    fread((char*)&conOffset, sizeof(off64_t), 1, gig);

    // make conName string
    string conName(conNameCstyle);

    // register offset
    conOffsetHash[conName] = conOffset;

    // register contig
    conNames.push_back(conName);
  }

  assemblyData.contigNames = conNames;

  // restore current position
  fseek64(gig, currentOffset, SEEK_SET);
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// member methods
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// getAssBlockOffset
//------------------------------------------------------------------------------
off64_t GigReader::getAssBlockOffset() {
  return assBlockOffset;
}

//------------------------------------------------------------------------------
// getConBlockOffset
//------------------------------------------------------------------------------
off64_t GigReader::getConBlockOffset() {
  return conBlockOffset;
}

//------------------------------------------------------------------------------
// getConIndexOffset
//------------------------------------------------------------------------------
off64_t GigReader::getConIndexOffset() {
  return conIndexOffset;
}

//------------------------------------------------------------------------------
// getAssemblyData
//------------------------------------------------------------------------------
AssemblyData GigReader::getAssemblyData() {
  return assemblyData;
}

//------------------------------------------------------------------------------
// getContigData
//------------------------------------------------------------------------------
ContigData GigReader::getContigData(string conName) {

  //----------------------------------------------------------------------------
  // make ContigData output
  //----------------------------------------------------------------------------
  ContigData contigData;
  contigData.valid = false;

  //----------------------------------------------------------------------------
  // look up offset for this contig
  //----------------------------------------------------------------------------
  if (! (conOffsetHash.count(conName) > 0)) {
    return contigData;
  }
  off64_t conOffset = conOffsetHash[conName];

  //----------------------------------------------------------------------------
  // store current offset
  //----------------------------------------------------------------------------
  off64_t currentOffset = ftell64(gig);

  //----------------------------------------------------------------------------
  // read contig info from contigOffset location
  //----------------------------------------------------------------------------

  // position to conOffset
  fseek64(gig, conOffset, SEEK_SET);
  
  // read number of assembled reads in contig
  int conReads; 
  fread((char*)&conReads, sizeof(int), 1, gig);

  // read offset for read index
  off64_t readIndexOffset;
  fread((char*)&readIndexOffset, sizeof(off64_t), 1, gig);

  // read offset for read block
  off64_t readBlockOffset;
  fread((char*)&readBlockOffset, sizeof(off64_t), 1, gig);

  // read length of conName (but we don't need it)
  int conNameLength; 
  fread((char*)&conNameLength, sizeof(int), 1, gig);

  // read con length
  int conLength; 
  fread((char*)&conLength, sizeof(int), 1, gig);

  // read number of base segments in contig
  int conBasesegments; 
  fread((char*)&conBasesegments, sizeof(int), 1, gig);

  // read contigName
  char conNameCstyle[conNameLength];
  fread(conNameCstyle, sizeof(char), conNameLength+1, gig);

  // read contigDna: cannot use char array because of index limitation. use string directly
  //  char conDnaCstyle[conLength];
  //  fread(conDnaCstyle, sizeof(char), conLength+1, gig);
  string conDna;
  for (int i=0; i<conLength; i++) {
    char b;
    fread(&b, sizeof(char), 1, gig);
    conDna += b;
  }
  // read last character (C style terminator) and ignore it
  char b;
  fread(&b, sizeof(char), 1, gig);

  // read base quality values
  vector<short> conQual;
  for (int i=0; i<conLength; i++) {
    short q;
    fread((char*)&q, sizeof(short), 1, gig);
    conQual.push_back(q);
  }

  // read base segments
  vector<BasesegmentData> conBasesegmentList;
  for (int i=0; i<conBasesegments; i++) {

    // read length of readName
    int readNameLength; 
    fread((char*)&readNameLength, sizeof(int), 1, gig);

    // read contigBegin position
    int conBegin; 
    fread((char*)&conBegin, sizeof(int), 1, gig);

    // read contigEnd position
    int conEnd; 
    fread((char*)&conEnd, sizeof(int), 1, gig);

    // read readName
    char readNameCstyle[readNameLength];
    fread(readNameCstyle, sizeof(char), readNameLength+1, gig);
  
    // make BasesegmentData
    BasesegmentData bs;
    bs.readName = readNameCstyle;
    bs.contigBegin = conBegin;
    bs.contigEnd = conEnd;

    // add base segment to contig data
    conBasesegmentList.push_back(bs);
  }

  // restore current position
  fseek64(gig, currentOffset, SEEK_SET);

  //----------------------------------------------------------------------------
  // load into ContigData
  //----------------------------------------------------------------------------
  contigData.name = conName;
  contigData.length = conLength;
  contigData.numReads = conReads;
  contigData.numBasesegments = conBasesegments;
  contigData.dna = conDna;
  contigData.qual = conQual;
  contigData.basesegments = conBasesegmentList;
  contigData.valid = true;

  //----------------------------------------------------------------------------
  // return contig data
  //----------------------------------------------------------------------------
  return contigData;
}

//------------------------------------------------------------------------------
// prepareReads
//------------------------------------------------------------------------------
bool GigReader::prepareReads(string conName) {

  //----------------------------------------------------------------------------
  // look up offset for this contig
  //----------------------------------------------------------------------------
  if (! (conOffsetHash.count(conName) > 0)) {
    return false;
  }
  off64_t conOffset = conOffsetHash[conName];

  //----------------------------------------------------------------------------
  // read contig info from contigOffset location
  //----------------------------------------------------------------------------

  // position to conOffset
  fseek64(gig, conOffset, SEEK_SET);
  
  // read number of assembled reads in contig
  int conReads; 
  fread((char*)&conReads, sizeof(int), 1, gig);

  // read offset for read index
  off64_t readIndexOffset;
  fread((char*)&readIndexOffset, sizeof(off64_t), 1, gig);

  // read offset for read block
  off64_t readBlockOffset;
  fread((char*)&readBlockOffset, sizeof(off64_t), 1, gig);

  //----------------------------------------------------------------------------
  // position to readBlockOffset
  //----------------------------------------------------------------------------
  fseek64(gig, readBlockOffset, SEEK_SET);

  //----------------------------------------------------------------------------
  // return true
  //----------------------------------------------------------------------------
  return true;
}

//------------------------------------------------------------------------------
// getReads
//------------------------------------------------------------------------------
unsigned int GigReader::getReads(int numReadsToRead, deque<ReadData> & reads) {

  //----------------------------------------------------------------------------
  // it is assumed that the current offset was set properly
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  // loop through all reads
  //----------------------------------------------------------------------------
  for (int readCount=1; readCount<=numReadsToRead; readCount++) {

    // read length of readName
    int readNameCstyleLength; 
    fread((char*)&readNameCstyleLength, sizeof(int), 1, gig);
    
    // read read length
    int readLength; 
    fread((char*)&readLength, sizeof(int), 1, gig);
    
    // read contig fit left coordinate
    int contigFitLeft;
    fread((char*)&contigFitLeft, sizeof(int), 1, gig);
    
    // read contig fit right coordinate
    int contigFitRight;
    fread((char*)&contigFitRight, sizeof(int), 1, gig);
    
    // read quality clip left coordinate
    int readQualClipLeft;
    fread((char*)&readQualClipLeft, sizeof(int), 1, gig);
    
    // read quality clip right coordinate
    int readQualClipRight;
    fread((char*)&readQualClipRight, sizeof(int), 1, gig);
    
    // read alignment clip left coordinate
    int readAlignmentClipLeft;
    fread((char*)&readAlignmentClipLeft, sizeof(int), 1, gig);
    
    // read alignment clip right coordinate
    int readAlignmentClipRight;
    fread((char*)&readAlignmentClipRight, sizeof(int), 1, gig);
    
    // read if alignment is complemented
    bool complemented;
    fread((char*)&complemented, sizeof(bool), 1, gig);
    
    // read readName
    char readNameCstyle[readNameCstyleLength];
    fread(readNameCstyle, sizeof(char), readNameCstyleLength+1, gig);
    
    // read readDna (plus C style terminator to be read but ignored)
    string readDna;
    char b;
    for (int i=0; i<readLength; i++) {
      fread(&b, sizeof(char), 1, gig);
      readDna += b;
    }
    fread(&b, sizeof(char), 1, gig);
    
    // read base quality values
    vector<short> readQual;
    for (int i=0; i<readLength; i++) {
      short q;
      fread((char*)&q, sizeof(short), 1, gig);
      readQual.push_back(q);
    }
    
    //--------------------------------------------------------------------------
    // load into readData
    //--------------------------------------------------------------------------
    ReadData readData;
    readData.name = readNameCstyle;
    readData.dna = readDna;
    readData.qual = readQual;
    readData.length = readLength;
    readData.contigFitLeft = contigFitLeft;
    readData.contigFitRight = contigFitRight;
    readData.readQualClipLeft = readQualClipLeft;
    readData.readQualClipRight = readQualClipRight;
    readData.readAlignmentClipLeft = readAlignmentClipLeft;
    readData.readAlignmentClipRight = readAlignmentClipRight;
    readData.complemented = complemented;

    //--------------------------------------------------------------------------
    // register readData in list
    //--------------------------------------------------------------------------
    reads.push_back(readData);
  }
  
  //----------------------------------------------------------------------------
  // return
  //----------------------------------------------------------------------------
  return numReadsToRead;
}

//------------------------------------------------------------------------------
// getNextRead
//------------------------------------------------------------------------------
bool GigReader::getNextRead(ReadData & readData) {

  //----------------------------------------------------------------------------
  // read read info 
  // it is assumed that the current offset was set properly
  //----------------------------------------------------------------------------

  // read length of readName
  int readNameCstyleLength; 
  fread((char*)&readNameCstyleLength, sizeof(int), 1, gig);
  
  // read read length
  int readLength; 
  fread((char*)&readLength, sizeof(int), 1, gig);
  
  // read contig fit left coordinate
  int contigFitLeft;
  fread((char*)&contigFitLeft, sizeof(int), 1, gig);
  
  // read contig fit right coordinate
  int contigFitRight;
  fread((char*)&contigFitRight, sizeof(int), 1, gig);
  
  // read quality clip left coordinate
  int readQualClipLeft;
  fread((char*)&readQualClipLeft, sizeof(int), 1, gig);
  
  // read quality clip right coordinate
  int readQualClipRight;
  fread((char*)&readQualClipRight, sizeof(int), 1, gig);
  
  // read alignment clip left coordinate
  int readAlignmentClipLeft;
  fread((char*)&readAlignmentClipLeft, sizeof(int), 1, gig);
  
  // read alignment clip right coordinate
  int readAlignmentClipRight;
  fread((char*)&readAlignmentClipRight, sizeof(int), 1, gig);
  
  // read if alignment is complemented
  bool complemented;
  fread((char*)&complemented, sizeof(bool), 1, gig);
  
  // read readName
  char readNameCstyle[readNameCstyleLength];
  fread(readNameCstyle, sizeof(char), readNameCstyleLength+1, gig);
  
  // read readDna (plus C style terminator to be read but ignored)
  string readDna;
  char b;
  for (int i=0; i<readLength; i++) {
    fread(&b, sizeof(char), 1, gig);
    readDna += b;
  }
  fread(&b, sizeof(char), 1, gig);

  // read base quality values
  vector<short> readQual;
  for (int i=0; i<readLength; i++) {
    short q;
    fread((char*)&q, sizeof(short), 1, gig);
    readQual.push_back(q);
  }

  //----------------------------------------------------------------------------
  // load into readData
  //----------------------------------------------------------------------------
  readData.name = readNameCstyle;
  readData.dna = readDna;
  readData.qual = readQual;
  readData.length = readLength;
  readData.contigFitLeft = contigFitLeft;
  readData.contigFitRight = contigFitRight;
  readData.readQualClipLeft = readQualClipLeft;
  readData.readQualClipRight = readQualClipRight;
  readData.readAlignmentClipLeft = readAlignmentClipLeft;
  readData.readAlignmentClipRight = readAlignmentClipRight;
  readData.complemented = complemented;
  
  //----------------------------------------------------------------------------
  // return
  //----------------------------------------------------------------------------
  return true;
}

//------------------------------------------------------------------------------
// getReadData
//------------------------------------------------------------------------------
ReadData GigReader::getReadData(string contigName, string readName) {

  //----------------------------------------------------------------------------
  // make ReadData output
  //----------------------------------------------------------------------------
  ReadData readData;
  readData.valid = false;

  //----------------------------------------------------------------------------
  // look up offset for this read
  //----------------------------------------------------------------------------
  if ((! conAssReadOffset.count(contigName) > 0) || (! conAssReadOffset[contigName].count(readName) > 0)) {
    return readData;
  }
  off64_t readOffset = conAssReadOffset[contigName][readName];

  //----------------------------------------------------------------------------
  // store current offset
  //----------------------------------------------------------------------------
  off64_t currentOffset = ftell64(gig);

  //----------------------------------------------------------------------------
  // read read info from readOffset location
  //----------------------------------------------------------------------------

  // position to offset in file
  fseek64(gig, readOffset, SEEK_SET);
  
  // read length of readName
  int readNameCstyleLength; 
  fread((char*)&readNameCstyleLength, sizeof(int), 1, gig);
  
  // read read length
  int readLength; 
  fread((char*)&readLength, sizeof(int), 1, gig);
  
  // read contig fit left coordinate
  int contigFitLeft;
  fread((char*)&contigFitLeft, sizeof(int), 1, gig);
  
  // read contig fit right coordinate
  int contigFitRight;
  fread((char*)&contigFitRight, sizeof(int), 1, gig);
  
  // read quality clip left coordinate
  int readQualClipLeft;
  fread((char*)&readQualClipLeft, sizeof(int), 1, gig);
  
  // read quality clip right coordinate
  int readQualClipRight;
  fread((char*)&readQualClipRight, sizeof(int), 1, gig);
  
  // read alignment clip left coordinate
  int readAlignmentClipLeft;
  fread((char*)&readAlignmentClipLeft, sizeof(int), 1, gig);
  
  // read alignment clip right coordinate
  int readAlignmentClipRight;
  fread((char*)&readAlignmentClipRight, sizeof(int), 1, gig);
  
  // read if alignment is complemented
  bool complemented;
  fread((char*)&complemented, sizeof(bool), 1, gig);
  
  // read readName
  char readNameCstyle[readNameCstyleLength];
  fread(readNameCstyle, sizeof(char), readNameCstyleLength+1, gig);
  
  // read readDna (plus C style terminator to be read but ignored)
  string readDna;
  char b;
  for (int i=0; i<readLength; i++) {
    fread(&b, sizeof(char), 1, gig);
    readDna += b;
  }
  fread(&b, sizeof(char), 1, gig);

  // read base quality values
  vector<short> readQual;
  for (int i=0; i<readLength; i++) {
    short q;
    fread((char*)&q, sizeof(short), 1, gig);
    readQual.push_back(q);
  }
  
  //----------------------------------------------------------------------------
  // load into readData
  //----------------------------------------------------------------------------
  readData.name = readNameCstyle;
  readData.dna = readDna;
  readData.qual = readQual;
  readData.length = readLength;
  readData.contigFitLeft = contigFitLeft;
  readData.contigFitRight = contigFitRight;
  readData.readQualClipLeft = readQualClipLeft;
  readData.readQualClipRight = readQualClipRight;
  readData.readAlignmentClipLeft = readAlignmentClipLeft;
  readData.readAlignmentClipRight = readAlignmentClipRight;
  readData.complemented = complemented;
  readData.valid = true;
  
  //----------------------------------------------------------------------------
  // restore current position
  //----------------------------------------------------------------------------
  fseek64(gig, currentOffset, SEEK_SET);
  
  // return
  return readData;
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// GigReader class methods
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

#endif

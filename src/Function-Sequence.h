//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// Function-Sequence
// DNA sequence manipulation methods
// Copyright 2006,2007 Gabor T. Marth, Boston College
// All rights reserved
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

#ifndef FUNCTION_SEQUENCE_H
#define FUNCTION_SEQUENCE_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <iterator>
#include <algorithm>

using std::ios;
using std::ifstream;
using std::ofstream;
using std::ostream;
using std::istream;
using std::cin;
using std::cout;
using std::clog;
using std::endl;
using std::string;
using std::vector;
using std::map;

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// type definitions
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// Sample: information describing a sample
//------------------------------------------------------------------------------
struct SampleInfo {
  string sampleId;
  string pedigreeId;
  string motherId;
  string fatherId;
  string sex;
};

//------------------------------------------------------------------------------
// unpadDna -- unpads dna sequence 
//------------------------------------------------------------------------------
string unpadDna(const string);

//------------------------------------------------------------------------------
// makeBaseMap -- makes base map between padded and unpadded dna 
//------------------------------------------------------------------------------
vector< vector<int> > makeBaseMap(const string);

//------------------------------------------------------------------------------
// padDna -- pads dna sequence according to base map
//------------------------------------------------------------------------------
string padDna(const string, const vector<int>);

//------------------------------------------------------------------------------
// padBaseQual -- pads base quality sequence according to base map
//------------------------------------------------------------------------------
vector<short> padQual(const vector<short>, const vector<int>);

//------------------------------------------------------------------------------
// revCompDna -- reverse complements dna sequence 
//------------------------------------------------------------------------------
string revCompDna(const string);

//------------------------------------------------------------------------------
// lowerCaseDna -- converts dna sequence into lower case 
//------------------------------------------------------------------------------
string lowerCaseDna(const string);

//------------------------------------------------------------------------------
// upperCaseDna -- converts dna sequence into upper case 
//------------------------------------------------------------------------------
string upperCaseDna(const string);

//------------------------------------------------------------------------------
// dnaGenotypes -- returns possible genotypes of given multiplicity
//------------------------------------------------------------------------------
vector<string> dnaGenotypes(const int, const bool);

//------------------------------------------------------------------------------
// dnaGenotypeIsValid -- returns bool according to whether or not input genotype
//                       is valid
//------------------------------------------------------------------------------
bool dnaGenotypeIsValid(const int, const bool, const string);

//------------------------------------------------------------------------------
// getDnaGenotypeIndex -- returns the int index of the input genotype
//                     -- return the int index of the complement if complement is true
//------------------------------------------------------------------------------
int getDnaGenotypeIndex(const int, string, bool);

//------------------------------------------------------------------------------
// extractSampleInfo -- extracts sample information from sequence name
//------------------------------------------------------------------------------
SampleInfo extractSampleInfo (const string, const string, const string);

//------------------------------------------------------------------------------
// extractTemplateName -- extracts template name from sequence name
//------------------------------------------------------------------------------
string extractTemplateName (const string, const string);

//------------------------------------------------------------------------------
// gapMapFromDna -- makes a map of gaps in padded dna
//------------------------------------------------------------------------------
map<int, int, std::less<int> > gapMapFromDna(const string);

//------------------------------------------------------------------------------
// unpaddedPosMap -- makes a map from padded to unpadded dna
//------------------------------------------------------------------------------
vector<int> unpaddedPosMap(const string);

//------------------------------------------------------------------------------
// unpaddedPosMapBegin -- makes a map from padded to unpadded dna (left)
//------------------------------------------------------------------------------
vector<int> unpaddedPosMapBegin(const string);

//------------------------------------------------------------------------------
// unpaddedPosMapEnd -- makes a map from padded to unpadded dna (right)
//------------------------------------------------------------------------------
vector<int> unpaddedPosMapEnd(const string);

//------------------------------------------------------------------------------
// paddedPosMap -- makes a map from unpadded to padded dna
//------------------------------------------------------------------------------
vector<int> paddedPosMap(const int, map<int, int, std::less<int> >);
  
//------------------------------------------------------------------------------
// padDnaGapMap -- pad dna sequence
//------------------------------------------------------------------------------
string padDnaGapMap(const string, map<int, int, std::less<int> >);

//------------------------------------------------------------------------------
// padQualGapMap -- pad qual sequence
//------------------------------------------------------------------------------
vector<short> padQualGapMap(const vector<short>, map<int, int, std::less<int> >);

//------------------------------------------------------------------------------
// printDna -- prints dna sequence to ostream
//------------------------------------------------------------------------------
void printDna(ostream &, const string, const int);

//------------------------------------------------------------------------------
// printDnaFasta -- prints dna sequence to ostream in FASTA format
//------------------------------------------------------------------------------
void printDnaFasta(ostream &, const string, const string, const int);

//------------------------------------------------------------------------------
// printQual -- prints qual sequence to ostream
//------------------------------------------------------------------------------
void printQual(ostream &, const vector<short>, const int);

//------------------------------------------------------------------------------
// printQual -- prints qual sequence to ostream
//------------------------------------------------------------------------------
void printQualFasta(ostream &, const string, const vector<short>, const int);

//------------------------------------------------------------------------------
// printBpos -- prints bpos sequence to ostream
//------------------------------------------------------------------------------
void printBpos(ostream &, const vector<int>, const int);

//------------------------------------------------------------------------------
// printQual -- prints bpos sequence to ostream
//------------------------------------------------------------------------------
void printBposFasta(ostream &, const string, const vector<int>, const int);

#endif

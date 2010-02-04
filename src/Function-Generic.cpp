//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// Function-Generic
// Generic methods
// Copyright 2006, 2007 Gabor T. Marth, Boston College
// All rights reserved
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

#ifndef FUNCTION_GENERIC_CPP
#define FUNCTION_GENERIC_CPP

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
using std::multimap;

#include "Type-Hash.h"
#include "Function-Generic.h"

//------------------------------------------------------------------------------
// prints a vector
//------------------------------------------------------------------------------
template < class T >
void printVector(const vector< T > &vectorRef, const bool separator, const bool bracket, ostream &output2) {
  if ( vectorRef.empty() ) {
    if (bracket) {
      output2 << "()";
    }
  }
  else {
    if (bracket) {
      output2 << "(";
    }
    if (separator) {
      std::ostream_iterator< T > output( output2, " " );
      std::copy( vectorRef.begin(), vectorRef.end(), output );
    }
    else {
      std::ostream_iterator< T > output( output2 );
      std::copy( vectorRef.begin(), vectorRef.end(), output );
    }
    if (bracket) {
      output2 << ")";
    }
  }
}

//------------------------------------------------------------------------------
// converts a string class string to a bool
//------------------------------------------------------------------------------
bool string2Bool (string sString) {

  if (sString == "1" || sString == "true") {
    return true;
  }
  else {
    return false;
  }
}

//------------------------------------------------------------------------------
// converts a string class string to a short
//------------------------------------------------------------------------------
short string2Short (string sString) {

  // convert it to C-style string
  const char *cString = sString.c_str();

  // convert it to integer
  short iString = atoi(cString);

  // return
  return iString;
}

//------------------------------------------------------------------------------
// converts a string class string to an int
//------------------------------------------------------------------------------
int string2Int (string sString) {

  // convert it to C-style string
  const char *cString = sString.c_str();

  // convert it to integer
  int iString = atoi(cString);

  // return
  return iString;
}

//------------------------------------------------------------------------------
// converts a string class string to a long int
//------------------------------------------------------------------------------
long string2Long (string sString) {

  // convert it to C-style string
  const char *cString = sString.c_str();

  // convert it to integer
  long lString = atol(cString);

  // return
  return lString;
}

//------------------------------------------------------------------------------
// converts a string class string to a long long
//------------------------------------------------------------------------------
long long string2LongLong (string sString) {

  // convert it to C-style string
  const char *cString = sString.c_str();

  // convert it to integer
  long long llString = atol(cString);

  // return
  return llString;
}

//------------------------------------------------------------------------------
// converts a string class string to a double
//------------------------------------------------------------------------------
double string2Double (string sString) {

  // convert it to C-style string
  const char *cString = sString.c_str();

  // convert it to integer
  double fString = atof(cString);

  // return
  return fString;
}

//------------------------------------------------------------------------------
// converts a string class string to an double
//------------------------------------------------------------------------------
long double string2LongDouble (string sString) {

  // convert it to C-style string
  const char *cString = sString.c_str();

  // convert it to integer
  long double fString = atof(cString);

  // return
  return fString;
}

//------------------------------------------------------------------------------
// sorts keys of a hash in order of associated value
//------------------------------------------------------------------------------
template< typename keyType, typename valueType >
vector<keyType> sortKeysByValue(map<keyType, valueType, std::less<keyType> > hash, bool descend) {

  // instantiate inverse hash as a multimap
  multimap<valueType, keyType, std::less<valueType> > inverseHash;

  // load elements of hash into inverseHash
  for (typename map<keyType, valueType, std::less<keyType> >::const_iterator iter = hash.begin(); iter != hash.end(); iter++) {
    keyType key = iter->first;
    valueType value = iter->second;
    inverseHash.insert(typename multimap<valueType, keyType, std::less<valueType> >::value_type(value, key));
  }

  // compose vector of original keys sorted by original values
  vector<keyType> sortedKeys;
  for(typename multimap<valueType, keyType, std::less<valueType> >::const_iterator iter = inverseHash.begin(); 
      iter != inverseHash.end(); iter++) {
    keyType key = iter->second;
    sortedKeys.push_back(key);
  }

  // reverse if descending order was required
  if (descend) {
    reverse(sortedKeys.begin(), sortedKeys.end());
  }

  // return
  return sortedKeys;
}

#endif

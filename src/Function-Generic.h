//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// Function-Generic
// Generic methods
// Copyright 2006 Gabor T. Marth, Boston College
// All rights reserved
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

#ifndef FUNCTION_GENERIC_H
#define FUNCTION_GENERIC_H

#include <iostream>
#include <ostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <iterator>
#include <algorithm>

// private
#include "Type-Hash.h"

using namespace std;

//------------------------------------------------------------------------------
// prints a vector
//------------------------------------------------------------------------------
template < class T >
void printVector(const vector< T > &, const bool, const bool, ostream &);

//------------------------------------------------------------------------------
// converts a string class string to a bool
//------------------------------------------------------------------------------
bool string2Bool (string);

//------------------------------------------------------------------------------
// converts a string class string to a short
//------------------------------------------------------------------------------
short string2Short (string);

//------------------------------------------------------------------------------
// converts a string class string to an int
//------------------------------------------------------------------------------
int string2Int (string);

//------------------------------------------------------------------------------
// converts a string class string to a long
//------------------------------------------------------------------------------
long string2Long (string);

//------------------------------------------------------------------------------
// converts a string class string to a long long
//------------------------------------------------------------------------------
long long string2LongLong (string);

//------------------------------------------------------------------------------
// converts a string class string to a double
//------------------------------------------------------------------------------
double string2Double (string);

//------------------------------------------------------------------------------
// converts a string class string to a long double
//------------------------------------------------------------------------------
long double string2LongDouble (string);

//------------------------------------------------------------------------------
// returns keys of a hash in order of associated value
//------------------------------------------------------------------------------
template< typename keyType, typename valueType >
vector<keyType> sortKeysByValue(map<keyType, valueType, less<keyType> >, bool);

#endif

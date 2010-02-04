//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// Type-Hash
// Type definitions
// Copyright 2004-2007 Gabor T. Marth, Boston College
// All rights reserved
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

#include <string>
#include <vector>
#include <map>

using std::string;
using std::vector;
using std::map;
using namespace __gnu_cxx;

// "hash_map" true hashes
#include <ext/hash_map>
/*
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
*/

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// single-key hash (map) definitions
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

// key: int; value: int
typedef std::map<int, int, std::less<int> > map_int_int;

// key: int; value: double
typedef std::map<int, double, std::less<int> > map_int_double;

// key: int; value: long double
typedef std::map<int, long double, std::less<int> > map_int_longDouble;

// key: int; value: bool
typedef std::map<int, bool, std::less<int> > map_int_bool;

// key: int; value: string
typedef std::map<int, string, std::less<int> > map_int_string;


// key: string; value: int
typedef std::map<string, int, std::less<string> > map_string_int;

// key: string; value: string
typedef std::map<string, string, std::less<string> > map_string_string;

// key: string; value: double
typedef std::map<string, double, std::less<string> > map_string_double;

// key: string; value: long double
typedef std::map<string, long double, std::less<string> > map_string_longDouble;

// key: string; value: long long
typedef std::map<string, long long, std::less<string> > map_string_longLong;

// key: string; value: bool
typedef std::map<string, bool, std::less<string> > map_string_bool;

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// multi-vallue hash (multimap) definitions
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
typedef std::multimap<int, string, std::less<int> > multimap_int_string;

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// multi-key hash (map) definitions
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

// key: vector<int>; value: int
typedef std::map< vector<int>, int, std::less< vector<int> > > mdim_map_int;

// key: vector<int>; value: double
typedef std::map< vector<int>, double, std::less< vector<int> > > mdim_map_double;

// key: vector<int>; value: long double
typedef std::map< vector<int>, long double, std::less< vector<int> > > mdim_map_longDouble;

// key: vector<int>; value: bool
typedef std::map< vector<int>, bool, std::less< vector<int> > > mdim_map_bool;

// key: vector<int>; value: string
typedef std::map< vector<int>, bool, std::less< vector<int> > > mdim_map_string;

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// complex value hash (map) definitions
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

// key: string; value: vector<short>
typedef std::map<string, vector<short>, std::less<string> > map_string_vectorShort;


#ifndef SEQLIB_UTILS_H
#define SEQLIB_UTILS_H

#include <string>
#include <time.h>
#include <ctime>
#include <vector>
#include <unistd.h>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <stdio.h>

#include "SeqLib/SeqLibCommon.h"

    
#if __cplusplus > 199711L
  #include <memory>
  #include <unordered_set>
  #include <unordered_map>
  #define SeqHashMap std::unordered_map
  #define SeqHashSet std::unordered_set
  #define SeqPointer std::shared_ptr
  #define HAVE_C11 1
#else

#ifdef __APPLE__
  #include <memory>
  #include <unordered_set>
  #include <unordered_map>
  #define SeqHashMap std::unordered_map
  #define SeqHashSet std::unordered_set
  #define SeqPointer std::shared_ptr
#else
  #include <tr1/memory>
  #include <tr1/unordered_set>
  #include <tr1/unordered_map>
  #define SeqHashMap std::tr1::unordered_map
  #define SeqHashSet std::tr1::unordered_set
  #define SeqPointer std::tr1::shared_ptr
#endif
#endif

namespace SeqLib {

  template<typename T> 
    inline std::string tostring(T d) { 
      std::stringstream ss;
      ss << d;
    return ss.str();
  }
  
  /** Check if a file is readable and exists
   * @param name Name of a file to test
   * @return File is readable and exists
   */
  inline bool read_access_test (const std::string& name) {
    return (access (name.c_str(), R_OK) == 0); 
  }

  /** Format an integer to include commas
   * @param data Number to format
   * @return String with formatted number containing commas
   */
  template <typename T> inline
    std::string AddCommas(T data) {
    std::stringstream ss; 
    ss << data; 
    std::string s = ss.str();
    if (s.length() > 3)
      for (int i = s.length()-3; i > 0; i -= 3)
	s.insert(i,",");
    return s;
  }

  /** Display the runtime (CPU and Wall)
   * 
   * @param start Running timer
   * @return Time formatted as "CPU: XmYs Wall: XmYs"
   * @note Does not work on OSX or Windows (returns "not configured")
   */
  inline std::string displayRuntime(
#ifndef __APPLE__
				    const timespec start
#endif
				    ) {
    
#ifndef __APPLE__
    struct timespec finish;
    clock_gettime(CLOCK_MONOTONIC, &finish);
    double elapsed = (finish.tv_sec - start.tv_sec);
    int t = clock()/CLOCKS_PER_SEC;
    int min = (int)std::floor(elapsed / 60.0);
    int sec = (int)(elapsed-min*60);
    char buffer[100];
    sprintf (buffer, "CPU: %4dm%02ds Wall: %4dm%02ds", 
	     (int)floor( ((double)t) /60.0), t % 60, min, sec);
    buffer[99] = '\0';
    return std::string(buffer);
#else
    return "--- time not configured for apple\n";
#endif
  }

  /** Reverse complement in-place sequence containg upper/lower case ACTGN
   * @param a Sequence to be reverse complmented
   */
  inline void rcomplement(std::string &a) {
    
    std::reverse(&a[0], &a[a.size()]);
    std::string::iterator it = a.begin();
    for (; it != a.end(); it++)
      *it = RCOMPLEMENT_TABLE[(unsigned char)*it];
  }
  

  /** Calculate the percentage and return as integer
   * @param numer Numerator
   * @param denom Denominator
   * @return Integer with the percentage floored
   */
  template <typename T> inline int percentCalc(T numer, T denom) {
    if (denom <= 0)
      return 0;
    int perc = numer * 100 / denom;
    //int perc  = static_cast<int>(floor((float)numer / (float)denom * 100.0));
    return perc;
  }

  /** Remove substrings from a string
   * @param toscrub Input string to clean
   * @param toremove Substring to remove from input
   * @return Scrubbed string
   */
  inline std::string scrubString(const std::string& toscrub, const std::string& toremove) 
  {
    if (toscrub.empty() || toremove.empty())
      return toscrub;
    
    std::string::size_type i = toscrub.find(toremove);
    if (i == std::string::npos)
      return toscrub;
    
    std::string ts = toscrub;
    while (i != std::string::npos) {
      ts.erase(i, toremove.length());
      i = ts.find(toremove);
    }
    return ts;
  }
  
  // Generate a weighed random integer 
  // @param cs Weighting for each integer (values must sum to one) 
  // @return Random integer bounded on [0,cs.size())
  //
  //int weightedRandom(const std::vector<double>& cs);

}

#endif

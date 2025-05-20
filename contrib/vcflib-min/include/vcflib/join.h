/*
    vcflib C++ library for parsing and manipulating VCF files

    Copyright © 2010-2020 Erik Garrison
    Copyright © 2020      Pjotr Prins

    This software is published under the MIT License. See the LICENSE file.
*/

#ifndef __JOIN_H
#define __JOIN_H

// functions to split a string by a specific delimiter
#include <string>
#include <vector>
#include <sstream>
#include <string.h>
#include <list>

// join a vector of elements by a delimiter object.  ostream<< must be defined
// for both class S and T and an ostream, as it is e.g. in the case of strings
// and character arrays
template<class S, class T>
std::string join(const std::vector<T>& elems, S& delim) {
    std::stringstream ss;
    typename std::vector<T>::const_iterator e = elems.begin();
    if (e != elems.end()) {
        ss << *e++;
        for (; e != elems.end(); ++e) {
            ss << delim << *e;
        }
    }
    return ss.str();
}

// same for lists
template<class S, class T>
std::string join(std::list<T>& elems, S& delim) {
    std::stringstream ss;
    typename std::list<T>::iterator e = elems.begin();
    if (e != elems.end()) {
        ss << *e++;
        for (; e != elems.end(); ++e) {
            ss << delim << *e;
        }
    }
    return ss.str();
}

#endif

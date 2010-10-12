// functions to split a string by a specific delimiter
#include <string>
#include <vector>
#include <sstream>
#include <string.h>

// join a vector of elements by a delimiter object.  ostream<< must be defined
// for both class S and T and an ostream, as it is e.g. in the case of strings
// and character arrays
template<class S, class T>
std::string join(std::vector<T>& elems, S& delim) {
    std::stringstream ss;
    typename std::vector<T>::iterator e = elems.begin();
    ss << *e++;
    for (; e != elems.end(); ++e) {
        ss << delim << *e;
    }
    return ss.str();
}

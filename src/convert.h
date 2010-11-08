#include <sstream>

// converts the string into the specified type, setting r to the converted
// value and returning true/false on success or failure
template<typename T>
bool convert(const string& s, T& r) {
    istringstream iss(s);
    iss >> r;
    return ((iss.fail() || iss.tellg() != s.size())) ? false : true;
}

#include <sstream>

// converts the string into the specified type, setting r to the converted
// value and returning true/false on success or failure
template<typename T>
bool convert(const std::string& s, T& r) {
    std::istringstream iss(s);
    iss >> r;
    return (iss.fail() || ((std::size_t) iss.tellg()) != s.size()) ? false : true;
}

#ifndef __SUM_H
#define __SUM_H

#include <vector>

template <class T>
T sum(const std::vector<T>& v) {
    T result = 0;
    for (typename std::vector<T>::const_iterator i = v.begin(); i != v.end(); ++i) {
        result += *i;
    }
    return result;
}

#endif

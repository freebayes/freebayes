#ifndef __PRODUCT_H
#define __PRODUCT_H

#include <vector>

template <class T>
T product(const std::vector<T>& v) {
    if (v.size() > 0) {
        T result = 1;
        for (typename std::vector<T>::const_iterator i = v.begin(); i != v.end(); ++i) {
            result *= *i;
        }
        return result;
    } else {
        return 0;
    }
}

#endif

/* 

multipermute.h  -- multiset permutations for generic vectors

Follows 'Algorithm 1' from "Loopless Generation of Multiset Permutations using
a Constant Number of Variables by Prefix Shifts."  Aaron Williams, 2009

author: Erik Garrison <erik.garrison@bc.edu>
last revised: 2010-04-16

Copyright (c) 2010 by Erik Garrison

Permission is hereby granted, free of charge, to any person
obtaining a copy of this software and associated documentation
files (the "Software"), to deal in the Software without
restriction, including without limitation the rights to use,
copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following
conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
OTHER DEALINGS IN THE SOFTWARE.

*/


#include <vector>
#include <algorithm>

template <class T>
class ListElement { 

public:
    T value;
    ListElement<T>* next;

    ListElement<T>() { }

    ListElement<T>(T val, ListElement<T>* n) {
        value = val;
        next = n;
    }

    ListElement<T>* nth(int n) {
        ListElement<T>* o = this;
        int i = 0;
        while (i < n && o->next != NULL) {
            o = o->next;
            ++i;
        }
        return o;
    }

    ~ListElement<T>() {
        if (next != NULL) {
            delete next;
        }
    }

};

template <class T>
ListElement<T>* list_init(std::vector<T>& multiset) {
    std::sort(multiset.begin(), multiset.end()); // ensures proper non-increasing order
    typename std::vector<T>::const_iterator item = multiset.begin();
    ListElement<T>* h = new ListElement<T>(*item, NULL);
    ++item;
    while (item != multiset.end()) {
        h = new ListElement<T>(*item, h);
        ++item;
    }
    return h;
}

template <class T>
std::vector<T> linked_list_to_vector(ListElement<T>* h) {
    ListElement<T>* o = h;
    std::vector<T> l;
    while (o != NULL) {
        l.push_back(o->value);
        o = o->next;
    }
    return l;
}

// provides multiset permutations out of the std::vector multiset
template <class T>
std::vector< std::vector<T> > multipermute(std::vector<T>& multiset) {

    std::vector< std::vector<T> > results;

    ListElement<T>* h = list_init(multiset);
    ListElement<T>* i = h->nth(multiset.size() - 2);
    ListElement<T>* j = h->nth(multiset.size() - 1);
    ListElement<T>* s;
    ListElement<T>* t;

    results.push_back(linked_list_to_vector(h));

    while (j->next != NULL || j->value < h->value) {
        if (j->next != NULL && i->value >= j->next->value) {
            s = j;
        } else {
            s = i;
        }
        t = s->next;
        s->next = t->next;
        t->next = h;
        if (t->value < h->value) {
            i = t;
        }
        j = i->next;
        h = t;
        results.push_back(linked_list_to_vector(h));
    }

    delete h;

    return results;

}


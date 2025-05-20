#ifndef __MULTICHOOSE_H
#define __MULTICHOOSE_H

/*

multichoose.h  -- n multichoose k for generic vectors

author: Erik Garrison <erik.garrison@bc.edu>
last revised: 2010-04-16
author: Pjotr Prins
last revised: 2024 (fixed template error and merged template into vcflib)

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


// provides multiset combinations out of the std::vector of objects
template <class T>
std::vector< std::vector<T> > multichoose(int k, std::vector<T>& objects) {

    std::vector< std::vector<T> > choices;

    int j,j_1,q,r;

    r = objects.size() - 1;

    // combination indexes
    std::vector<T*> a, b;

    for (int i=0;i<k;i++) {
        a.push_back(&objects[0]); b.push_back(&objects[r]);
    }

    j=k;
    while(true){
        std::vector<T> multiset;
        multiset.reserve(k);
        for(const auto elm : a)
            multiset.push_back(*elm);
        choices.push_back(multiset);
        j=k;
        do {
	    j--;
	    if (j<0) break;
	} while(a[j]==b[j]);
        if (j<0) break;
        j_1=j;
        while(j_1<=k-1){
            a[j_1]=a[j_1]+1;
            q=j_1;
            while(q<k-1) {
                a[q+1]=a[q];
                q++;
            }
            q++;
            j_1=q;
        }
    }

    return choices;
}

#endif

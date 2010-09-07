// alleles.cpp
// outputs a json-formatted stream of alleles over target regions
//
// Erik Garrison <erik.garrison@bc.edu>
// Marth Lab, Boston College
// July 14, 2010
//
// standard includes
//#include <cstdio>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <iterator>
#include <algorithm>
#include <cmath>
#include <time.h>

// "boost" regular expression library
#include <boost/regex.hpp>

// "boost" string manipulation
#include <boost/algorithm/string/join.hpp>
// tuple
#include <boost/tuple/tuple.hpp>
// bind
#include <boost/bind.hpp>

// private libraries
#include "BamReader.h"
#include "Fasta.h"
#include "TryCatch.h"
#include "Parameters.h"
#include "Allele.h"
#include "AlleleParser.h"

#include "multichoose.h"
#include "multipermute.h"

using namespace std; 

using boost::tuple;
using boost::make_tuple;

// Allele object recycling:
//
// We use the Allele freelist for performance reasons.  When an Allele object
// is destroyed, it is pushed onto this freelist.  When a new Allele object is
// created, new first checks if we have a free Allele object on the freelist.
// Because we are dynamically linked, we have to declare the freelist here,
// although it exists as a static member of the Allele class.
//
AlleleFreeList Allele::_freeList;

int main (int argc, char *argv[]) {

    AlleleParser* parser = new AlleleParser(argc, argv);
    list<Allele*> alleles;

    int allowedAlleleTypes = ALLELE_REFERENCE | ALLELE_SNP;

    while (parser->getNextAlleles(alleles)) {

        filterAlleles(alleles, allowedAlleleTypes);

        // skips 0-coverage regions
        if (alleles.size() == 0)
            continue;

        map<string, vector<Allele*> > sampleGroups = groupAllelesBySample(alleles);

        // report in json-formatted stream
        //
        cout << "{\"sequence\":\"" << parser->currentTarget->seq << "\","
            << "\"total coverage\":" << alleles.size() << ","
            << "\"position\":" << parser->currentPosition + 1 << ","  /// XXX basing somehow is 1-off... 
            << "\"reference base\":\"" << parser->currentReferenceBase << "\","
            //<< "\"raDepth\":" << parser->registeredAlleles.size() << ","
            << "\"samples\":{";  // TODO ... quality (~pSnp)

        bool suppressComma = true; // output flag
        for (map<string, vector<Allele*> >::iterator 
                sample = sampleGroups.begin();
                sample != sampleGroups.end(); ++sample) {

            if (!suppressComma) { cout << ","; } else { suppressComma = false; }

            cout << "\"" << sample->first << "\":{"
                << "\"coverage\":" << sample->second.size()
                << ",\"alleles\":" << json(sample->second)
                << "}";

        }

        cout << "}}" << endl;

    }

    delete parser;

    return 0;

}

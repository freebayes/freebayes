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

// private libraries
#include "BamReader.h"
#include "Fasta.h"
#include "TryCatch.h"
#include "Parameters.h"
#include "Allele.h"
#include "AlleleParser.h"
#include "split.h"

#include "multichoose.h"
#include "multipermute.h"

using namespace std; 

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

    int allowedAlleleTypes = ALLELE_REFERENCE | ALLELE_SNP | ALLELE_INSERTION | ALLELE_DELETION;

    map<string, vector<Allele*> > sampleGroups;

    while (parser->getNextAlleles(sampleGroups, allowedAlleleTypes)) {

        int coverage = countAlleles(sampleGroups);

        // skips 0-coverage regions
        if (coverage == 0)
            continue;

        // report in json-formatted stream
        //
        cout << "{\"sequence\":\"" << parser->currentTarget->seq << "\","
            << "\"total coverage\":" << coverage << ","
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

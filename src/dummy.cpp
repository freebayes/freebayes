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

#include "multichoose.h"
#include "multipermute.h"

using namespace std; 

AlleleFreeList Allele::_freeList;

int main (int argc, char *argv[]) {

    AlleleParser* parser = new AlleleParser(argc, argv);

    while (parser->dummyProcessNextTarget()) {
    }

    delete parser;

    return 0;

}

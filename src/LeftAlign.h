#ifndef __LEFTALIGN_H
#define __LEFTALIGN_H

#include <iostream>
#include <getopt.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <signal.h>
#include <stdlib.h>
#include <cmath>
#include <algorithm>
#include <map>
#include <vector>

#include "Fasta.h"
#include "api/BamAlignment.h"
#include "api/BamReader.h"
#include "api/BamWriter.h"

#include "IndelAllele.h"

#ifdef VERBOSE_DEBUG
#define LEFTALIGN_DEBUG(msg) \
    if (debug) { cerr << msg; }
#else
#define LEFTALIGN_DEBUG(msg)
#endif

using namespace std;
using namespace BamTools;

bool leftAlign(BamAlignment& alignment, string& referenceSequence, bool debug = false);
bool stablyLeftAlign(BamAlignment& alignment, string referenceSequence, int maxiterations = 20, bool debug = false);
int countMismatches(BamAlignment& alignment, string referenceSequence);

#endif

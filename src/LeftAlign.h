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

#ifdef HAVE_BAMTOOLS
#include "api/BamAlignment.h"
#include "api/BamReader.h"
#include "api/BamWriter.h"
using namespace BamTools;
#else
#include "SeqLib/BamReader.h"
#include "SeqLib/BamWriter.h"
#endif

#include "IndelAllele.h"

#ifdef VERBOSE_DEBUG
#define LEFTALIGN_DEBUG(msg) \
    if (debug) { cerr << msg; }
#else
#define LEFTALIGN_DEBUG(msg)
#endif

using namespace std;

#ifdef HAVE_BAMTOOLS
bool leftAlign(BamAlignment& alignment, string& referenceSequence, bool debug = false);
bool stablyLeftAlign(BamAlignment& alignment, string referenceSequence, int maxiterations = 20, bool debug = false);
int countMismatches(BamAlignment& alignment, string referenceSequence);
#else
bool leftAlign(SeqLib::BamRecord& alignment, string& referenceSequence, bool debug = false);
bool stablyLeftAlign(SeqLib::BamRecord& alignment, string referenceSequence, int maxiterations = 20, bool debug = false);
int countMismatches(SeqLib::BamRecord& alignment, string referenceSequence);
#endif


#endif

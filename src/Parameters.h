#ifndef _PARAMETERS_H
#define _PARAMETERS_H

#include <cstdio>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <getopt.h>
#include <stdlib.h>
#include "split.h"
#include "Version.h"
#include "Utility.h"

using namespace std;

// Encapsulates tracking and parsing of command line program parameters
class Parameters {

    friend ostream &operator<<(ostream &out, const Parameters &p);

public:

    // i/o parameters:
    string bam;                  // -b --bam
    vector<string> bams;
    bool useStdin;               // -c --stdin
    string fasta;                // -f --fasta-reference
    string targets;              // -t --targets
    string region;               // -r --region
    string samples;              // -s --samples
    //string log;
    string output;               // -v --vcf
    string outputFile;
    string traceFile;

    // operation parameters
    bool outputAlleles;          // -O --output-alleles
    bool trace;                  // -L --trace
    bool useDuplicateReads;      // -E --use-duplicate-reads
    bool suppressOutput;         // -S --suppress-output
    bool bamBayesDataLikelihoods;// -G --factorial-data-likelihoods
    int useBestNAlleles;         // -n --use-best-n-alleles
    bool forceRefAllele;         // -F --force-reference-allele
    bool useRefAllele;           // -U --use-reference-allele
    bool allowIndels;            // -I --allow-indels
    bool pooled;                 // -J --pooled
    int MQR;                     // -M --reference-mapping-quality
    int BQR;                     // -B --reference-base-quality
    int ploidy;                  // -p --ploidy
    int MQL0;                    // -m --min-mapping-quality
    int BQL0;                    // -q --min-base-quality
    int MQL1;                    // -R --min-supporting-mapping-quality
    int BQL1;                    // -S --min-supporting-base-quality
    int BQL2;                    // -Q --mismatch-base-quality-threshold
    int RMU;                     // -U --read-mismatch-limit
    int IDW;                     // -I --indel-exclusion-window
    long double TH;              // -T --theta
    long double PVL;             // -P --pvar
                                 // -K --posterior-integration-depth
    int posteriorIntegrationDepth;
    string algorithm;
    long double RDF;             // -D --read-dependence-factor
    int WB;                      // -W --posterior-integration-bandwidth
    // XXX adjusting this to anything other than 1 may have bad consequences
    // for large numbers of samples
    int TB;                      // -Y --posterior-integration-depth
    bool includeMonoB;
    int TR;
    int I;
    long double minAltFraction;  // -F --min-alternate-fraction
    int minAltCount;             // -C --min-alternate-count
    int debuglevel;              // -d --debug increments
    bool debug; // set if debuglevel >=1
    bool debug2; // set if debuglevel >=2

    // functions
    Parameters(int argc, char** argv);
    void usage(char **argv);

};

#endif

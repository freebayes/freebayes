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
#include "version_git.h"
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
    vector<string> regions;               // -r --region
    string samples;              // -s --samples
    string populationsFile;
    string cnvFile;
    //string log;
    string output;               // -v --vcf
    string outputFile;
    bool gVCFout;    // -l --gvcf
    int gVCFchunk;
    string variantPriorsFile;
    string haplotypeVariantFile;
    bool reportAllHaplotypeAlleles;
    bool reportMonomorphic;
    bool boundIndels;
    bool onlyUseInputAlleles;
    string alleleObservationBiasFile;
    bool standardGLs;
    int baseQualityCap;
    double probContamination;
    string contaminationEstimateFile;

    // operation parameters
    bool useDuplicateReads;      // -E --use-duplicate-reads
    bool suppressOutput;         // -S --suppress-output
    int useBestNAlleles;         // -n --use-best-n-alleles
    bool forceRefAllele;         // -F --force-reference-allele
    bool useRefAllele;           // -U --use-reference-allele
    bool diploidReference;       // -H --haploid-reference
    bool allowIndels;            // -I --allow-indels
    bool leftAlignIndels;        // -O --left-align-indels
    bool allowMNPs;              // -X --allow-mnps
    bool allowComplex;           // -X --allow-complex
    bool strictVCF;
    int maxComplexGap;
    //int maxHaplotypeLength;
    int minRepeatSize;
    double minRepeatEntropy;
    bool usePartialObservations;
    bool allowSNPs;              // -I --no-snps
    bool pooledDiscrete;
    bool pooledContinuous;
    bool ewensPriors;
    bool permute;                //    --permute
    bool useMappingQuality;      //
    bool useMinIndelQuality;
    bool obsBinomialPriors;
    bool alleleBalancePriors;
    bool hwePriors;
    bool reportGenotypeLikelihoodMax;
    int genotypingMaxIterations;
    int genotypingMaxBandDepth;
    bool excludePartiallyObservedGenotypes;
    bool excludeUnobservedGenotypes;
    float genotypeVariantThreshold;
    int siteSelectionMaxIterations;
    bool allSites; // TODO
    double limitGL; // minimum GL that is output
    int minPairedAltCount;
    double minAltMeanMapQ;
    int minAltQSum;
    int MQR;                     // -M --reference-mapping-quality
    int BQR;                     // -B --reference-base-quality
    int ploidy;                  // -p --ploidy
    int MQL0;                    // -m --min-mapping-quality
    int BQL0;                    // -q --min-base-quality
    int minSupportingMappingQualitySum;                    // -R --min-supporting-mapping-quality
    int minSupportingAlleleQualitySum;                    // -S --min-supporting-base-quality
    int BQL2;                    // -Q --mismatch-base-quality-threshold
    int RMU;                     // -U --read-mismatch-limit
    float readMaxMismatchFraction;  // -z --read-max-mismatch-fraction
    int readSnpLimit;            // -$ --read-snp-limit
    int readIndelLimit;          // -e --read-indel-limit
    int IDW;                     // -I --indel-exclusion-window
    long double TH;              // -T --theta
    long double PVL;             // -P --pvar
                                 // -K --posterior-integration-depth
    int posteriorIntegrationDepth;
    bool calculateMarginals;
    string algorithm;
    double RDF;             // -D --read-dependence-factor
    long double diffusionPriorScalar; // -V --diffusion-prior-scalar
    int WB;                      // -W --posterior-integration-bandwidth
    // XXX adjusting this to anything other than 1 may have bad consequences
    // for large numbers of samples
    int TB;                      // -Y --posterior-integration-depth
    bool includeMonoB;
    int TR;
    int I;
    long double minAltFraction;  // -F --min-alternate-fraction
    int minAltCount;             // -C --min-alternate-count
    int minAltTotal;             // -G --min-alternate-total
    int minCoverage;             // -! --min-coverage
    int maxCoverage;             // -+ --max-coverage
    int debuglevel;              // -d --debug increments
    bool debug; // set if debuglevel >=1
    bool debug2; // set if debuglevel >=2

    bool showReferenceRepeats;

    // functions
    Parameters(int argc, char** argv);
    void usage(char **argv);
    void simpleUsage(char **argv);

    // reporting
    string commandline;

};

#endif

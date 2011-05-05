#include "Parameters.h"
#include "convert.h"

using namespace std;

void Parameters::usage(char** argv) {
    cout 
         << "usage: " << argv[0] << " [OPTION] ... [BAM FILE] ... " << endl
         << endl
         << "Bayesian SNP and short INDEL polymorphism discovery." << endl
         << endl
         << "parameters:" << endl
         << endl
         << "   -h --help       Prints this help dialog." << endl
         << endl
         << "input and output:" << endl
         << endl
         << "   -b --bam FILE   Add FILE to the set of BAM files to be analyzed." << endl
         << "   -c --stdin      Read BAM input on stdin." << endl
         << "   -v --vcf FILE   Output VCF-format results to FILE." << endl
         << "   -f --fasta-reference FILE" << endl
         << "                   Use FILE as the reference sequence for analysis." << endl
         << "                   An index file (FILE.fai) will be created if none exists." << endl
         << "                   If neither --targets nor --region are specified, FreeBayes" << endl
         << "                   will analyze every position in this reference." << endl
         << "   -t --targets FILE" << endl
         << "                   Limit analysis to targets listed in the BED-format FILE." << endl
         << "   -r --region <chrom>:<start_position>..<end_position>" << endl
         << "                   Limit analysis to the specified region, 0-base coordinates," << endl
         << "                   end_position not included (same as BED format)." << endl
         << "   -s --samples FILE" << endl
         << "                   Limit analysis to samples listed (one per line) in the FILE." << endl
         << "                   By default FreeBayes will analyze all samples in its input" << endl
         << "                   BAM files." << endl
         << "   -A --cnv-map FILE" << endl
         << "                   Read a copy number map from the BED file FILE, which has" << endl
         << "                   the format:" << endl
         << "                      reference sequence, start, end, sample name, copy number" << endl
         << "                   ... for each region in each sample which does not have the" << endl
         << "                   default copy number as set by --ploidy." << endl
         << "   -L --trace FILE  Output an algorithmic trace to FILE." << endl
         << "   -l --failed-alleles FILE" << endl
         << "                   Write a BED file of the analyzed positions which do not" << endl
         << "                   pass --pvar to FILE." << endl
         << endl
         << "reporting:" << endl
         << endl
         << "   -P --pvar N     Report sites if the probability that there is a polymorphism" << endl
         << "                   at the site is greater than N.  default: 0.0001" << endl
         << "   -@ --report-all-alternates" << endl
         << "                   Report (in non-standard VCF format) each alternate allele" << endl
         << "                   at a site on its own line of VCF." << endl
         << "   -_ --show-reference-repeats" << endl
         << "                   Calculate and show information about reference repeats in" << endl
         << "                   the VCF output." << endl
         << endl
         << "population model:" << endl
         << endl
         << "   -T --theta N    The expected mutation rate or pairwise nucleotide diversity" << endl
         << "                   among the population under analysis.  This serves as the" << endl
         << "                   single parameter to the Ewens Sampling Formula prior model" << endl
         << "                   default: 0.001" << endl
         << "   -p --ploidy N   Sets the default ploidy for the analysis to N.  default: 2" << endl
         << "   -J --pooled     Assume that samples result from pooled sequencing." << endl
         << "                   When using this flag, set --ploidy to the number of" << endl
         << "                   alleles in each sample." << endl
         << endl
         << "reference allele:" << endl
         << endl
         << "   -Z --ignore-reference-allele" << endl
         << "                   By default, the reference allele is considered as another" << endl
         << "                   sample.  This flag excludes it from the analysis." << endl
         << "   -H --haploid-reference" << endl
         << "                   If using the reference sequence as a sample, consider it" << endl
         << "                   to be haploid.  default: false" << endl
         << "   --reference-quality MQ,BQ" << endl
         << "                   Assign mapping quality of MQ to the reference allele at each" << endl
         << "                   site and base quality of BQ.  default: 100,60" << endl
         << endl
         << "allele scope:" << endl
         << endl
         << "   -i --indels     Include insertion and deletion alleles in the analysis." << endl
         << "                   default: only analyze SNP alleles." << endl
         << "   -X --mnps       Include multi-nuceotide polymorphisms, MNPs, in the analysis." << endl
         << "                   default: only analyze SNP alleles." << endl
         << "   -I --no-snps    Ignore SNP alleles.  default: only analyze SNP alleles." << endl
         << "   -n --use-best-n-alleles N" << endl
         << "                   Evaluate only the best N SNP alleles, ranked by sum of" << endl
         << "                   supporting quality scores.  Set to 0 to use all.  default: 2" << endl
         << endl
         << "indel realignment:" << endl
         << endl
         << "   -O --left-align-indels" << endl
         << "                   Left-realign and merge gaps embedded in reads. default: false" << endl
         << endl
         << "input filters:" << endl
         << endl
         << "   -E --use-duplicate-reads" << endl
         << "                   Include duplicate-marked alignments in the analysis." << endl
         << "                   default: exclude duplicates" << endl
         << "   -m --min-mapping-quality Q" << endl
         << "                   Exclude alignments from analysis if they have a mapping" << endl
         << "                   quality less than Q.  default: 30" << endl
         << "   -q --min-base-quality Q" << endl
         << "                   Exclude alleles from analysis if their supporting base" << endl
         << "                   quality is less than Q.  default: 20" << endl
         << "   -R --min-supporting-quality MQ,BQ" << endl
         << "                   In order to consider an alternate allele, at least one supporting" << endl
         << "                   alignment must have mapping quality MQ, and one supporting " << endl
         << "                   allele must have base quality BQ. default: 0,0, unset" << endl
         << "   -Q --mismatch-base-quality-threshold Q" << endl
         << "                   Count mismatches toward --read-mismatch-limit if the base" << endl
         << "                   quality of the mismatch is >= Q.  default: 10" << endl
         << "   -U --read-mismatch-limit N" << endl
         << "                   Exclude reads with more than N mismatches where each mismatch" << endl
         << "                   has base quality >= mismatch-base-quality-threshold." << endl
         << "                   default: ~unbounded" << endl
         << "   -z --read-max-mismatch-fraction N" << endl
         << "                   Exclude reads with more than N [0,1] fraction of mismatches where" << endl
         << "                   each mismatch has base quality >= mismatch-base-quality-threshold" << endl
         << "                   default: 1.0" << endl
         << "   -$ --read-snp-limit N" << endl
         << "                   Exclude reads with more than N base mismatches, ignoring gaps" << endl
         << "                   with quality >= mismatch-base-quality-threshold." << endl
         << "                   default: ~unbounded" << endl
         << "   -e --read-indel-limit N" << endl
         << "                   Exclude reads with more than N separate gaps." << endl
         << "                   default: ~unbounded" << endl
         << "   -0 --no-filters Do not use any input base and mapping quality filters" << endl
         << "                   Equivalent to -m 0 -q 0 -R 0 -S 0" << endl
         << "   -x --indel-exclusion-window" << endl
         << "                   Ignore portions of alignments this many bases from a" << endl
         << "                   putative insertion or deletion allele.  default: 0" << endl
         << "   -F --min-alternate-fraction N" << endl
         << "                   Require at least this fraction of observations supporting" << endl
         << "                   an alternate allele within a single individual in the" << endl
         << "                   in order to evaluate the position.  default: 0.0" << endl
         << "   -C --min-alternate-count N" << endl
         << "                   Require at least this count of observations supporting" << endl
         << "                   an alternate allele within a single individual in order" << endl
         << "                   to evaluate the position.  default: 1" << endl
         << "   -3 --min-alternate-qsum N" << endl
         << "                   Require at least this sum of quality of observations supporting" << endl
         << "                   an alternate allele within a single individual in order" << endl
         << "                   to evaluate the position.  default: 0" << endl
         << "   -G --min-alternate-total N" << endl
         << "                   Require at least this count of observations supporting" << endl
         << "                   an alternate allele within the total population in order" << endl
         << "                   to use the allele in analysis.  default: 1" << endl
         << "   -! --min-coverage N" << endl
         << "                   Require at least this coverage to process a site.  default: 0" << endl
         << endl
         << "bayesian priors:" << endl
         << endl
         << "   -V --binomial-obs-priors" << endl
         << "                   Incorporate expectations about osbervations into the priors," << endl
         << "                   Uses read placement probability, strand balance probability," << endl
         << "                   and read position (5'-3') probability." << endl
         << "   -a --allele-balance-priors" << endl
         << "                   Use aggregate probability of observation balance between alleles" << endl
         << "                   as a component of the priors.  Best for observations with minimal" << endl
         << "                   inherent reference bias." << endl
         << "   -w --hwe-priors Use the probability of the combination arising under HWE given" << endl
         << "                   the allele frequency as estimated by observation frequency." << endl
         << "   -Y --no-ewens-priors" << endl
         << "                   Turns off the Ewens' Sampling Formula component of the priors." << endl
         << "   -k --no-population-priors" << endl
         << "                   Equivalent to --pooled --no-ewens-priors" << endl
         << endl
         << "algorithmic features:" << endl
         << endl
         << "   -M --expectation-maximization" << endl
         << "                   Use expectation maximization algorithm to integrate over posterior" << endl
         << "                   and obtain maximum likelihood genotype combination across all samples." << endl
         << "   -u --expectation-maximization-max-iterations N" << endl
         << "                   Iterate no more than this many times during expectation" << endl
         << "                   maximization step.  default: 5" << endl
         << "   -B --genotyping-max-iterations N" << endl
         << "                   Iterate no more than N times during genotyping step. default: 200." << endl
         << "   -W --posterior-integration-limits N,M" << endl
         << "                   Integrate all genotype combinations in our posterior space" << endl
         << "                   which include no more than N samples with their Mth best" << endl
         << "                   data likelihood. default: 1,3." << endl
         << "   -K --no-permute" << endl
         << "                   Do not scale prior probability of genotype combination given allele" << endl
         << "                   frequency by the number of permutations of included genotypes." << endl
         << "   -^ --genotype-combo-step-max N" << endl
         << "                   When generating genotype combinations, do not include genotypes" << endl
         << "                   where the genotype data likelihood is log(N) from the highest" << endl
         << "                   likelihood genotype for that individual.  default: ~unbounded" << endl
         << "   -N --exclude-unobserved-genotypes" << endl
         << "                   Skip sample genotypings for which the sample has no supporting reads." << endl
         << "   -S --genotype-variant-threshold N" << endl
         << "                   Limit posterior integration to samples where the second-best" << endl
         << "                   genotype likelihood is no more than log(N) from the highest" << endl
         << "                   genotype likelihood for the sample.  default: ~unbounded" << endl
         << "   -j --use-mapping-quality" << endl
         << "                   Use mapping quality of alleles when calculating data likelihoods." << endl
         << "   -D --read-dependence-factor N" << endl
         << "                   Incorporate non-independence of reads by scaling successive" << endl
         << "                   observations by this factor during data likelihood" << endl
         << "                   calculations.  default: 0.9" << endl
         << "   -= --no-marginals" << endl
         << "                   Do not calculate the marginal probability of genotypes.  Saves" << endl
         << "                   time and improves scaling performance in large populations." << endl
         << endl
         << "debugging:" << endl
         << endl
         << "   -d --debug      Print debugging output." << endl
         << "   -dd             Print more verbose debugging output (requires \"make DEBUG\")" << endl
         << endl
         << endl
         << "author:  Erik Garrison <erik.garrison@bc.edu>, Marth Lab, Boston College, 2010" << endl
         << "date:    " << FREEBAYES_COMPILE_DATE << endl
         << "version: " <<  FREEBAYES_VERSION << endl;
}


Parameters::Parameters(int argc, char** argv) {

    if (argc == 1) {
        usage(argv);
        exit(1);
    }

    // record command line parameters
    commandline = argv[0];
    for (int i = 1; i < argc; ++i) {
        commandline += " ";
        commandline += argv[i];
    }

    // set defaults

    // i/o parameters:
    useStdin = false;               // -c --stdin
    fasta = "";                // -f --fasta-reference
    targets = "";              // -t --targets
    region = "";               // -r --region
    samples = "";              // -s --samples
    cnvFile = "";
    output = "vcf";               // -v --vcf
    outputFile = "";
    traceFile = "";
    failedFile = "";

    // operation parameters
    outputAlleles = false;          //
    trace = false;                  // -L --trace
    useDuplicateReads = false;      // -E --use-duplicate-reads
    suppressOutput = false;         // -N --suppress-output
    reportAllAlternates = false;
    useBestNAlleles = 2;         // -n --use-best-n-alleles
    forceRefAllele = true;         // -Z --ignore-reference-allele
    useRefAllele = true;           // .....
    haploidReference = false;      // -H --haploid-reference
    allowIndels = false;            // -i --indels
    leftAlignIndels = false;       // -O --left-align-indels
    allowMNPs = false;            // -X --mnps
    allowSNPs = true;          // -I --no-snps
    pooled = false;                 // -J --pooled
    ewensPriors = true;
    permute = true;                // -K --permute
    useMappingQuality = false;
    obsBinomialPriors = false; // TODO
    hwePriors = false;
    alleleBalancePriors = false;
    excludeUnobservedGenotypes = false;
    excludePartiallyObservedGenotypes = false;
    genotypeVariantThreshold = 0;
    expectationMaximization = false;
    expectationMaximizationMaxIterations = 5;
    genotypingMaxIterations = 200;
    minPairedAltCount = 0;
    minAltMeanMapQ = 0;
    MQR = 100;                     // -M --reference-mapping-quality
    BQR = 60;                     // -B --reference-base-quality
    ploidy = 2;                  // -p --ploidy
    MQL0 = 30;                    // -m --min-mapping-quality
    BQL0 = 20;                    // -q --min-base-quality
    MQL1 = 0;                    // -R --min-supporting-quality MQ,BQ
    BQL1 = 0;                    //
    BQL2 = 10;                    // -Q --mismatch-base-quality-threshold
    RMU = 10000000;                     // -U --read-mismatch-limit
    readMaxMismatchFraction = 1.0;    //  -z --read-max-mismatch-fraction
    readSnpLimit = 10000000;       // -$ --read-snp-limit
    readIndelLimit = 10000000;     // -e --read-indel-limit
    genotypeComboStepMax = -1.0;                // -^ --genotype-combo-step-max
    IDW = -1;                     // -x --indel-exclusion-window
    TH = 10e-3;              // -T --theta
    PVL = 0.0001;             // -P --pvar
    RDF = 0.9;             // -D --read-dependence-factor
    diffusionPriorScalar = 1.0;     // -V --diffusion-prior-scalar
    WB = 1;                      // -W --posterior-integration-limits
    TB = 3;
    posteriorIntegrationDepth = 0;
    calculateMarginals = true;
    minAltFraction = 0.0;
    minAltCount = 1;
    minAltTotal = 1;
    minAltQSum = 0;
    minCoverage = 0;
    debuglevel = 0;
    debug = false;
    debug2 = false;

    showReferenceRepeats = false;

    int c; // counter for getopt

    static struct option long_options[] =
    {
        {"help", no_argument, 0, 'h'},
        {"bam", required_argument, 0, 'b'},
        {"stdin", no_argument, 0, 'c'},
        {"fasta-reference", required_argument, 0, 'f'},
        {"targets", required_argument, 0, 't'},
        {"region", required_argument, 0, 'r'},
        {"samples", required_argument, 0, 's'},
        {"cnv-map", required_argument, 0, 'A'},
        {"vcf", required_argument, 0, 'v'},
        {"trace", required_argument, 0, 'L'},
        {"failed-alleles", required_argument, 0, 'l'},
        {"use-duplicate-reads", no_argument, 0, 'E'},
        {"use-best-n-alleles", required_argument, 0, 'n'},
        {"ignore-reference-allele", no_argument, 0, 'Z'},
        {"haploid-reference", no_argument, 0, 'H'},
        {"no-filters", no_argument, 0, '0'},
        {"reference-quality", required_argument, 0, '1'},
        {"ploidy", required_argument, 0, 'p'},
        {"pooled", no_argument, 0, 'J'},
        {"no-ewens-priors", no_argument, 0, 'Y'}, // TODO
        {"no-population-priors", no_argument, 0, 'k'},
        {"use-mapping-quality", no_argument, 0, 'j'},
        {"min-mapping-quality", required_argument, 0, 'm'},
        {"min-base-quality", required_argument, 0, 'q'},
        {"min-supporting-quality", required_argument, 0, 'R'},
        {"mismatch-base-quality-threshold", required_argument, 0, 'Q'},
        {"read-mismatch-limit", required_argument, 0, 'U'},
        {"read-max-mismatch-fraction", required_argument, 0, 'z'},
        {"read-snp-limit", required_argument, 0, '$'},
        {"read-indel-limit", required_argument, 0, 'e'},
        {"genotype-combo-step-max", required_argument, 0, '^'},
        {"indels", no_argument, 0, 'i'},
        {"left-align-indels", no_argument, 0, 'O'},
        {"mnps", no_argument, 0, 'X'},
        {"no-snps", no_argument, 0, 'I'},
        {"indel-exclusion-window", required_argument, 0, 'x'},
        {"theta", required_argument, 0, 'T'},
        {"pvar", required_argument, 0, 'P'},
        {"read-dependence-factor", required_argument, 0, 'D'},
        {"binomial-obs-priors", no_argument, 0, 'V'},
        {"allele-balance-priors", no_argument, 0, 'a'},
        {"hwe-priors", no_argument, 0, 'w'},
        {"posterior-integration-limits", required_argument, 0, 'W'},
        {"min-alternate-fraction", required_argument, 0, 'F'},
        {"min-alternate-count", required_argument, 0, 'C'},
        //{"min-paired-alternate-count", required_argument, 0, 'Y'},
        {"min-alternate-total", required_argument, 0, 'G'},
        //{"min-alternate-mean-mapq", required_argument, 0, 'k'},
        {"min-alternate-qsum", required_argument, 0, '3'},
        {"min-coverage", required_argument, 0, '!'},
        {"no-permute", no_argument, 0, 'K'},
        {"no-marginals", no_argument, 0, '='},
        {"report-all-alternates", no_argument, 0, '@'},
        {"show-reference-repeats", no_argument, 0, '_'},
        {"exclude-unobserved-genotypes", no_argument, 0, 'N'},
        {"genotype-variant-threshold", required_argument, 0, 'S'},
        {"expectation-maximization", no_argument, 0, 'M'},
        {"expectation-maximization-max-iterations", required_argument, 0, 'u'},
        {"genotyping-max-iterations", required_argument, 0, 'B'},
        {"debug", no_argument, 0, 'd'},

        {0, 0, 0, 0}

    };

    while (true) {

        int option_index = 0;
        c = getopt_long(argc, argv, "hcOEZKjH0diNaI@_YkM=wVXJb:G:x:A:f:t:r:s:v:n:u:B:p:m:q:R:Q:U:$:e:T:P:D:^:S:W:F:C:L:l:z:1:3:",
                        long_options, &option_index);

        if (c == -1) // end of options
            break;

        switch (c) {

            // i/o parameters:
            // -b --bam
            case 'b':
                bams.push_back(optarg);
                break;

            // -c --stdin
            case 'c':
                useStdin = true;
                bams.push_back("stdin");
                break;

            // -f --fasta-reference
            case 'f':
                fasta = optarg;
                break;

            // -t --targets
            case 't':
                targets = optarg;
                break;

            // -r --region
            case 'r':
                region = optarg;
                break;

            // -s --samples
            case 's':
                samples = optarg;
                break;

            // -A --cnv-file
            case 'A':
                cnvFile = optarg;
                break;

            // -j --use-mapping-quality
            case 'j':
                useMappingQuality = true;
                break;

            // -v --vcf
            case 'v':
                output = "vcf";
                outputFile = optarg;
                break;

            // -O --left-align-indels
            case 'O':
                leftAlignIndels = true;
                break;

            // -L --trace
            case 'L':
                traceFile = optarg;
                trace = true;
                break;

            // -l --failed-alleles
            case 'l':
                failedFile = optarg;
                break;

            // -E --use-duplicate-reads
            case 'E':
                useDuplicateReads = true;
                break;

            // -3 --min-alternate-qsum
            case '3':
                if (!convert(optarg, minAltQSum)) {
                    cerr << "could not parse min-alternate-qsum" << endl;
                    exit(1);
                }
                break;

            // -G --min-alternate-total
            case 'G':
                if (!convert(optarg, minAltTotal)) {
                    cerr << "could not parse min-alternate-total" << endl;
                    exit(1);
                }
                break;

            // -! --min-coverage
            case '!':
                if (!convert(optarg, minCoverage)) {
                    cerr << "could not parse min-coverage" << endl;
                    exit(1);
                }
                break;

            // -n --use-best-n-alleles
            case 'n':
                if (!convert(optarg, useBestNAlleles)) {
                    cerr << "could not parse use-best-n-alleles" << endl;
                    exit(1);
                }
                break;

            // -Z --ignore-reference-allele
            case 'Z':
                forceRefAllele = false;
                useRefAllele = false;
                break;

            // -H --haploid-reference
            case 'H':
                haploidReference = true;
                break;

            // -0 --no-filters
            case '0':
                MQL0 = 0;
                MQL1 = 0;
                BQL0 = 0;
                BQL1 = 0;
                break;

            // -M --expectation-maximization
            case 'M':
                expectationMaximization = true;
                break;

            // -u --expectation-maximization-max-iterations
            case 'u':
                if (!convert(optarg, expectationMaximizationMaxIterations)) {
                    cerr << "could not parse expectation-maximization-max-iterations" << endl;
                    exit(1);
                }
                break;

            // -B --genotyping-max-iterations
            case 'B':
                if (!convert(optarg, genotypingMaxIterations)) {
                    cerr << "could not parse genotyping-max-iterations" << endl;
                    exit(1);
                }
                break;

            // -1 --reference-quality
            case '1':
                if (!convert(split(optarg, ",").front(), MQR)) {
                    cerr << "could not parse reference mapping quality" << endl;
                    exit(1);
                }
                if (!convert(split(optarg, ",").back(), BQR)) {
                    cerr << "could not parse reference base quality" << endl;
                    exit(1);
                }
                break;

            // -p --ploidy
            case 'p':
                if (!convert(optarg, ploidy)) {
                    cerr << "could not parse ploidy" << endl;
                    exit(1);
                }
                if (ploidy <= 0) {
                    cerr << "cannot set ploidy to less than 1" << endl;
                    exit(1);
                }
                break;

            // -J --pooled
            case 'J':
                pooled = true;
                break;

            // -m --min-mapping-quality
            case 'm':
                if (!convert(optarg, MQL0)) {
                    cerr << "could not parse min-base-quality" << endl;
                    exit(1);
                }
                break;

            // -q --min-base-quality
            case 'q':
                if (!convert(optarg, BQL0)) {
                    cerr << "could not parse min-base-quality" << endl;
                    exit(1);
                }
                break;

            // -R --min-supporting-quality
            case 'R':
                if (!convert(split(optarg, ",").front(), MQL1)) {
                    cerr << "could not parse min-supporting-quality MQ" << endl;
                    exit(1);
                }
                if (!convert(split(optarg, ",").back(), BQL1)) {
                    cerr << "could not parse min-supporting-quality MQ" << endl;
                    exit(1);
                }
                break;

            // -N --exclude-unobserved-genotypes
            case 'N':
                excludeUnobservedGenotypes = true;
                break;

            // -S --genotype-variant-threshold
            case 'S':
                if (!convert(optarg, genotypeVariantThreshold)) {
                    cerr << "could not parse genotype-variant-threshold" << endl;
                    exit(1);
                }
                break;

            // -Q --mismatch-base-quality-threshold
            case 'Q':
                if (!convert(optarg, BQL2)) {
                    cerr << "could not parse mismatch-base-quality-threshold" << endl;
                    exit(1);
                }
                break;

            // -U --read-mismatch-limit
            case 'U':
                if (!convert(optarg, RMU)) {
                    cerr << "could not parse read-mismatch-limit" << endl;
                    exit(1);
                }
                break;

            // -z --read-max-mismatch-fraction
            case 'z':
                if (!convert(optarg, readMaxMismatchFraction)) {
                    cerr << "could not parse read-mismatch-limit" << endl;
                    exit(1);
                }
                break;

            // -$ --read-snp-limit
            case '$':
                if (!convert(optarg, readSnpLimit)) {
                    cerr << "could not parse read-snp-limit" << endl;
                    exit(1);
                }
                break;

            // -e --read-indel-limit
            case 'e':
                if (!convert(optarg, readIndelLimit)) {
                    cerr << "could not parse read-indel-limit" << endl;
                    exit(1);
                }
                break;

            case '^':
                if (!convert(optarg, genotypeComboStepMax)) {
                    cerr << "could not parse genotype-combo-step-max" << endl;
                    exit(1);
                }
                break;

            // -x --indel-exclusion-window
            case 'x':
                if (!convert(optarg, IDW)) {
                    cerr << "could not parse indel-exclusion-window" << endl;
                    exit(1);
                }
                break;

            // -i --indels
            case 'i':
                allowIndels = true;
                reportAllAlternates = true;
                break;

            // -X --mnps
            case 'X':
                allowMNPs = true;
                break;

            // -I --no-snps
            case 'I':
                allowSNPs = false;
                break;

            // -T --theta
            case 'T':
                if (!convert(optarg, TH)) {
                    cerr << "could not parse theta" << endl;
                    exit(1);
                }
                break;

            // -P --pvar
            case 'P':
                if (!convert(optarg, PVL)) {
                    cerr << "could not parse pvar" << endl;
                    exit(1);
                }
                break;

            // -D --read-dependence-factor
            case 'D':
                if (!convert(optarg, RDF)) {
                    cerr << "could not parse read-dependence-factor" << endl;
                    exit(1);
                }
                break;

            // observation priors
            case 'V':
                obsBinomialPriors = true;
                break;

            // allele balance
            case 'a':
                alleleBalancePriors = true;
                break;

            case 'w':
                hwePriors = true;
                break;

            // -W --posterior-integration-limits
            case 'W':
                if (!convert(split(optarg, ",").front(), WB)) {
                    cerr << "could not parse posterior-integration-limits (bandwidth)" << endl;
                    exit(1);
                }
                if (!convert(split(optarg, ",").back(), TB)) {
                    cerr << "could not parse posterior-integration-limits (banddepth)" << endl;
                    exit(1);
                }
                break;

            // -F --min-alternate-fraction
            case 'F':
                if (!convert(optarg, minAltFraction)) {
                    cerr << "could not parse min-alternate-fraction" << endl;
                    exit(1);
                }
                break;

            // -C --min-alternate-count
            case 'C':
                if (!convert(optarg, minAltCount)) {
                    cerr << "could not parse min-alternate-count" << endl;
                    exit(1);
                }
                break;

            // -Y --no-ewens-priors
            case 'Y':
                ewensPriors = false;
                break;

            // -k --no-population-priors
            case 'k':
                pooled = true;
                ewensPriors = false;
                break;

            // -K --no-permute
            case 'K':
                permute = false;
                break;

            case '=':
                calculateMarginals = false;
                break;

            case '@':
                reportAllAlternates = true;
                break;

            case '_':
                showReferenceRepeats = true;
                break;

            // -d --debug
            case 'd':
                ++debuglevel;
                break;

            case 'h':
                usage(argv);
                exit(0);
                break;
 
            case '?': // print a suggestion about the most-likely long option which the argument matches
                {
                    string bad_arg(argv[optind - 1]);
                    option* opt = &long_options[0];
                    option* closest_opt = opt;
                    int shortest_distance = levenshteinDistance(opt->name, bad_arg);
                    ++opt;
                    while (opt->name != 0) {
                        int distance = levenshteinDistance(opt->name, bad_arg);
                        if (distance < shortest_distance) {
                            shortest_distance = distance;
                            closest_opt = opt;
                        }
                        ++opt;
                    }
                    cerr << "did you mean --" << closest_opt->name << " ?" << endl;
                    exit(1);
                }
                break;

            default:
                abort ();
        }

    }
 
    // any remaining arguments are considered as bam files
    if (optind < argc) {
        if (useStdin) {
            cerr << "--stdin flag specified, but a list of BAM files given.  Jumping disabled." << endl;
        }
        while (optind < argc) {
            bams.push_back(argv[optind++]);
        }
    }

    if (debuglevel >= 1) {
        debug = true;
    }
    if (debuglevel >= 2) {
        debug2 = true;
    }

    if (bams.size() == 0) {
        cerr << "Please specify a BAM file or files." << endl;
        exit(1);
    }

    if (fasta == "") {
        cerr << "Please specify a fasta reference file." << endl;
        exit(1);
    }

}

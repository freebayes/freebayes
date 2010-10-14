#include "Parameters.h"
#include "convert.h"

using namespace std;

void Parameters::usage(char** argv) {
    cout 
         << "usage: " << argv[0] << " [OPTION] ... [BAM FILE] ... " << endl
         << endl
         << "Bayesian SNP and short INDEL polymorphism discovery." << endl
         << endl
         << "options:" << endl
         << endl
         << "   -h --help       Prints this help dialog." << endl
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
         << "                   Limit analysis to the specified region." << endl
         << "   -s --samples FILE" << endl
         << "                   Limit analysis to samples listed (one per line) in the FILE." << endl
         << "                   By default FreeBayes will analyze all samples in its input" << endl
         << "                   BAM files." << endl
         << "   -j --json       Toggle JSON output of results on stdout." << endl
         << "   -O --output-alleles" << endl
         << "                   When --json is set, add records to the JSON output stream" << endl
         << "                   describing each allele in the input." << endl
         << "   -L --trace FILE  Output an algorithmic trace to FILE." << endl
         << "   -X --suppress-output" << endl
         << "                   Suppress output.  Used for debugging." << endl
         << "   -P --pvar N     Report sites if the probability that there is a polymorphism" << endl
         << "                   at the site is greater than N.  default: 0.0" << endl
         << "   -G --factorial-data-likelihoods" << endl
         << "                   Use the factorial/gamma distribution for calculating the" << endl
         << "                   sampling component of data likelihoods.  Default is to use" << endl
         << "                   the multinomial distribution." << endl
         << "   -T --theta N    The expected mutation rate or pairwise nucleotide diversity" << endl
         << "                   among the population under analysis.  This serves as the" << endl
         << "                   single parameter to the Ewens Sampling Formula prior model" << endl
         << "                   default: 0.001" << endl
         << "   -p --ploidy N   Sets the default ploidy for the analysis to N.  default: 2" << endl
         << "   -i --allow-indels" << endl
         << "                   Include insertion and deletion alleles in the analysis." << endl
         << "                   default: only analyze SNP alleles." << endl
         << "   -n --use-best-n-alleles N" << endl
         << "                   Evaluate only the best N alleles, ranked by sum of" << endl
         << "                   supporting quality scores.  By default, evaluate all" << endl
         << "                   possible alleles." << endl
         << "   -E --use-duplicate-reads" << endl
         << "                   Include duplicate-marked alignments in the analysis." << endl
         << "                   default: exclude duplicates" << endl
         << "   -A --use-reference-allele" << endl
         << "                   Include the reference allele in the analysis as if it is an" << endl
         << "                   additional sample." << endl
         << "   -M --reference-mapping-quality Q" << endl
         << "                   Assign mapping quality of Q to the reference allele at each" << endl
         << "                   site.  default: 100" << endl
         << "   -B --reference-base-quality Q" << endl
         << "                   Assign a base quality of Q to the reference allele at each" << endl
         << "                   site.  default: 60" << endl
         << "   -Z --force-reference-allele" << endl
         << "                   Force the use of the reference allele in analysis, even if" << endl
         << "                   it isn't one of the best N alleles provided its assigned" << endl
         << "                   base quality." << endl
         << "   -m --min-mapping-quality Q" << endl
         << "                   Exclude alignments from analysis if they have a mapping" << endl
         << "                   quality less than Q.  default: 30" << endl
         << "   -q --min-base-quality Q" << endl
         << "                   Exclude alleles from analysis if their supporting base" << endl
         << "                   quality is less than Q.  default: 20" << endl
         << "   -R --min-supporting-mapping-quality Q" << endl
         << "                   In order to consider an alternate allele, at least one" << endl
         << "                   supporting alignment must have this mapping quality." << endl
         << "                   default: 40" << endl
         << "   -S --min-supporting-base-quality Q" << endl
         << "                   In order to consider an alternate allele, at least one" << endl
         << "                   supporting alignment must have this base quality at the" << endl
         << "                   site of the allele.  default: 40" << endl
         << "   -U --read-mismatch-limit N" << endl
         << "                   The maximum allowable number of mismatches between an" << endl
         << "                   alignment and the reference where each mismatch has" << endl
         << "                   base quality >= the mismatch-base-quality-threshold." << endl
         << "                   default: ~unbounded" << endl
         << "   -Q --mismatch-base-quality-threshold Q" << endl
         << "                   Count mismatches toward --read-mismatch-limit if the base" << endl
         << "                   quality of the mismatch is >= Q.  default: 10" << endl
         << "   -x --indel-exclusion-window" << endl
         << "                   Ignore portions of alignments this many bases from a" << endl
         << "                   putative insertion or deletion allele.  default: 0" << endl
         << "   -D --read-dependence-factor N" << endl
         << "                   Incorporate non-independence of reads by scaling successive" << endl
         << "                   observations by this factor during data likelihood" << endl
         << "                   calculations.  default: 0.9" << endl
         << "   -W --posterior-integration-bandwidth N" << endl
         << "                   Integrate all genotype combinations in our posterior space" << endl
         << "                   which lie no more than N steps from the most likely" << endl
         << "                   combination in terms of data likelihoods, taking the N" << endl
         << "                   steps from the most to least likely genotype for each" << endl
         << "                   individual.  default: 2" << endl
         << "   -F --min-alternate-fraction N" << endl
         << "                   Require at least this fraction of observations supporting" << endl
         << "                   an alternate allele within a single individual in the" << endl
         << "                   in order to evaluate the position.  default: 0.0" << endl
         << "   -C --min-alternate-count N" << endl
         << "                   Require at least this count of observations supporting" << endl
         << "                   an alternate allele within a single individual in order" << endl
         << "                   to evaluate the position.  default: 1" << endl
         << "   -d --debug      Print debugging output." << endl
         << "   -dd             Print more verbose debugging output" << endl
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

    // set defaults

    // i/o parameters:
    useStdin = false;               // -c --stdin
    fasta = "";                // -f --fasta-reference
    targets = "";              // -t --targets
    region = "";               // -r --region
    samples = "";              // -s --samples
    output = "vcf";               // -v --vcf
    outputFile = "";
    traceFile = "";

    // operation parameters
    outputAlleles = false;          // -O --output-alleles
    trace = false;                  // -L --trace
    useDuplicateReads = false;      // -E --use-duplicate-reads
    suppressOutput = false;         // -S --suppress-output
    bamBayesDataLikelihoods = false;// -G --factorial-data-likelihoods
    useBestNAlleles = 0;         // -n --use-best-n-alleles
    forceRefAllele = false;         // -F --force-reference-allele
    useRefAllele = false;           // -U --use-reference-allele
    allowIndels = false;            // -i --allow-indels
    MQR = 100;                     // -M --reference-mapping-quality
    BQR = 60;                     // -B --reference-base-quality
    ploidy = 2;                  // -p --ploidy
    MQL0 = 30;                    // -m --min-mapping-quality
    BQL0 = 20;                    // -q --min-base-quality
    MQL1 = 40;                    // -R --min-supporting-mapping-quality
    BQL1 = 30;                    // -S --min-supporting-base-quality
    BQL2 = 10;                    // -Q --mismatch-base-quality-threshold
    RMU = 10000000;                     // -U --read-mismatch-limit
    IDW = -1;                     // -x --indel-exclusion-window
    TH = 10e-3;              // -T --theta
    PVL = 0.0;             // -P --pvar
    RDF = 0.9;             // -D --read-dependence-factor
    WB = 2;                      // -W --posterior-integration-bandwidth
    minAltFraction = 0.0;
    minAltCount = 1;
    debuglevel = 0;
    debug = false;
    debug2 = false;


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
        {"json", no_argument, 0, 'v'},
        {"vcf", required_argument, 0, 'v'},
        {"output-alleles", no_argument, 0, 'O'},
        {"trace", required_argument, 0, 'L'},
        {"use-duplicate-reads", no_argument, 0, 'E'},
        {"suppress-output", no_argument, 0, 'X'},
        {"factorial-data-likelihoods", no_argument, 0, 'G'},
        {"use-best-n-alleles", required_argument, 0, 'n'},
        {"force-reference-allele", no_argument, 0, 'Z'},
        {"use-reference-allele", no_argument, 0, 'A'},
        {"reference-mapping-quality", required_argument, 0, 'M'},
        {"reference-base-quality", required_argument, 0, 'B'},
        {"ploidy", required_argument, 0, 'p'},
        {"min-mapping-quality", required_argument, 0, 'm'},
        {"min-base-quality", required_argument, 0, 'q'},
        {"min-supporting-mapping-quality", required_argument, 0, 'R'},
        {"min-supporting-base-quality", required_argument, 0, 'S'},
        {"mismatch-base-quality-threshold", required_argument, 0, 'Q'},
        {"read-mismatch-limit", required_argument, 0, 'U'},
        {"allow-indels", no_argument, 0, 'i'},
        {"indel-exclusion-window", required_argument, 0, 'x'},
        {"theta", required_argument, 0, 'T'},
        {"pvar", required_argument, 0, 'P'},
        {"read-dependence-factor", required_argument, 0, 'D'},
        {"posterior-integration-bandwidth", required_argument, 0, 'W'},
        {"min-alternate-fraction", required_argument, 0, 'F'},
        {"min-alternate-count", required_argument, 0, 'C'},
        {"debug", no_argument, 0, 'd'},

        {0, 0, 0, 0}

    };

    while (true) {

        int option_index = 0;
        c = getopt_long(argc, argv, "hcOEXGZAdDb:f:t:r:s:v:j:n:M:B:p:m:q:R:S:Q:U:I:T:P:D:W:F:C:",
                        long_options, &option_index);

        if (c == -1) // end of options
            break;

        switch (c) {

            // i/o parameters:
            // -b --bam
            case 'b':
                bams.push_back(optarg);
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

            // -j --json
            case 'j':
                output = "json";
                break;

            // -v --vcf
            case 'v':
                output = "vcf";
                outputFile = optarg;
                break;

            // -O --output-alleles
            case 'O':
                outputAlleles = true;
                break;

            // -L --trace
            case 'L':
                traceFile = optarg;
                trace = true;
                break;

            // -E --use-duplicate-reads
            case 'E':
                useDuplicateReads = true;
                break;

            // -X --suppress-output
            case 'X':
                suppressOutput = true;
                break;

            // -G --factorial-data-likelihoods
            case 'G':
                bamBayesDataLikelihoods = true;
                break;

            // -n --use-best-n-alleles
            case 'n':
                if (!convert(optarg, useBestNAlleles)) {
                    cerr << "could not parse reference-mapping-quality" << endl;
                    exit(1);
                }
                break;

            // -Z --force-reference-allele
            case 'Z':
                forceRefAllele = true;
                break;

            // -A --use-reference-allele
            case 'A':
                useRefAllele = true;
                break;

            // -M --reference-mapping-quality
            case 'M':
                if (!convert(optarg, MQR)) {
                    cerr << "could not parse reference-mapping-quality" << endl;
                    exit(1);
                }
                break;

            // -B --reference-base-quality
            case 'B':
                if (!convert(optarg, BQR)) {
                    cerr << "could not parse reference-base-quality" << endl;
                    exit(1);
                }
                break;

            // -p --ploidy
            case 'p':
                if (!convert(optarg, ploidy)) {
                    cerr << "could not parse ploidy" << endl;
                    exit(1);
                }
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

            // -R --min-supporting-mapping-quality
            case 'R':
                if (!convert(optarg, MQL1)) {
                    cerr << "could not parse min-supporting-mapping-quality" << endl;
                    exit(1);
                }
                break;

            // -S --min-supporting-base-quality
            case 'S':
                if (!convert(optarg, BQL1)) {
                    cerr << "could not parse min-supporting-base-quality" << endl;
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

            // -x --indel-exclusion-window
            case 'x':
                if (!convert(optarg, IDW)) {
                    cerr << "could not parse indel-exclusion-window" << endl;
                    exit(1);
                }
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

            // -W --posterior-integration-bandwidth
            case 'W':
                if (!convert(optarg, WB)) {
                    cerr << "could not parse posterior-integration-bandwidth" << endl;
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

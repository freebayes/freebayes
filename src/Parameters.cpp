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
         << "                   Limit analysis to the specified region, 1-base coordinates," << endl
         << "                   end_position inclusive." << endl
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
         << "   -j --json       Toggle JSON output of results on stdout." << endl
         << "   -O --output-alleles" << endl
         << "                   When --json is set, add records to the JSON output stream" << endl
         << "                   describing each allele in the input." << endl
         << "   -@ --report-all-alternates" << endl
         << "                   Report (in non-standard VCF format) each alternate allele" << endl
         << "                   at a site on its own line of VCF." << endl
         << "   -L --trace FILE  Output an algorithmic trace to FILE." << endl
         << "   -l --failed-alleles FILE" << endl
         << "                   Write a BED file of the analyzed positions which do not" << endl
         << "                   pass --pvar to FILE." << endl
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
         << "   -J --pooled     Assume that samples result from pooled sequencing." << endl
         << "                   When using this flag, set --ploidy to the number of" << endl
         << "                   alleles in each sample." << endl
         << "   -i --indels     Include insertion and deletion alleles in the analysis." << endl
         << "                   default: only analyze SNP alleles." << endl
         << "   -X --mnps       Include multi-nuceotide polymorphisms, MNPs, in the analysis." << endl
         << "                   default: only analyze SNP alleles." << endl
         << "   -I --no-snps    Ignore SNP alleles.  default: only analyze SNP alleles." << endl
         << "   -n --use-best-n-alleles N" << endl
         << "                   Evaluate only the best N alleles, ranked by sum of" << endl
         << "                   supporting quality scores.  default: 2" << endl
         << "   -N --use-all-alleles" << endl
         << "                   Evaluate all possible alleles." << endl
         << "   -E --use-duplicate-reads" << endl
         << "                   Include duplicate-marked alignments in the analysis." << endl
         << "                   default: exclude duplicates" << endl
         << "   -M --reference-mapping-quality Q" << endl
         << "                   Assign mapping quality of Q to the reference allele at each" << endl
         << "                   site.  default: 100" << endl
         << "   -B --reference-base-quality Q" << endl
         << "                   Assign a base quality of Q to the reference allele at each" << endl
         << "                   site.  default: 60" << endl
         << "   -Z --ignore-reference-allele" << endl
         << "                   By default, the reference allele is considered as another" << endl
         << "                   sample.  This flag excludes it from the analysis." << endl
         << "   -H --haploid-reference" << endl
         << "                   If using the reference sequence as a sample, consider it" << endl
         << "                   to be haploid.  default: false" << endl
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
         << "                   site of the allele.  default: 30" << endl
         << "   -U --read-mismatch-limit N" << endl
         << "                   The maximum allowable number of mismatches between an" << endl
         << "                   alignment and the reference where each mismatch has" << endl
         << "                   base quality >= the mismatch-base-quality-threshold." << endl
         << "                   default: ~unbounded" << endl
         << "   -Q --mismatch-base-quality-threshold Q" << endl
         << "                   Count mismatches toward --read-mismatch-limit if the base" << endl
         << "                   quality of the mismatch is >= Q.  default: 10" << endl
         << "   -0 --no-filters Do not use any input base and mapping quality filters" << endl
         << "                   Equivalent to -m 0 -q 0 -R 0 -S 0" << endl
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
         << "   -Y --posterior-integration-banddepth N" << endl
         << "                   Generate all genotype combinations for which up to this" << endl
         << "                   number of samples have up to their -W'th worst genotype" << endl
         << "                   according to data likelihood.  default: 2" << endl
         << "   -K --posterior-integration-depth N" << endl
         << "                   Keep this many genotype combinations for calculating genotype" << endl
         << "                   marginal probabilities for each sample and overall variant" << endl
         << "                   quality.  Default behavior is to keep all.  For the default" << endl
         << "                   value of -W (2) each variant site will require 3N^2 calculations" << endl
         << "                   to establish marginal genotype probabilities where N is the" << endl
         << "                   number of individuals, so setting a sensible number here" << endl
         << "                   (0.1N to 0.25N) will help improve scalability in large" << endl
         << "                   datasets at the cost of accuracy in calcuating marginal" << endl
         << "                   genotype probabilites." << endl
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
    outputAlleles = false;          // -O --output-alleles
    trace = false;                  // -L --trace
    useDuplicateReads = false;      // -E --use-duplicate-reads
    suppressOutput = false;         // -N --suppress-output
    reportAllAlternates = false;
    bamBayesDataLikelihoods = false;// -G --factorial-data-likelihoods
    useBestNAlleles = 2;         // -n --use-best-n-alleles
    forceRefAllele = true;         // -Z --ignore-reference-allele
    useRefAllele = true;           // .....
    haploidReference = false;      // -H --haploid-reference
    allowIndels = false;            // -i --indels
    allowMNPs = false;            // -X --mnps
    allowSNPs = true;          // -I --no-snps
    pooled = false;                 // -J --pooled
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
    TB = 1;                 // -Y --posterior-integration-banddepth
    posteriorIntegrationDepth = 0;
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
        {"cnv-map", required_argument, 0, 'A'},
        {"json", no_argument, 0, 'j'},
        {"vcf", required_argument, 0, 'v'},
        {"output-alleles", no_argument, 0, 'O'},
        {"trace", required_argument, 0, 'L'},
        {"failed-alleles", required_argument, 0, 'l'},
        {"use-duplicate-reads", no_argument, 0, 'E'},
        {"factorial-data-likelihoods", no_argument, 0, 'G'},
        {"use-best-n-alleles", required_argument, 0, 'n'},
        {"use-all-alleles", no_argument, 0, 'N'},
        {"ignore-reference-allele", no_argument, 0, 'Z'},
        {"haploid-reference", no_argument, 0, 'H'},
        {"no-filters", no_argument, 0, '0'},
        {"reference-mapping-quality", required_argument, 0, 'M'},
        {"reference-base-quality", required_argument, 0, 'B'},
        {"ploidy", required_argument, 0, 'p'},
        {"pooled", no_argument, 0, 'J'},
        {"min-mapping-quality", required_argument, 0, 'm'},
        {"min-base-quality", required_argument, 0, 'q'},
        {"min-supporting-mapping-quality", required_argument, 0, 'R'},
        {"min-supporting-base-quality", required_argument, 0, 'S'},
        {"mismatch-base-quality-threshold", required_argument, 0, 'Q'},
        {"read-mismatch-limit", required_argument, 0, 'U'},
        {"indels", no_argument, 0, 'i'},
        {"mnps", no_argument, 0, 'X'},
        {"no-snps", no_argument, 0, 'I'},
        {"indel-exclusion-window", required_argument, 0, 'x'},
        {"theta", required_argument, 0, 'T'},
        {"pvar", required_argument, 0, 'P'},
        {"read-dependence-factor", required_argument, 0, 'D'},
        {"posterior-integration-bandwidth", required_argument, 0, 'W'},
        {"min-alternate-fraction", required_argument, 0, 'F'},
        {"min-alternate-count", required_argument, 0, 'C'},
        {"posterior-integration-depth", required_argument, 0, 'K'},
        {"report-all-alternates", no_argument, 0, '@'},
        {"debug", no_argument, 0, 'd'},

        {0, 0, 0, 0}

    };

    while (true) {

        int option_index = 0;
        c = getopt_long(argc, argv, "hcOENGZH0dDiI@XJb:x:A:f:t:r:s:v:j:n:M:B:p:m:q:R:S:Q:U:T:P:D:W:F:C:K:Y:L:l:",
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

            // -l --failed-alleles
            case 'l':
                failedFile = optarg;
                break;

            // -E --use-duplicate-reads
            case 'E':
                useDuplicateReads = true;
                break;

            // -N --use-all-alleles
            case 'N':
                useBestNAlleles = 0;
                break;

            // -G --factorial-data-likelihoods
            case 'G':
                bamBayesDataLikelihoods = true;
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

            // -i --indels
            case 'i':
                allowIndels = true;
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

            // -W --posterior-integration-bandwidth
            case 'W':
                if (!convert(optarg, WB)) {
                    cerr << "could not parse posterior-integration-bandwidth" << endl;
                    exit(1);
                }
                break;

            // -Y --posterior-integration-banddepth
            case 'Y':
                if (!convert(optarg, TB)) {
                    cerr << "could not parse posterior-integration-banddepth" << endl;
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

            // -K --posterior-marginal-depth
            case 'K':
                if (!convert(optarg, posteriorIntegrationDepth)) {
                    cerr << "could not parse posterior-integration-depth" << endl;
                    exit(1);
                }
                break;

            case '@':
                reportAllAlternates = true;
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

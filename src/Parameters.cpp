#include "Parameters.h"
#include "convert.h"

using namespace std;

void Parameters::simpleUsage(char ** argv) {
    cout
        << "usage: " << argv[0] << " -f [REFERENCE] [OPTIONS] [BAM FILES] >[OUTPUT]" << endl
        << endl
        << "Bayesian haplotype-based polymorphism discovery." << endl
        << endl
        << "parameters:" << endl
        << endl
        << "   -h --help       For a complete description of options." << endl
        << endl
        << "citation: Erik Garrison, Gabor Marth" << endl
        << "          \"Haplotype-based variant detection from short-read sequencing\"" << endl
        << "          arXiv:1207.3907 (http://arxiv.org/abs/1207.3907)" << endl
        << endl
        << "author:   Erik Garrison <erik.garrison@bc.edu>, Marth Lab, Boston College, 2010-2014" << endl
        << "version:  " << VERSION_GIT << endl;

}

void Parameters::usage(char** argv) {
    cout 
        << "usage: " << argv[0] << " [OPTION] ... [BAM FILE] ... " << endl
        << endl
        << "Bayesian haplotype-based polymorphism discovery." << endl
        << endl
        << "citation: Erik Garrison, Gabor Marth" << endl
        << "          \"Haplotype-based variant detection from short-read sequencing\"" << endl
        << "          arXiv:1207.3907 (http://arxiv.org/abs/1207.3907)" << endl
        << endl
        << "overview:" << endl
        << endl
        << "    To call variants from aligned short-read sequencing data, supply BAM files and" << endl
        << "    a reference.  FreeBayes will provide VCF output on standard out describing SNPs," << endl
        << "    indels, and complex variants in samples in the input alignments." << endl
        << endl
        << "    By default, FreeBayes will consider variants supported by at least 2" << endl
        << "    observations in a single sample (-C) and also by at least 20% of the reads from" << endl
        << "    a single sample (-F).  These settings are suitable to low to high depth" << endl
        << "    sequencing in haploid and diploid samples, but users working with polyploid or" << endl
        << "    pooled samples may wish to adjust them depending on the characteristics of" << endl
        << "    their sequencing data." << endl
        << endl
        << "    FreeBayes is capable of calling variant haplotypes shorter than a read length" << endl
        << "    where multiple polymorphisms segregate on the same read.  The maximum distance" << endl
        << "    between polymorphisms phased in this way is determined by the" << endl
        << "    --max-complex-gap, which defaults to 3bp.  In practice, this can comfortably be" << endl
        << "    set to half the read length." << endl
        << endl
        << "    Ploidy may be set to any level (-p), but by default all samples are assumed to" << endl
        << "    be diploid.  FreeBayes can model per-sample and per-region variation in" << endl
        << "    copy-number (-A) using a copy-number variation map." << endl
        << endl
        << "    FreeBayes can act as a frequency-based pooled caller and describe variants" << endl
        << "    and haplotypes in terms of observation frequency rather than called genotypes." << endl
        << "    To do so, use --pooled-continuous and set input filters to a suitable level." << endl
        << "    Allele observation counts will be described by AO and RO fields in the VCF output." << endl
        << endl
        << endl
        << "examples:" << endl
        << endl
        << "    # call variants assuming a diploid sample" << endl
        << "    freebayes -f ref.fa aln.bam >var.vcf" << endl
        << endl
        << "    # call variants assuming a diploid sample, providing gVCF output" << endl
        << "    freebayes -f ref.fa --gvcf aln.bam >var.gvcf" << endl
        << endl
        << "    # require at least 5 supporting observations to consider a variant" << endl
        << "    freebayes -f ref.fa -C 5 aln.bam >var.vcf" << endl
        << endl
        << "    # use a different ploidy" << endl
        << "    freebayes -f ref.fa -p 4 aln.bam >var.vcf" << endl
        << endl
        << "    # assume a pooled sample with a known number of genome copies" << endl
        << "    freebayes -f ref.fa -p 20 --pooled-discrete aln.bam >var.vcf" << endl
        << endl
        << "    # generate frequency-based calls for all variants passing input thresholds" << endl
        << "    freebayes -f ref.fa -F 0.01 -C 1 --pooled-continuous aln.bam >var.vcf" << endl
        << endl
        << "    # use an input VCF (bgzipped + tabix indexed) to force calls at particular alleles" << endl
        << "    freebayes -f ref.fa -@ in.vcf.gz aln.bam >var.vcf" << endl
        << endl
        << "    # generate long haplotype calls over known variants" << endl
        << "    freebayes -f ref.fa --haplotype-basis-alleles in.vcf.gz \\ " << endl
        << "                        --haplotype-length 50 aln.bam" << endl
        << endl
        << "    # naive variant calling: simply annotate observation counts of SNPs and indels" << endl
        << "    freebayes -f ref.fa --haplotype-length 0 --min-alternate-count 1 \\ " << endl
        << "        --min-alternate-fraction 0 --pooled-continuous --report-monomorphic >var.vcf" << endl
        << endl
        << endl
        << "parameters:" << endl
        << endl
        << "   -h --help       Prints this help dialog." << endl
        << "   --version       Prints the release number and the git commit id." << endl
        << endl
        << "input:" << endl
        << endl
        << "   -b --bam FILE   Add FILE to the set of BAM files to be analyzed." << endl
        << "   -L --bam-list FILE" << endl
        << "                   A file containing a list of BAM files to be analyzed." << endl
        << "   -c --stdin      Read BAM input on stdin." << endl
        << "   -f --fasta-reference FILE" << endl
        << "                   Use FILE as the reference sequence for analysis." << endl
        << "                   An index file (FILE.fai) will be created if none exists." << endl
        << "                   If neither --targets nor --region are specified, FreeBayes" << endl
        << "                   will analyze every position in this reference." << endl
        << "   -t --targets FILE" << endl
        << "                   Limit analysis to targets listed in the BED-format FILE." << endl
        << "   -r --region <chrom>:<start_position>-<end_position>" << endl
        << "                   Limit analysis to the specified region, 0-base coordinates," << endl
        << "                   end_position not included (same as BED format)." << endl
        << "                   Either '-' or '..' maybe used as a separator." << endl
        << "   -s --samples FILE" << endl
        << "                   Limit analysis to samples listed (one per line) in the FILE." << endl
        << "                   By default FreeBayes will analyze all samples in its input" << endl
        << "                   BAM files." << endl
        << "   --populations FILE" << endl
        << "                   Each line of FILE should list a sample and a population which" << endl
        << "                   it is part of.  The population-based bayesian inference model" << endl
        << "                   will then be partitioned on the basis of the populations." << endl
        << "   -A --cnv-map FILE" << endl
        << "                   Read a copy number map from the BED file FILE, which has" << endl
        << "                   the format:" << endl
        << "                      reference sequence, start, end, sample name, copy number" << endl
        << "                   ... for each region in each sample which does not have the" << endl
        << "                   default copy number as set by --ploidy." << endl
        << endl
        << "output:" << endl
        << endl
        << "   -v --vcf FILE   Output VCF-format results to FILE. (default: stdout)" << endl
        << "   --gvcf" << endl
        << "                   Write gVCF output, which indicates coverage in uncalled regions." << endl
        << "   --gvcf-chunk NUM" << endl
        << "                   When writing gVCF output emit a record for every NUM bases." << endl
        << "   -@ --variant-input VCF" << endl
        << "                   Use variants reported in VCF file as input to the algorithm." << endl
        << "                   Variants in this file will included in the output even if" << endl
        << "                   there is not enough support in the data to pass input filters." << endl
        << "   -l --only-use-input-alleles" << endl
        << "                   Only provide variant calls and genotype likelihoods for sites" << endl
        << "                   and alleles which are provided in the VCF input, and provide" << endl
        << "                   output in the VCF for all input alleles, not just those which" << endl
        << "                   have support in the data." << endl
        << "   --haplotype-basis-alleles VCF" << endl
        << "                   When specified, only variant alleles provided in this input" << endl
        << "                   VCF will be used for the construction of complex or haplotype" << endl
        << "                   alleles." << endl
        << "   --report-all-haplotype-alleles" << endl
        << "                   At sites where genotypes are made over haplotype alleles," << endl
        << "                   provide information about all alleles in output, not only" << endl
        << "                   those which are called." << endl
        << "   --report-monomorphic" << endl
        << "                   Report even loci which appear to be monomorphic, and report all" << endl
        << "                   considered alleles, even those which are not in called genotypes." << endl
        << "                   Loci which do not have any potential alternates have '.' for ALT." << endl
        << "   -P --pvar N     Report sites if the probability that there is a polymorphism" << endl
        << "                   at the site is greater than N.  default: 0.0.  Note that post-" << endl
        << "                   filtering is generally recommended over the use of this parameter." << endl
        << "   --strict-vcf" << endl
        << "                   Generate strict VCF format (FORMAT/GQ will be an int)" << endl
        << endl
        << "population model:" << endl
        << endl
        << "   -T --theta N    The expected mutation rate or pairwise nucleotide diversity" << endl
        << "                   among the population under analysis.  This serves as the" << endl
        << "                   single parameter to the Ewens Sampling Formula prior model" << endl
        << "                   default: 0.001" << endl
        << "   -p --ploidy N   Sets the default ploidy for the analysis to N.  default: 2" << endl
        << "   -J --pooled-discrete" << endl
        << "                   Assume that samples result from pooled sequencing." << endl
        << "                   Model pooled samples using discrete genotypes across pools." << endl
        << "                   When using this flag, set --ploidy to the number of" << endl
        << "                   alleles in each sample or use the --cnv-map to define" << endl
        << "                   per-sample ploidy." << endl
        << "   -K --pooled-continuous" << endl
        << "                   Output all alleles which pass input filters, regardles of" << endl
        << "                   genotyping outcome or model." << endl
        << endl
        << "reference allele:" << endl
        << endl
        << "   -Z --use-reference-allele" << endl
        << "                   This flag includes the reference allele in the analysis as" << endl
        << "                   if it is another sample from the same population." << endl
        << "   --reference-quality MQ,BQ" << endl
        << "                   Assign mapping quality of MQ to the reference allele at each" << endl
        << "                   site and base quality of BQ.  default: 100,60" << endl
        << endl
        << "allele scope:" << endl
        << endl
        << "   -I --no-snps    Ignore SNP alleles." << endl
        << "   -i --no-indels  Ignore insertion and deletion alleles." << endl
        << "   -X --no-mnps    Ignore multi-nuceotide polymorphisms, MNPs." << endl
        << "   -u --no-complex Ignore complex events (composites of other classes)." << endl
        << "   -n --use-best-n-alleles N" << endl
        << "                   Evaluate only the best N SNP alleles, ranked by sum of" << endl
        << "                   supporting quality scores.  (Set to 0 to use all; default: all)" << endl
        << "   -E --max-complex-gap N" << endl
        << "      --haplotype-length N" << endl
        << "                   Allow haplotype calls with contiguous embedded matches of up" << endl
        << "                   to this length.  (default: 3)" << endl
        << "   --min-repeat-size N" << endl
        << "                   When assembling observations across repeats, require the total repeat" << endl
        << "                   length at least this many bp.  (default: 5)" << endl
        << "   --min-repeat-entropy N" << endl
        << "                   To detect interrupted repeats, build across sequence until it has" << endl
        << "                   entropy > N bits per bp.  (default: 0, off)" << endl
        << "   --no-partial-observations" << endl
        << "                   Exclude observations which do not fully span the dynamically-determined" << endl
        << "                   detection window.  (default, use all observations, dividing partial" << endl
        << "                   support across matching haplotypes when generating haplotypes.)" << endl
        << endl
        << "indel realignment:" << endl
        << endl
        << "   -O --dont-left-align-indels" << endl
        << "                   Turn off left-alignment of indels, which is enabled by default." << endl
        << endl
        << "input filters:" << endl
        << endl
        << "   -4 --use-duplicate-reads" << endl
        << "                   Include duplicate-marked alignments in the analysis." << endl
        << "                   default: exclude duplicates marked as such in alignments" << endl
        << "   -m --min-mapping-quality Q" << endl
        << "                   Exclude alignments from analysis if they have a mapping" << endl
        << "                   quality less than Q.  default: 1" << endl
        << "   -q --min-base-quality Q" << endl
        << "                   Exclude alleles from analysis if their supporting base" << endl
        << "                   quality is less than Q.  default: 0" << endl
        << "   -R --min-supporting-allele-qsum Q" << endl
        << "                   Consider any allele in which the sum of qualities of supporting" << endl
        << "                   observations is at least Q.  default: 0" << endl
        << "   -Y --min-supporting-mapping-qsum Q" << endl
        << "                   Consider any allele in which and the sum of mapping qualities of" << endl
        << "                   supporting reads is at least Q.  default: 0" << endl
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
        << "   -0 --standard-filters  Use stringent input base and mapping quality filters" << endl
        << "                   Equivalent to -m 30 -q 20 -R 0 -S 0" << endl
        << "   -F --min-alternate-fraction N" << endl
        << "                   Require at least this fraction of observations supporting" << endl
        << "                   an alternate allele within a single individual in the" << endl
        << "                   in order to evaluate the position.  default: 0.2" << endl
        << "   -C --min-alternate-count N" << endl
        << "                   Require at least this count of observations supporting" << endl
        << "                   an alternate allele within a single individual in order" << endl
        << "                   to evaluate the position.  default: 2" << endl
        << "   -3 --min-alternate-qsum N" << endl
        << "                   Require at least this sum of quality of observations supporting" << endl
        << "                   an alternate allele within a single individual in order" << endl
        << "                   to evaluate the position.  default: 0" << endl
        << "   -G --min-alternate-total N" << endl
        << "                   Require at least this count of observations supporting" << endl
        << "                   an alternate allele within the total population in order" << endl
        << "                   to use the allele in analysis.  default: 1" << endl
        << "   --min-coverage N" << endl
        << "                   Require at least this coverage to process a site. default: 0" << endl
        << "   --max-coverage N" << endl
        << "                   Do not process sites with greater than this coverage. default: no limit" << endl
        << endl
        << "population priors:" << endl
        << endl
        << "   -k --no-population-priors" << endl
        << "                   Equivalent to --pooled-discrete --hwe-priors-off and removal of" << endl
        << "                   Ewens Sampling Formula component of priors." << endl
        << endl
        << "mappability priors:" << endl
        << endl
        << "   -w --hwe-priors-off" << endl
        << "                   Disable estimation of the probability of the combination" << endl
        << "                   arising under HWE given the allele frequency as estimated" << endl
        << "                   by observation frequency." << endl
        << "   -V --binomial-obs-priors-off" << endl
        << "                   Disable incorporation of prior expectations about observations." << endl
        << "                   Uses read placement probability, strand balance probability," << endl
        << "                   and read position (5'-3') probability." << endl
        << "   -a --allele-balance-priors-off" << endl
        << "                   Disable use of aggregate probability of observation balance between alleles" << endl
        << "                   as a component of the priors." << endl
        << endl
        << "genotype likelihoods:" << endl
        << endl
        << "   --observation-bias FILE" << endl
        << "                   Read length-dependent allele observation biases from FILE." << endl
        << "                   The format is [length] [alignment efficiency relative to reference]" << endl
        << "                   where the efficiency is 1 if there is no relative observation bias." << endl
        << "   --base-quality-cap Q" << endl
        << "                   Limit estimated observation quality by capping base quality at Q." << endl
        << "   --prob-contamination F" << endl
        << "                   An estimate of contamination to use for all samples.  default: 10e-9" << endl
        << "   --legacy-gls    Use legacy (polybayes equivalent) genotype likelihood calculations" << endl
        << "   --contamination-estimates FILE" << endl
        << "                   A file containing per-sample estimates of contamination, such as" << endl
        << "                   those generated by VerifyBamID.  The format should be:" << endl
        << "                       sample p(read=R|genotype=AR) p(read=A|genotype=AA)" << endl
        << "                   Sample '*' can be used to set default contamination estimates." << endl
        << endl
        << "algorithmic features:" << endl
        << endl
        << "   --report-genotype-likelihood-max" << endl
        << "                   Report genotypes using the maximum-likelihood estimate provided" << endl
        << "                   from genotype likelihoods." << endl
        << "   -B --genotyping-max-iterations N" << endl
        << "                   Iterate no more than N times during genotyping step. default: 1000." << endl
        << "   --genotyping-max-banddepth N" << endl
        << "                   Integrate no deeper than the Nth best genotype by likelihood when" << endl
        << "                   genotyping. default: 6." << endl
        << "   -W --posterior-integration-limits N,M" << endl
        << "                   Integrate all genotype combinations in our posterior space" << endl
        << "                   which include no more than N samples with their Mth best" << endl
        << "                   data likelihood. default: 1,3." << endl
        << "   -N --exclude-unobserved-genotypes" << endl
        << "                   Skip sample genotypings for which the sample has no supporting reads." << endl
        << "   -S --genotype-variant-threshold N" << endl
        << "                   Limit posterior integration to samples where the second-best" << endl
        << "                   genotype likelihood is no more than log(N) from the highest" << endl
        << "                   genotype likelihood for the sample.  default: ~unbounded" << endl
        << "   -j --use-mapping-quality" << endl
        << "                   Use mapping quality of alleles when calculating data likelihoods." << endl
        << "   -H --harmonic-indel-quality" << endl
        << "                   Use a weighted sum of base qualities around an indel, scaled by the" << endl
        << "                   distance from the indel.  By default use a minimum BQ in flanking sequence." << endl
        << "   -D --read-dependence-factor N" << endl
        << "                   Incorporate non-independence of reads by scaling successive" << endl
        << "                   observations by this factor during data likelihood" << endl
        << "                   calculations.  default: 0.9" << endl
        << "   -= --genotype-qualities" << endl
        << "                   Calculate the marginal probability of genotypes and report as GQ in" << endl
        << "                   each sample field in the VCF output." << endl
        << endl
        << "debugging:" << endl
        << endl
        << "   -d --debug      Print debugging output." << endl
        << "   -dd             Print more verbose debugging output (requires \"make DEBUG\")" << endl
        << endl
        << endl
        << "author:   Erik Garrison <erik.garrison@bc.edu>, Marth Lab, Boston College, 2010-2014" << endl
        << "version:  " << VERSION_GIT << endl;

}


Parameters::Parameters(int argc, char** argv) {

    if (argc == 1) {
        simpleUsage(argv);
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
    samples = "";              // -s --samples
    populationsFile = "";
    cnvFile = "";
    output = "vcf";               // -v --vcf
    outputFile = "";
    gVCFout = false;
    gVCFchunk = 0;
    alleleObservationBiasFile = "";

    // operation parameters
    useDuplicateReads = false;      // -E --use-duplicate-reads
    suppressOutput = false;         // -N --suppress-output
    useBestNAlleles = 0;         // -n --use-best-n-alleles
    forceRefAllele = false;         // -Z --use-reference-allele
    useRefAllele = false;           // .....
    diploidReference = false;      // -H --diploid-reference
    allowIndels = true;            // -i --no-indels
    leftAlignIndels = true;       // -O --dont-left-align-indels
    allowMNPs = true;            // -X --no-mnps
    allowSNPs = true;          // -I --no-snps
    allowComplex = true;
    strictVCF = false;
    maxComplexGap = 3;
    //maxHaplotypeLength = 100;
    minRepeatSize = 5;
    minRepeatEntropy = 0;
    usePartialObservations = true;
    pooledDiscrete = false;                 // -J --pooled
    pooledContinuous = false;
    ewensPriors = true;
    permute = true;                // -K --permute
    useMappingQuality = false;
    useMinIndelQuality = true;
    obsBinomialPriors = true;
    hwePriors = true;
    alleleBalancePriors = true;
    excludeUnobservedGenotypes = false;
    excludePartiallyObservedGenotypes = false;
    genotypeVariantThreshold = 0;
    siteSelectionMaxIterations = 5;
    reportGenotypeLikelihoodMax = false;
    genotypingMaxIterations = 1000;
    genotypingMaxBandDepth = 7;
    minPairedAltCount = 0;
    minAltMeanMapQ = 0;
    limitGL = 0;
    reportAllHaplotypeAlleles = false;
    reportMonomorphic = false;
    boundIndels = true; // ignore indels at ends of reads
    onlyUseInputAlleles = false;
    standardGLs = false; // use experimental gls by default // XXX
    MQR = 100;                     // -M --reference-mapping-quality
    BQR = 60;                     // -B --reference-base-quality
    ploidy = 2;                  // -p --ploidy
    MQL0 = 1;                    // -m --min-mapping-quality
    BQL0 = 0;                    // -q --min-base-quality
    minSupportingAlleleQualitySum = 0;
    minSupportingMappingQualitySum = 0;
    BQL2 = 10;                    // -Q --mismatch-base-quality-threshold
    RMU = 10000000;                     // -U --read-mismatch-limit
    readMaxMismatchFraction = 1.0;    //  -z --read-max-mismatch-fraction
    readSnpLimit = 10000000;       // -$ --read-snp-limit
    readIndelLimit = 10000000;     // -e --read-indel-limit
    IDW = -1;                     // -x --indel-exclusion-window
    TH = 10e-3;              // -T --theta
    PVL = 0.0;             // -P --pvar
    RDF = 0.9;             // -D --read-dependence-factor
    diffusionPriorScalar = 1.0;     // -V --diffusion-prior-scalar
    WB = 1;                      // -W --posterior-integration-limits
    TB = 3;
    posteriorIntegrationDepth = 0;
    calculateMarginals = false;
    minAltFraction = 0.2;  // require 20% of reads from sample to be supporting the same alternate to consider
    minAltCount = 2; // require 2 reads in same sample call
    minAltTotal = 1;
    minAltQSum = 0;
    baseQualityCap = 0;
    probContamination = 10e-9;
    //minAltQSumTotal = 0;
    minCoverage = 0;
    maxCoverage = 0;
    debuglevel = 0;
    debug = false;
    debug2 = false;

    showReferenceRepeats = false;

    int c; // counter for getopt

    static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"version", no_argument, 0, '#'},
            {"bam", required_argument, 0, 'b'},
            {"bam-list", required_argument, 0, 'L'},
            {"stdin", no_argument, 0, 'c'},
            {"fasta-reference", required_argument, 0, 'f'},
            {"targets", required_argument, 0, 't'},
            {"region", required_argument, 0, 'r'},
            {"samples", required_argument, 0, 's'},
            {"populations", required_argument, 0, '2'},
            {"cnv-map", required_argument, 0, 'A'},
            {"vcf", required_argument, 0, 'v'},
            {"gvcf", no_argument, 0, '8'},
            {"gvcf-chunk", required_argument, 0, '&'},
            {"use-duplicate-reads", no_argument, 0, '4'},
            {"no-partial-observations", no_argument, 0, '['},
            {"use-best-n-alleles", required_argument, 0, 'n'},
            {"use-reference-allele", no_argument, 0, 'Z'},
            {"harmonic-indel-quality", no_argument, 0, 'H'},
            {"standard-filters", no_argument, 0, '0'},
            {"reference-quality", required_argument, 0, '1'},
            {"ploidy", required_argument, 0, 'p'},
            {"pooled-discrete", no_argument, 0, 'J'},
            {"pooled-continuous", no_argument, 0, 'K'},
            {"no-population-priors", no_argument, 0, 'k'},
            {"use-mapping-quality", no_argument, 0, 'j'},
            {"min-mapping-quality", required_argument, 0, 'm'},
            {"min-base-quality", required_argument, 0, 'q'},
            {"min-supporting-allele-qsum", required_argument, 0, 'R'},
            {"min-supporting-mapping-qsum", required_argument, 0, 'Y'},
            {"mismatch-base-quality-threshold", required_argument, 0, 'Q'},
            {"read-mismatch-limit", required_argument, 0, 'U'},
            {"read-max-mismatch-fraction", required_argument, 0, 'z'},
            {"read-snp-limit", required_argument, 0, '$'},
            {"read-indel-limit", required_argument, 0, 'e'},
            {"no-indels", no_argument, 0, 'i'},
            {"dont-left-align-indels", no_argument, 0, 'O'},
            {"no-mnps", no_argument, 0, 'X'},
            {"no-complex", no_argument, 0, 'u'},
            {"max-complex-gap", required_argument, 0, 'E'},
            {"haplotype-length", required_argument, 0, 'E'},
            {"min-repeat-size", required_argument, 0, 'E'},
            {"min-repeat-entropy", required_argument, 0, 'E'},
            {"no-snps", no_argument, 0, 'I'},
            {"indel-exclusion-window", required_argument, 0, 'x'},
            {"theta", required_argument, 0, 'T'},
            {"pvar", required_argument, 0, 'P'},
            {"strict-vcf", no_argument, 0, '/'},
            {"read-dependence-factor", required_argument, 0, 'D'},
            {"binomial-obs-priors-off", no_argument, 0, 'V'},
            {"allele-balance-priors-off", no_argument, 0, 'a'},
            {"hwe-priors-off", no_argument, 0, 'w'},
            {"posterior-integration-limits", required_argument, 0, 'W'},
            {"min-alternate-fraction", required_argument, 0, 'F'},
            {"min-alternate-count", required_argument, 0, 'C'},
            //{"min-paired-alternate-count", required_argument, 0, 'Y'},
            {"observation-bias", required_argument, 0, '%'},
            {"min-alternate-total", required_argument, 0, 'G'},
            //{"min-alternate-mean-mapq", required_argument, 0, 'k'},
            {"min-alternate-qsum", required_argument, 0, '3'},
            {"min-coverage", required_argument, 0, '!'},
            {"max-coverage", required_argument, 0, '+'},
            {"genotype-qualities", no_argument, 0, '='},
            {"variant-input", required_argument, 0, '@'},
            {"only-use-input-alleles", no_argument, 0, 'l'},
            //{"show-reference-repeats", no_argument, 0, '_'},
            {"exclude-unobserved-genotypes", no_argument, 0, 'N'},
            {"genotype-variant-threshold", required_argument, 0, 'S'},
            {"site-selection-max-iterations", required_argument, 0, 'M'},
            {"genotyping-max-iterations", required_argument, 0, 'B'},
            {"genotyping-max-banddepth", required_argument, 0, '7'},
            {"haplotype-basis-alleles", required_argument, 0, '9'},
            {"report-genotype-likelihood-max", no_argument, 0, '5'},
            {"report-all-haplotype-alleles", no_argument, 0, '6'},
            {"base-quality-cap", required_argument, 0, '('},
            {"legacy-gls", no_argument, 0, ')'},
            {"prob-contamination", required_argument, 0, '_'},
            {"contamination-estimates", required_argument, 0, ','},
            {"report-monomorphic", no_argument, 0, '6'},
            {"debug", no_argument, 0, 'd'},
            {0, 0, 0, 0}

        };

    while (true) {

        int option_index = 0;
        c = getopt_long(argc, argv, "hcO4ZKjH[0diN5a)Ik=wl6#uVXJY:b:G:M:x:@:A:f:t:r:s:v:n:B:p:m:q:R:Q:U:$:e:T:P:D:^:S:W:F:C:&:L:8z:1:3:E:7:2:9:%:_:,:(:!:+:",
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
            regions.push_back(optarg);
            break;

            // -s --samples
        case 's':
            samples = optarg;
            break;

            // --populations
        case '2':
            populationsFile = optarg;
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

            // -O --dont-left-align-indels
        case 'O':
            leftAlignIndels = false;
            break;

            // --bam-list
        case 'L':
            addLinesFromFile(bams, string(optarg));
            break;

            // -8 --gvcf
        case '8':
            gVCFout = true;
            break;

        case '&':
            gVCFchunk = atoi(optarg);
            break;

            // -4 --use-duplicate-reads
        case '4':
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

            // -+ --max-coverage
        case '+':
            if (!convert(optarg, maxCoverage)) {
                cerr << "could not parse max-coverage" << endl;
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

            // -Z --use-reference-allele
        case 'Z':
            forceRefAllele = true;
            useRefAllele = true;
            break;

            // -H --harmonic-indel-quality
        case 'H':
            useMinIndelQuality = false;
            break;

            // -0 --standard-filters
        case '0':
            MQL0 = 30;
            BQL0 = 20;
            break;

            // -M --expectation-maximization
        case 'M':
            if (!convert(optarg, siteSelectionMaxIterations)) {
                cerr << "could not parse site-selection-max-iterations" << endl;
                exit(1);
            }
            break;

        case '/':
            strictVCF = true;
            break;

        case 'u':
            allowComplex = false;
            break;

        case 'E':
        {
            string arg(argv[optind - 2]);
            if (arg == "--min-repeat-size") {
                if (!convert(optarg, minRepeatSize)) {
                    cerr << "could not parse " << arg << endl;
                    exit(1);
                }
            } else if (arg == "--min-repeat-entropy") {
                if (!convert(optarg, minRepeatEntropy)) {
                    cerr << "could not parse " << arg << endl;
                    exit(1);
                }
            } else {
                if (!convert(optarg, maxComplexGap)) {
                    cerr << "could not parse maxComplexGap" << endl;
                    exit(1);
                }
            }
            break;
        }

            // -B --genotyping-max-iterations
        case 'B':
            if (!convert(optarg, genotypingMaxIterations)) {
                cerr << "could not parse genotyping-max-iterations" << endl;
                exit(1);
            }
            break;

            // -7 --genotyping-max-banddepth
        case '7':
            if (!convert(optarg, genotypingMaxBandDepth)) {
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

        case 'J':
            pooledDiscrete = true;
            hwePriors = false;  // disable hwe sampling prob when using discrete pooling
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

            // -R --min-supporting-allele-qsum
        case 'R':
            if (!convert(optarg, minSupportingAlleleQualitySum)) {
                cerr << "could not parse min-supporting-allele-qsum" << endl;
                exit(1);
            }
            break;

            // -Y --min-supporting-mapping-quality
        case 'Y':
            if (!convert(optarg, minSupportingMappingQualitySum)) {
                cerr << "could not parse min-supporting-mapping-qsum" << endl;
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

        case '5':
            reportGenotypeLikelihoodMax = true;
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

            // -x --indel-exclusion-window
        case 'x':
            if (!convert(optarg, IDW)) {
                cerr << "could not parse indel-exclusion-window" << endl;
                exit(1);
            }
            break;

            // -i --indels
        case 'i':
            allowIndels = false;
            break;

            // -X --mnps
        case 'X':
            allowMNPs = false;
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

            // -% --observation-bias
        case '%':
            alleleObservationBiasFile = optarg;
            break;

            // observation priors
        case 'V':
            obsBinomialPriors = false;
            break;

            // allele balance
        case 'a':
            alleleBalancePriors = false;
            break;

            // hwe expectations
        case 'w':
            hwePriors = false;
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

            // -k --no-population-priors
        case 'k':
            pooledDiscrete = true;
            ewensPriors = false;
            hwePriors = false;
            break;

        case 'K':
            pooledContinuous = true;
            break;

        case '=':
            calculateMarginals = true;
            break;

        case '@':
            variantPriorsFile = optarg;
            break;

        case '9':
            haplotypeVariantFile = optarg;
            break;

        case 'l':
            onlyUseInputAlleles = true;
            break;

        case '6':
        {
            string arg(argv[optind - 1]);
            if (arg == "--report-monomorphic") {
                reportMonomorphic = true;
            }
            reportAllHaplotypeAlleles = true;
        }
            break;

        case '[':
            usePartialObservations = false;
            break;

        case '_':
            if (!convert(optarg, probContamination)) {
                cerr << "could not parse prob-contamination" << endl;
                exit(1);
            }
            break;

        case ',':
            contaminationEstimateFile = optarg;
            break;

        case ')':
            standardGLs = true;
            break;

        case '(':
            if (!convert(optarg, baseQualityCap)) {
                cerr << "could not parse base-quality-cap" << endl;
                exit(1);
            }
            break;

            // -d --debug
        case 'd':
            ++debuglevel;
            break;

    case '#':
        
        // --version
            cout << "version:  " << VERSION_GIT << endl;
        exit(0);
        break;

        case 'h':
            usage(argv);
            exit(0);
            break;
 
            // either catch "long options" or
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

    // check that there aren't duplicates in the bams list
    for( int i=1; i<bams.size(); ++i ){
        for( int j=0; j<i; ++j ){
            if( bams[i] == bams[j] ){
                cerr << "Error: Duplicate bam file '" << bams[i] << "'" << endl;
                exit(1);
            }
        }
    }

}

#include "AlleleParser.h"
#include "multichoose.h" // includes generic functions, so it must be included here
                         // otherwise we will get a linker error
                         // see: http://stackoverflow.com/questions/36039/templates-spread-across-multiple-files
                         // http://www.cplusplus.com/doc/tutorial/templates/ "Templates and Multi-file projects"
#include "multipermute.h"

// local helper debugging macros to improve code readability
#define DEBUG(msg) \
    if (parameters.debug) { cerr << msg << endl; }

// lower-priority messages
#ifdef VERBOSE_DEBUG
#define DEBUG2(msg) \
    if (parameters.debug2) { cerr << msg << endl; }
#else
#define DEBUG2(msg)
#endif

// must-see error messages
#define ERROR(msg) \
    cerr << msg << endl;

using namespace std;


// open BAM input file
void AlleleParser::openBams(void) {

    // report differently if we have one or many bam files
    if (parameters.bams.size() == 1) {
        DEBUG("Opening BAM fomat alignment input file: " << parameters.bams.front() << " ...");
    } else if (parameters.bams.size() > 1) {
        DEBUG("Opening " << parameters.bams.size() << " BAM fomat alignment input files");
        for (vector<string>::const_iterator b = parameters.bams.begin(); 
                b != parameters.bams.end(); ++b) {
            DEBUG2(*b);
        }
    }
    
    // set no index caching if we are only making one jump
    if (targets.size() == 1) {
        bamMultiReader.SetIndexCacheMode(BamIndex::NoIndexCaching);
    }

    if (parameters.useStdin) {
        if (!bamMultiReader.Open(parameters.bams)) {
            ERROR("Could not read BAM data from stdin");
            exit(1);
        }
    } else {
        if (!bamMultiReader.Open(parameters.bams)) {
            ERROR("Could not open input BAM files");
            exit(1);
        } else {
            if (!bamMultiReader.LocateIndexes(BamIndex::BAMTOOLS)) {
                ERROR("Opened BAM reader without index file, jumping is disabled.");
                if (!targets.empty()) {
                    ERROR("Targets specified but no BAM index file provided.");
                    ERROR("FreeBayes cannot jump through targets in BAM files without BAM index files, exiting.");
                    ERROR("Please generate a BAM index file eithe, e.g.:");
                    ERROR("    \% bamtools index -in <bam_file>");
                    ERROR("    \% samtools index <bam_file>");
                    exit(1);
                }
            }
        }
    }


    // retrieve header information
    bamHeader = bamMultiReader.GetHeaderText();
    bamHeaderLines = split(bamHeader, '\n');

    DEBUG(" done");

}

void AlleleParser::openTraceFile(void) {
    if (parameters.trace) {
        traceFile.open(parameters.traceFile.c_str(), ios::out);
        DEBUG("Opening trace file: " << parameters.traceFile << " ...");
        if (!traceFile) {
            ERROR(" unable to open trace file: " << parameters.traceFile );
            exit(1);
        }
    }
}

void AlleleParser::openFailedFile(void) {
    if (!parameters.failedFile.empty()) {
        failedFile.open(parameters.failedFile.c_str(), ios::out);
        DEBUG("Opening failed alleles file: " << parameters.failedFile << " ...");
        if (!failedFile) {
            ERROR(" unable to open failed alleles file: " << parameters.failedFile );
            exit(1);
        }
    }
}

void AlleleParser::openOutputFile(void) {
    if (parameters.outputFile != "") {
        outputFile.open(parameters.outputFile.c_str(), ios::out);
        DEBUG("Opening output file: " << parameters.outputFile << " ...");
        if (!outputFile) {
            ERROR(" unable to open output file: " << parameters.outputFile);
            exit(1);
        }
        output = &outputFile;
    } else {
        output = &cout;
    }
}

void AlleleParser::getSequencingTechnologies(void) {

    map<string, bool> technologies;

    for (vector<string>::const_iterator it = bamHeaderLines.begin(); it != bamHeaderLines.end(); ++it) {

        // get next line from header, skip if empty
        string headerLine = *it;
        if ( headerLine.empty() ) { continue; }

        // lines of the header look like:
        // "@RG     ID:-    SM:NA11832      CN:BCM  PL:454"
        //                     ^^^^^^^\ is our sample name
        if ( headerLine.find("@RG") == 0 ) {
            vector<string> readGroupParts = split(headerLine, "\t ");
            string tech;
            string readGroupID;
            for (vector<string>::const_iterator r = readGroupParts.begin(); r != readGroupParts.end(); ++r) {
                vector<string> nameParts = split(*r, ":");
                if (nameParts.at(0) == "PL") {
                   tech = nameParts.at(1);
                } else if (nameParts.at(0) == "ID") {
                   readGroupID = nameParts.at(1);
                }
            }
            if (tech.empty()) {
                cerr << " could not find PL: in @RG tag " << endl << headerLine << endl;
                continue;
            }
            if (readGroupID.empty()) {
                cerr << " could not find ID: in @RG tag " << endl << headerLine << endl;
                continue;
            }
            //string name = nameParts.back();
            //mergedHeader.append(1, '\n');
            //cerr << "found read group id " << readGroupID << " containing sample " << name << endl;
            readGroupToTechnology[readGroupID] = tech;
            technologies[tech] = true;
        }
    }

    for (map<string, bool>::iterator st = technologies.begin(); st != technologies.end(); ++st) {
        sequencingTechnologies.push_back(st->first);
    }

}

void AlleleParser::getPopulations(void) {

    map<string, string> allSamplePopulation;

    if (!parameters.populationsFile.empty()) {
        ifstream populationsFile(parameters.populationsFile.c_str(), ios::in);
        if (!populationsFile) {
            cerr << "unable to open population file: " << parameters.populationsFile << endl;
            exit(1);
        }
        string line;
        while (getline(populationsFile, line)) {
            DEBUG2("found sample-population mapping: " << line);
            vector<string> popsample = split(line, "\t ");
            if (popsample.size() == 2) {
                string& sample = popsample.front();
                string& population = popsample.back();
                DEBUG2("sample: " << sample << " population: " << population);
                allSamplePopulation[sample] = population;
            } else {
                cerr << "malformed population/sample pair, " << line << endl;
                exit(1);
            }
        }
    }

    // XXX
    // TODO now, assign a default population to all the rest of the samples...
    // XXX
    for (vector<string>::iterator s = sampleList.begin(); s != sampleList.end(); ++s) {
        if (!allSamplePopulation.count(*s)) {
            samplePopulation[*s] = "DEFAULT";
        } else {
            samplePopulation[*s] = allSamplePopulation[*s];
        }
    }

    // now, only keep the samples we are using for processing


    for (map<string, string>::iterator s = samplePopulation.begin(); s != samplePopulation.end(); ++s) {
        populationSamples[s->second].push_back(s->first);
    }

}

// read sample list file or get sample names from bam file header
void AlleleParser::getSampleNames(void) {

    // If a sample file is given, use it.  But otherwise process the bam file
    // header to get the sample names.
    //

    if (!parameters.samples.empty()) {
        ifstream sampleFile(parameters.samples.c_str(), ios::in);
        if (! sampleFile) {
            cerr << "unable to open file: " << parameters.samples << endl;
            exit(1);
        }
        string line;
        while (getline(sampleFile, line)) {
            DEBUG2("found sample " << line);
            sampleList.push_back(line);
        }
    }

    for (vector<string>::const_iterator it = bamHeaderLines.begin(); it != bamHeaderLines.end(); ++it) {

        // get next line from header, skip if empty
        string headerLine = *it;
        if ( headerLine.empty() ) { continue; }

        // lines of the header look like:
        // "@RG     ID:-    SM:NA11832      CN:BCM  PL:454"
        //                     ^^^^^^^\ is our sample name
        if ( headerLine.find("@RG") == 0 ) {
            vector<string> readGroupParts = split(headerLine, "\t ");
            string name = "";
            string readGroupID = "";
            for (vector<string>::const_iterator r = readGroupParts.begin(); r != readGroupParts.end(); ++r) {
                vector<string> nameParts = split(*r, ":");
                if (nameParts.at(0) == "SM") {
                   name = nameParts.at(1);
                } else if (nameParts.at(0) == "ID") {
                   readGroupID = nameParts.at(1);
                }
            }
            if (name == "") {
                ERROR(" could not find SM: in @RG tag " << endl << headerLine);
                exit(1);
            }
            if (readGroupID == "") {
                ERROR(" could not find ID: in @RG tag " << endl << headerLine);
                exit(1);
            }
            //string name = nameParts.back();
            //mergedHeader.append(1, '\n');
            DEBUG2("found read group id " << readGroupID << " containing sample " << name);
            sampleListFromBam.push_back(name);

            map<string, string>::iterator s = readGroupToSampleNames.find(readGroupID);
            if (s != readGroupToSampleNames.end()) {
                if (s->second != name) {
                    ERROR("ERROR: multiple samples (SM) map to the same read group (RG)" << endl
                       << endl
                       << "samples " << name << " and " << s->second << " map to " << readGroupID << endl
                       << endl
                       << "As freebayes operates on a virtually merged stream of its input files," << endl 
                       << "it will not be possible to determine what sample an alignment belongs to" << endl
                       << "at runtime." << endl
                       << endl
                       << "To resolve the issue, ensure that RG ids are unique to one sample" << endl
                       << "across all the input files to freebayes." << endl
                       << endl
                       << "See bamaddrg (https://github.com/ekg/bamaddrg) for a method which can" << endl
                       << "add RG tags to alignments." << endl);
                    exit(1);
                }
                // if it's the same sample name and RG combo, no worries
            }
            readGroupToSampleNames[readGroupID] = name;
        }
    }
    //cout << sampleListFromBam.size() << endl;
     // no samples file given, read from BAM file header for sample names
    if (sampleList.empty()) {
        DEBUG("no sample list file given, reading sample names from bam file");
        for (vector<string>::const_iterator s = sampleListFromBam.begin(); s != sampleListFromBam.end(); ++s) {
            DEBUG2("found sample " << *s);
            if (!stringInVector(*s, sampleList)) {
                sampleList.push_back(*s);
            }
        }
        DEBUG("found " << sampleList.size() << " samples in BAM file");
    } else {
        // verify that the samples in the sample list are present in the bam,
        // and raise an error and exit if not
        for (vector<string>::const_iterator s = sampleList.begin(); s != sampleList.end(); ++s) {
            bool inBam = false;
            bool inReadGroup = false;
            //cout << "checking sample from sample file " << *s << endl;
            for (vector<string>::const_iterator b = sampleListFromBam.begin(); b != sampleListFromBam.end(); ++b) {
                //cout << *s << " against " << *b << endl;
                if (*s == *b) { inBam = true; break; }
            }
            for (map<string, string>::const_iterator p = readGroupToSampleNames.begin(); p != readGroupToSampleNames.end(); ++p) {
                if (*s == p->second) { inReadGroup = true; break; }
            }
            if (!inBam) {
                ERROR("sample " << *s << " listed in sample file "
                    << parameters.samples.c_str() << " is not listed in the header of BAM file(s) "
                    << parameters.bam);
                exit(1);
            }
            if (!inReadGroup) {
                ERROR("sample " << *s << " listed in sample file "
                    << parameters.samples.c_str() << " is not associated with any read group in the header of BAM file(s) "
                    << parameters.bam);
                exit(1);
            }
        }
    }

    if (sampleList.empty()) {
        /*
        ERROR(string(80, '-') << endl
             //--------------------------------------------------------------------------------
           << "Warning: No sample file given, and no @RG tags found in BAM header." << endl
           << "All alignments from all input files will be assumed to come from the same" << endl
           << "individual.  To group alignments by sample, you must add read groups and sample" << endl 
           << "names to your alignments.  You can do this using ./scripts/sam_add_rg.pl in the" << endl
           << "freebayes source tree, or by specifying read groups and sample names when you" << endl
           << "prepare your sequencing data for alignment." << endl
           << string(80, '-'));
           */
        sampleList.push_back("unknown");
        readGroupToSampleNames["unknown"] = "unknown";
        oneSampleAnalysis = true;
    }

}

string AlleleParser::vcfHeader() {

    stringstream headerss;
    headerss << "##fileformat=VCFv4.1" << endl
        << "##fileDate=" << dateStr() << endl
        << "##source=freeBayes version " << FREEBAYES_VERSION << endl
        << "##reference=" << reference.filename << endl
        << "##phasing=none" << endl
        << "##commandline=\"" << parameters.commandline << "\"" << endl
        << "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of samples with data\">" << endl
        << "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total read depth at the locus\">" << endl

        // allele frequency metrics
        << "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Total number of alternate alleles in called genotypes\">" << endl
        << "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">" << endl
        << "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Estimated allele frequency in the range (0,1]\">" << endl
        //<< "##INFO=<ID=HETAR,Number=1,Type=Integer,Description=\"Number of individuals heterozygous alternate / reference\">" << endl
        //<< "##INFO=<ID=HOMA,Number=1,Type=Integer,Description=\"Number of individuals homozygous for the alternate\">" << endl
        //<< "##INFO=<ID=HOMR,Number=1,Type=Integer,Description=\"Number of individuals homozygous for the reference\">" << endl

        // binomial balance metrics
        << "##INFO=<ID=RO,Number=1,Type=Integer,Description=\"Reference allele observations\">" << endl
        << "##INFO=<ID=AO,Number=A,Type=Integer,Description=\"Alternate allele observations\">" << endl
        //<< "##INFO=<ID=SRF,Number=1,Type=Integer,Description=\"Number of reference observations on the forward strand\">" << endl
        //<< "##INFO=<ID=SRR,Number=1,Type=Integer,Description=\"Number of reference observations on the reverse strand\">" << endl
        //<< "##INFO=<ID=SAF,Number=1,Type=Integer,Description=\"Number of alternate observations on the forward strand\">" << endl
        //<< "##INFO=<ID=SAR,Number=1,Type=Integer,Description=\"Number of alternate observations on the reverse strand\">" << endl
        //<< "##INFO=<ID=SRB,Number=1,Type=Float,Description=\"Strand bias for the reference allele: SRF / ( SRF + SRR )\">" << endl
        //<< "##INFO=<ID=SAB,Number=1,Type=Float,Description=\"Strand bias for the alternate allele: SAF / ( SAF + SAR )\">" << endl
        << "##INFO=<ID=SRP,Number=1,Type=Float,Description=\"Strand balance probability for the reference allele: Phred-scaled upper-bounds estimate of the probability of observing the deviation between SRF and SRR given E(SRF/SRR) ~ 0.5, derived using Hoeffding's inequality\">" << endl
        << "##INFO=<ID=SAP,Number=A,Type=Float,Description=\"Strand balance probability for the alternate allele: Phred-scaled upper-bounds estimate of the probability of observing the deviation between SAF and SAR given E(SAF/SAR) ~ 0.5, derived using Hoeffding's inequality\">" << endl
        //<< "##INFO=<ID=ABR,Number=1,Type=Integer,Description=\"Reference allele balance count: the number of sequence reads from apparent heterozygotes supporting the reference allele\">" << endl
        //<< "##INFO=<ID=ABA,Number=1,Type=Integer,Description=\"Alternate allele balance count: the number of sequence reads from apparent heterozygotes supporting the alternate allele\">" << endl
        << "##INFO=<ID=AB,Number=A,Type=Float,Description=\"Allele balance at heterozygous sites: a number between 0 and 1 representing the ratio of reads showing the reference allele to all reads, considering only reads from individuals called as heterozygous\">" << endl
        << "##INFO=<ID=ABP,Number=A,Type=Float,Description=\"Allele balance probability at heterozygous sites: Phred-scaled upper-bounds estimate of the probability of observing the deviation between ABR and ABA given E(ABR/ABA) ~ 0.5, derived using Hoeffding's inequality\">" << endl
        << "##INFO=<ID=RUN,Number=A,Type=Integer,Description=\"Run length: the number of consecutive repeats of the alternate allele in the reference genome\">" << endl
        //<< "##INFO=<ID=RL,Number=1,Type=Integer,Description=\"Reads Placed Left: number of reads supporting the alternate balanced to the left (5') of the alternate allele\">" << endl
        //<< "##INFO=<ID=RR,Number=1,Type=Integer,Description=\"Reads Placed Right: number of reads supporting the alternate balanced to the right (3') of the alternate allele\">" << endl
        << "##INFO=<ID=RPP,Number=A,Type=Float,Description=\"Read Placement Probability: Phred-scaled upper-bounds estimate of the probability of observing the deviation between RPL and RPR given E(RPL/RPR) ~ 0.5, derived using Hoeffding's inequality\">" << endl
        << "##INFO=<ID=RPPR,Number=1,Type=Float,Description=\"Read Placement Probability for reference observations: Phred-scaled upper-bounds estimate of the probability of observing the deviation between RPL and RPR given E(RPL/RPR) ~ 0.5, derived using Hoeffding's inequality\">" << endl
        //<< "##INFO=<ID=EL,Number=1,Type=Integer,Description=\"Allele End Left: number of observations of the alternate where the alternate occurs in the left end of the read\">" << endl
        //<< "##INFO=<ID=ER,Number=1,Type=Integer,Description=\"Allele End Right: number of observations of the alternate where the alternate occurs in the right end of the read\">" << endl
        << "##INFO=<ID=EPP,Number=A,Type=Float,Description=\"End Placement Probability: Phred-scaled upper-bounds estimate of the probability of observing the deviation between EL and ER given E(EL/ER) ~ 0.5, derived using Hoeffding's inequality\">" << endl
        << "##INFO=<ID=EPPR,Number=1,Type=Float,Description=\"End Placement Probability for reference observations: Phred-scaled upper-bounds estimate of the probability of observing the deviation between EL and ER given E(EL/ER) ~ 0.5, derived using Hoeffding's inequality\">" << endl
        //<< "##INFO=<ID=BL,Number=1,Type=Integer,Description=\"Base Pairs Left: number of base pairs in reads supporting the alternate to the left (5') of the alternate allele\">" << endl
        //<< "##INFO=<ID=BR,Number=1,Type=Integer,Description=\"Base Pairs Right: number of base pairs in reads supporting the alternate to the right (3') of the alternate allele\">" << endl
        //<< "##INFO=<ID=LRB,Number=1,Type=Float,Description=\"((max(BR, BL) / (BR + BL)) - 0.5) * 2 : The proportion of base pairs in reads on one side of the alternate allele relative to total bases, scaled from [0.5,1] to [0,1]\">" << endl
        //<< "##INFO=<ID=LRBP,Number=1,Type=Float,Description=\"Left-Right Balance Probability: Phred-scaled upper-bounds estimate of the probability of observing the deviation between BL and BR given E(BR/BL) ~ 0.5, derived using Hoeffding's inequality\">" << endl
        << "##INFO=<ID=DPRA,Number=A,Type=Float,Description=\"Alternate allele depth ratio.  Ratio between depth in samples with each called alternate allele and those without.\">" << endl

        // error rates
        << "##INFO=<ID=XRM,Number=1,Type=Float,Description=\"Reference allele read mismatch rate: The rate of SNPs + MNPs + INDELs in reads supporting the reference allele.\">" << endl
        << "##INFO=<ID=XRS,Number=1,Type=Float,Description=\"Reference allele read SNP rate: The rate of per-base mismatches (SNPs + MNPs) in reads supporting the reference allele.\">" << endl
        << "##INFO=<ID=XRI,Number=1,Type=Float,Description=\"Reference allele read INDEL rate: The rate of INDELs (gaps) in reads supporting the reference allele.\">" << endl
        << "##INFO=<ID=XAM,Number=A,Type=Float,Description=\"Alternate allele read mismatch rate: The rate of SNPs + MNPs + INDELs in reads supporting the alternate allele, excluding the called variant.\">" << endl
        << "##INFO=<ID=XAS,Number=A,Type=Float,Description=\"Alternate allele read SNP rate: The rate of per-base mismatches (SNPs + MNPs) in reads supporting the alternate allele, excluding the called variant.\">" << endl
        << "##INFO=<ID=XAI,Number=A,Type=Float,Description=\"Alternate allele read INDEL rate: The rate of INDELs (gaps) in reads supporting the alternate allele, excluding the called variant.\">" << endl
        // error rate ratios
        //<< "##INFO=<ID=ARM,Number=A,Type=Float,Description=\"Alternate allele / reference allele read mismatch ratio: The rate of SNPs + MNPs + INDELs in reads supporting the alternate allele versus reads supporting the reference allele, excluding the called variant.\">" << endl
        //<< "##INFO=<ID=ARS,Number=A,Type=Float,Description=\"Alternate allele / reference allele read SNP ratio: The rate of per-base mismatches (SNPs + MNPs) in reads supporting the alternate allele versus reads supporting the reference allele, excluding the called variant.\">" << endl
        //<< "##INFO=<ID=ARI,Number=A,Type=Float,Description=\"Alternate allele / reference allele read INDEL ratio: The ratio in rate rate of INDELs (gaps) in reads supporting the alternate allele versus reads supporting the reference allele, excluding the called variant.\">" << endl

        // supplementary information about the site
        << "##INFO=<ID=ODDS,Number=1,Type=Float,Description=\"The log odds ratio of the best genotype combination to the second-best.\">" << endl
        << "##INFO=<ID=BVAR,Number=0,Type=Flag,Description=\"The best genotype combination in the posterior is variant (non homozygous).\">" << endl
        //<< "##INFO=<ID=TS,Number=0,Type=Flag,Description=\"site has transition SNP\">" << endl
        //<< "##INFO=<ID=TV,Number=0,Type=Flag,Description=\"site has transversion SNP\">" << endl
        << "##INFO=<ID=CpG,Number=0,Type=Flag,Description=\"CpG site (either CpG, TpG or CpA)\">" << endl
        << "##INFO=<ID=TYPE,Number=A,Type=String,Description=\"The type of allele, either snp, mnp, ins, del, or complex.\">" << endl
        << "##INFO=<ID=CIGAR,Number=A,Type=String,Description=\"The extended CIGAR representation of each alternate allele, with the exception that '=' is replaced by 'M' to ease VCF parsing.  Note that INDEL alleles do not have the first matched base (which is provided by default, per the spec) referred to by the CIGAR.\">" << endl
        //<< "##INFO=<ID=SNP,Number=0,Type=Flag,Description=\"SNP allele at site\">" << endl
        //<< "##INFO=<ID=MNP,Number=0,Type=Flag,Description=\"MNP allele at site\">" << endl
        //<< "##INFO=<ID=INS,Number=0,Type=Flag,Description=\"insertion allele at site\">" << endl
        //<< "##INFO=<ID=DEL,Number=0,Type=Flag,Description=\"deletion allele at site\">" << endl
        //<< "##INFO=<ID=COMPLEX,Number=0,Type=Flag,Description=\"complex allele (insertion/deletion/substitution composite) at site\">" << endl
        << "##INFO=<ID=NUMALT,Number=1,Type=Integer,Description=\"Number of unique non-reference alleles in called genotypes at this position.\">" << endl
        << "##INFO=<ID=MEANALT,Number=A,Type=Float,Description=\"Mean number of unique non-reference allele observations per sample with the corresponding alternate alleles.\">" << endl
        << "##INFO=<ID=HWE,Number=1,Type=Float,Description=\"Phred-scaled discrete HWE prior probability of the genotyping across all samples.\">" << endl
        << "##INFO=<ID=LEN,Number=A,Type=Integer,Description=\"allele length\">" << endl
        << "##INFO=<ID=MQM,Number=A,Type=Float,Description=\"Mean mapping quality of observed alternate alleles\">" << endl
        << "##INFO=<ID=MQMR,Number=1,Type=Float,Description=\"Mean mapping quality of observed reference alleles\">" << endl
        << "##INFO=<ID=PAIRED,Number=A,Type=Float,Description=\"Proportion of observed alternate alleles which are supported by properly paired read fragments\">" << endl
        << "##INFO=<ID=PAIREDR,Number=1,Type=Float,Description=\"Proportion of observed reference alleles which are supported by properly paired read fragments\">" << endl;

    // sequencing technology tags, which vary according to input data
    for (vector<string>::iterator st = sequencingTechnologies.begin(); st != sequencingTechnologies.end(); ++st) {
        string& tech = *st;
        headerss << "##INFO=<ID=technology." << tech << ",Number=A,Type=Float,Description=\"Fraction of observations supporting the alternate observed in reads from " << tech << "\">" << endl;
    }

    if (parameters.showReferenceRepeats) {
        headerss << "##INFO=<ID=REPEAT,Number=1,Type=String,Description=\"Description of the local repeat structures flanking the current position\">" << endl;
    }

        // format fields for genotypes
    headerss << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << endl
        << "##FORMAT=<ID=GQ,Number=1,Type=Float,Description=\"Genotype Quality, the Phred-scaled marginal (or unconditional) probability of the called genotype\">" << endl
        // this can be regenerated with RA, AA, QR, QA
        << "##FORMAT=<ID=GL,Number=G,Type=Float,Description=\"Genotype Likelihood, log10-scaled likelihoods of the data given the called genotype for each possible genotype generated from the reference and alternate alleles given the sample ploidy\">" << endl
        << "##FORMAT=<ID=GLE,Number=1,Type=String,Description=\"Genotype Likelihood Explicit, same as GL, but with tags to indicate the specific genotype.  For instance, 0^-75.22|1^-223.42|0/0^-323.03|1/0^-99.29|1/1^-802.53 represents both haploid and diploid genotype likilehoods in a biallelic context\">" << endl
        << "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">" << endl
        << "##FORMAT=<ID=RO,Number=1,Type=Integer,Description=\"Reference allele observation count\">" << endl
        << "##FORMAT=<ID=QR,Number=1,Type=Integer,Description=\"Sum of quality of the reference observations\">" << endl
        << "##FORMAT=<ID=AO,Number=A,Type=Integer,Description=\"Alternate allele observation count\">" << endl
        << "##FORMAT=<ID=QA,Number=A,Type=Integer,Description=\"Sum of quality of the alternate observations\">" << endl
        // TODO (?)
        //<< "##FORMAT=<ID=SRF,Number=1,Type=Integer,Description=\"Number of reference observations on the forward strand\">" << endl
        //<< "##FORMAT=<ID=SRR,Number=1,Type=Integer,Description=\"Number of reference observations on the reverse strand\">" << endl
        //<< "##FORMAT=<ID=SAF,Number=1,Type=Integer,Description=\"Number of alternate observations on the forward strand\">" << endl
        //<< "##FORMAT=<ID=SAR,Number=1,Type=Integer,Description=\"Number of alternate observations on the reverse strand\">" << endl
        //<< "##FORMAT=<ID=LR,Number=1,Type=Integer,Description=\"Number of reference observations placed left of the loci\">" << endl
        //<< "##FORMAT=<ID=LA,Number=1,Type=Integer,Description=\"Number of alternate observations placed left of the loci\">" << endl
        //<< "##FORMAT=<ID=ER,Number=1,Type=Integer,Description=\"Number of reference observations overlapping the loci in their '3 end\">" << endl
        //<< "##FORMAT=<ID=EA,Number=1,Type=Integer,Description=\"Number of alternate observations overlapping the loci in their '3 end\">" << endl
        << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
        << join(sampleList, "\t") << endl;

    return headerss.str();

}


void AlleleParser::setupVCFOutput(void) {
    string vcfheader = vcfHeader();
    variantCallFile.openForOutput(vcfheader);
}

void AlleleParser::setupVCFInput(void) {
    if (!parameters.variantPriorsFile.empty()) {
        variantCallInputFile.open(parameters.variantPriorsFile);
        currentVariant = new vcf::Variant(variantCallInputFile);

        // get sample names from VCF input file
        //
        // NB, adding this stanza will change the way that the VCF output
        // describes alternates, present observations, etc. so that the samples
        // in the VCF input are also included.  the result is confusing output,
        // but it could be useful in some situations.
        //
        // TODO optionally include this (via command-line parameter)
        //
        //for (vector<string>::iterator s = variantCallInputFile.sampleNames.begin(); s != variantCallInputFile.sampleNames.end(); ++s) {
        //    sampleList.push_back(*s);
        //}

    }
}

void AlleleParser::loadBamReferenceSequenceNames(void) {

    //--------------------------------------------------------------------------
    // read reference sequences from input file
    //--------------------------------------------------------------------------

    // store the names of all the reference sequences in the BAM file
    referenceSequences = bamMultiReader.GetReferenceData();
    int i = 0;
    for (RefVector::iterator r = referenceSequences.begin(); r != referenceSequences.end(); ++r) {
        referenceIDToName[i] = r->RefName;
        ++i;
    }

    DEBUG("Number of ref seqs: " << bamMultiReader.GetReferenceCount());

}


void AlleleParser::loadFastaReference(void) {

    DEBUG("loading fasta reference " << parameters.fasta);

    // This call loads the reference and reads any index file it can find.  If
    // it can't find an index file for the reference, it will attempt to
    // generate one alongside it.  Note that this only loads the reference.
    // Sequence data is obtained by progressive calls to
    // reference.getSubSequence(..), thus keeping our memory requirements low.

    reference.open(parameters.fasta);

}

// alignment-based method for loading the first bit of our reference sequence
void AlleleParser::loadReferenceSequence(BamAlignment& alignment) {
    DEBUG2("loading reference sequence overlapping first alignment");
    currentPosition = alignment.Position;
    currentSequenceStart = alignment.Position;
    currentSequenceName = referenceIDToName[alignment.RefID];
    currentRefID = alignment.RefID;
    DEBUG2("reference.getSubSequence("<< currentSequenceName << ", " << currentSequenceStart << ", " << alignment.AlignedBases.length() << ")");
    currentSequence = uppercase(reference.getSubSequence(currentSequenceName, currentSequenceStart, alignment.Length));
}

// intended to load all the sequence covered by reads which overlap our current target
// this lets us process the reads fully, checking for suspicious reads, etc.
// but does not require us to load the whole sequence
void AlleleParser::loadReferenceSequence(BedTarget* target, int before, int after) {
    basesBeforeCurrentTarget = before;
    basesAfterCurrentTarget = after;
    DEBUG2("loading reference subsequence " << target->seq << " from " << target->left << " - " << before << " to " << target->right << " + " << after << " + before");
    string name = reference.sequenceNameStartingWith(target->seq);
    currentSequence = uppercase(reference.getSubSequence(name, (target->left - 1) - before, (target->right - target->left) + after + before));
    currentReferenceBase = currentReferenceBaseChar();
}

// used to extend the cached reference subsequence when we encounter a read which extends beyond its right bound
void AlleleParser::extendReferenceSequence(int rightExtension) {
    currentSequence += uppercase(reference.getSubSequence(reference.sequenceNameStartingWith(currentSequenceName), 
                                                 currentTarget->right + basesAfterCurrentTarget,
                                                 rightExtension));
    basesAfterCurrentTarget += rightExtension;
}

// maintain a 10bp window around the curent position
// to guarantee our ability to process sequence
void AlleleParser::preserveReferenceSequenceWindow(int bp) {

    // establish left difference between ideal window size and current cached sequence
    int leftdiff = currentSequenceStart - (floor(currentPosition) - bp);
    // guard against falling off the left end of our sequence
    //leftdiff = (floor(currentPosition) - leftdiff < 0) ? floor(currentPosition) : leftdiff;
    leftdiff = (currentSequenceStart - leftdiff < 0) ? currentSequenceStart : leftdiff;
    // right guard is not needed due to the fact that we are just attempting to
    // append sequence substr will simply return fewer additional characters if
    // we go past the right end of the ref
    int rightdiff = (floor(currentPosition) + bp) - (currentSequenceStart + currentSequence.size());

    if (leftdiff > 0) {
        //cerr << currentSequenceStart << endl;
        currentSequence.insert(0, uppercase(reference.getSubSequence(currentSequenceName, currentSequenceStart, leftdiff)));
        currentSequenceStart -= leftdiff;
    }
    if (rightdiff > 0) {
        currentSequence += uppercase(reference.getSubSequence(
                currentSequenceName,
                (currentSequenceStart + currentSequence.size()),
                rightdiff));  // always go 10bp past the end of what we need for alignment registration
    }
}

// ensure we have cached reference sequence according to the current alignment
void AlleleParser::extendReferenceSequence(BamAlignment& alignment) {

    int leftdiff = currentSequenceStart - alignment.Position;
    leftdiff = (currentSequenceStart - leftdiff < 0) ? currentSequenceStart : leftdiff;
    if (leftdiff > 0) {
        currentSequenceStart -= leftdiff;
        if (currentSequenceStart < 0) currentSequenceStart = 0;
        currentSequence.insert(0, uppercase(reference.getSubSequence(currentSequenceName, currentSequenceStart, leftdiff)));
    }

    int rightdiff = (alignment.Position + alignment.AlignedBases.size()) - (currentSequenceStart + currentSequence.size());
    if (rightdiff > 0) {
        currentSequence += uppercase(reference.getSubSequence(
                currentSequenceName,
                (currentSequenceStart + currentSequence.size()),
                rightdiff));
    }

}

void AlleleParser::eraseReferenceSequence(int leftErasure) {
    //cerr << "erasing leftmost " << leftErasure << "bp of cached reference sequence" << endl;
    currentSequence.erase(0, leftErasure);
    currentSequenceStart += leftErasure;
}

void AlleleParser::loadTargets(void) {

    // if we have a targets file, use it...
    // if target file specified use targets from file
    if (!parameters.targets.empty()) {

        DEBUG("Making BedReader object for target file: " << parameters.targets << " ...");

        bedReader.openFile(parameters.targets);

        if (!bedReader.is_open()) {
            ERROR("Unable to open target file: " << parameters.targets << "... terminating.");
            exit(1);
        }

        targets = bedReader.targets;

        if (targets.empty()) {
            ERROR("Could not load any targets from " << parameters.targets);
            exit(1);
        }

        bedReader.close();

        DEBUG("done");

    }

    // if we have a region specified, use it to generate a target
    for (vector<string>::iterator r = parameters.regions.begin(); r != parameters.regions.end(); ++r) {
        // drawn from bamtools_utilities.cpp, modified to suit 1-based context, no end sequence

        string region = *r;
        string startSeq;
        int startPos;
        int stopPos;

        size_t foundFirstColon = region.find(":");

        // we only have a single string, use the whole sequence as the target
        if (foundFirstColon == string::npos) {
            startSeq = region;
            startPos = 0;
            stopPos = -1;
        } else {
            startSeq = region.substr(0, foundFirstColon);
            size_t foundRangeDots = region.find("..", foundFirstColon);
            if (foundRangeDots == string::npos) {
                startPos = atoi(region.substr(foundFirstColon + 1).c_str());
                // differ from bamtools in this regard, in that we process only
                // the specified position if a range isn't given
                stopPos = startPos + 1;
            } else {
                startPos = atoi(region.substr(foundFirstColon + 1, foundRangeDots - foundRangeDots - 1).c_str());
                // if we have range dots specified, but no second number, read to the end of sequence
                if (foundRangeDots + 2 != region.size()) {
                    stopPos = atoi(region.substr(foundRangeDots + 2).c_str()); // end-exclusive, bed-format
                } else {
                    stopPos = reference.sequenceLength(startSeq);
                }
            }
        }

        //DEBUG("startPos == " << startPos);
        //DEBUG("stopPos == " << stopPos);

        // REAL BED format is 0 based, half open (end base not included)
        BedTarget bd(startSeq,
                    (startPos == 0) ? 0 : startPos,
                    (stopPos == -1) ? reference.sequenceLength(startSeq) : stopPos - 1); // internally, we use 0-base inclusive end
        DEBUG("will process reference sequence " << startSeq << ":" << bd.left << ".." << bd.right);
        targets.push_back(bd);
        bedReader.targets.push_back(bd);

    }

    // check validity of targets wrt. reference
    for (vector<BedTarget>::iterator e = targets.begin(); e != targets.end(); ++e) {
        BedTarget& bd = *e;
        if (bd.left < 0 || bd.right - 1 > reference.sequenceLength(bd.seq)) {
            ERROR("Target region coordinates (" << bd.seq << " "
                    << bd.left << " " << bd.right
                    << ") outside of reference sequence bounds ("
                    << bd.seq << " " << reference.sequenceLength(bd.seq) << ") terminating.");
            exit(1);
        }
        if (bd.right < bd.left) {
            ERROR("Invalid target region coordinates (" << bd.seq << " " << bd.left << " " << bd.right << ")"
                    << " right bound is lower than left bound!");
            exit(1);
        }
    }

    bedReader.buildIntervals(); // set up interval tree in the bedreader

    DEBUG("Number of target regions: " << targets.size());

}

void AlleleParser::loadTargetsFromBams(void) {
    // otherwise, if we weren't given a region string or targets file, analyze
    // all reference sequences from BAM file
    DEBUG2("no targets specified, using all targets from BAM files");
    RefVector::iterator refIter = referenceSequences.begin();
    RefVector::iterator refEnd  = referenceSequences.end();
    for( ; refIter != refEnd; ++refIter) {
        RefData refData = *refIter;
        string refName = refData.RefName;
        BedTarget bd(refName, 0, refData.RefLength); // 0-based half open
        DEBUG2("will process reference sequence " << refName << ":" << bd.left << ".." << bd.right);
        targets.push_back(bd);
    }
}

void AlleleParser::loadSampleCNVMap(void) {
    // set default ploidy
    sampleCNV.setDefaultPloidy(parameters.ploidy);

    // load CNV map if provided
    if (!parameters.cnvFile.empty()) {
        if (!sampleCNV.load(parameters.cnvFile)) {
            ERROR("could not load sample map " << parameters.cnvFile << " ... exiting!");
            exit(1);
        }
    }

    // to assert that the reference is haploid, we can iterate through the BAM
    // header to get the reference names and sizes, and then setPloidy on them
    // in the sampleCNV map.  note that the reference "sample" is named after
    // the current reference sequence.
    if (!parameters.diploidReference) {
        for (RefVector::iterator r = referenceSequences.begin(); r != referenceSequences.end(); ++r) {
            sampleCNV.setPloidy(r->RefName, r->RefName, 0, r->RefLength, 1);
        }
    }

}

int AlleleParser::currentSamplePloidy(string const& sample) {
    return sampleCNV.ploidy(sample, currentSequenceName, currentPosition);
}

vector<int> AlleleParser::currentPloidies(Samples& samples) {
    map<int, bool> ploidiesMap;
    vector<int> ploidies;
    for (Samples::iterator s = samples.begin(); s != samples.end(); ++s) {
        string const& name = s->first;
        int samplePloidy = currentSamplePloidy(name);
        ploidiesMap[samplePloidy] = true;
    }
    ploidiesMap[parameters.ploidy] = true;
    for (map<int, bool>::iterator p = ploidiesMap.begin(); p != ploidiesMap.end(); ++p) {
        ploidies.push_back(p->first);
    }
    return ploidies;
}

// meant to be used when we are reading from stdin, to check if we are within targets
bool AlleleParser::inTarget(void) {
    if (targets.empty()) {
        return true;  // everything is in target if we don't have targets
    } else {
        if (bedReader.targetsOverlap(currentSequenceName, currentPosition, currentPosition + 1)) {
            return true;
        } else {
            return false;
        }
    }
}

// initialization function
// sets up environment so we can start registering alleles
AlleleParser::AlleleParser(int argc, char** argv) : parameters(Parameters(argc, argv))
{

    oneSampleAnalysis = false;
    currentRefID = 0; // will get set properly via toNextRefID
    currentTarget = NULL; // to be initialized on first call to getNextAlleles
    currentReferenceAllele = NULL; // same, NULL is brazenly used as an initialization flag
    justSwitchedTargets = false;  // flag to trigger cleanup of Allele*'s and objects after jumping targets
    hasMoreAlignments = true; // flag to track when we run out of alignments in the current target or BAM files
    currentSequenceStart = 0;
    lastHaplotypeLength = 1;
    nullSample = new Sample();


    // initialization
    openTraceFile();
    openFailedFile();
    openOutputFile();

    loadFastaReference();
    // when we open the bam files we can use the number of targets to decide if
    // we should load the indexes
    openBams();
    loadBamReferenceSequenceNames();
    // check how many targets we have specified
    loadTargets();
    getSampleNames();
    getPopulations();
    getSequencingTechnologies();

    // sample CNV
    loadSampleCNVMap();

    // output
    setupVCFOutput();

    // input
    // (now that the VCF file is set up with the samples which are in the input alignments
    // add the samples from the input VCF to the mix)
    setupVCFInput();


}

AlleleParser::~AlleleParser(void) {

    delete nullSample;

    // close trace file?  seems to get closed properly on object deletion...
    if (currentReferenceAllele) delete currentReferenceAllele;

    if (variantCallInputFile.is_open()) delete currentVariant;

}

// position of alignment relative to current sequence
int AlleleParser::currentSequencePosition(const BamAlignment& alignment) {
    return alignment.Position - currentSequenceStart;
}

// TODO clean up this.... just use a string
char AlleleParser::currentReferenceBaseChar(void) {
    return toupper(*currentReferenceBaseIterator());
}

string AlleleParser::currentReferenceBaseString(void) {
    //cerr << currentPosition << " >= " << currentSequenceStart << endl;
    return currentSequence.substr(floor(currentPosition) - currentSequenceStart, 1);
}

string::iterator AlleleParser::currentReferenceBaseIterator(void) {
    return currentSequence.begin() + (floor(currentPosition) - currentSequenceStart);
}

string AlleleParser::referenceSubstr(long int pos, unsigned int len) {
    return uppercase(reference.getSubSequence(currentSequenceName, floor(pos), len));
}

bool AlleleParser::isCpG(string& altbase) {
    // bounds check
    if (floor(currentPosition) - currentSequenceStart - 1 < 0
            || floor(currentPosition) - currentSequenceStart + 1 > currentSequence.size()) {
        return false;
    }
    string prevb = currentSequence.substr(floor(currentPosition) - currentSequenceStart - 1, 1);
    string currb = currentSequence.substr(floor(currentPosition) - currentSequenceStart, 1);
    string nextb = currentSequence.substr(floor(currentPosition) - currentSequenceStart + 1, 1);
    // 5'-3' CpG <-> TpG is represented as CpG <-> CpA in on the opposite strand
    if ((nextb == "G" && ((currb == "C" && altbase == "T") || (currb == "T" && altbase == "C")))
        ||
        (prevb == "C" && ((currb == "G" && altbase == "A") || (currb == "A" && altbase == "G"))))
    {
        return true;
    } else {
        return false;
    }
}

void RegisteredAlignment::addAllele(Allele newAllele, bool mergeComplex, int maxComplexGap) {

    // allele combination rules.  combine the last allele in the list of allele
    // observations according to the following rules
    // 0) reference + SNP, MNP
    // 1) INDEL + (REF <= maxComplexGap) + MNP, INDEL + (REF <= maxComplexGap) + SNP -> complex
    // 2) MNP + SNP, SNP + SNP -> MNP
    // 2) reference + INDEL -> reference.substr(0, reference.size() - 1), reference.at(reference.size()) + INDEL
    assert(newAllele.alternateSequence.size() == newAllele.baseQualities.size());

    alleleTypes |= newAllele.type;

    if (alleles.empty()) {

        // presently, it's unclear how to handle insertions and deletions
        // reported at the beginning of the read.  are these events actually
        // indicative of longer alleles?
        if (!newAllele.isInsertion() && !newAllele.isDeletion() && !newAllele.isNull()) {
            alleles.push_back(newAllele);
        }
        // the same goes for insertions and deletions at the end of reads,
        // these must be dealt with elsewhere

    } else {

        Allele& lastAllele = alleles.back();

        if (newAllele.isReference() && ( newAllele.referenceLength > maxComplexGap || ( newAllele.referenceLength <= maxComplexGap && newAllele.basesRight <= maxComplexGap ))) {
            alleles.push_back(newAllele);
        } else if (lastAllele.isReference()) {
            if (newAllele.isSNP() || newAllele.isMNP() || newAllele.isComplex()) {
                alleles.push_back(newAllele);
            } else if (newAllele.isInsertion() || newAllele.isDeletion()) {
                int p = newAllele.position - 1;
                string seq; vector<pair<int, string> > cig; vector<short> quals;
                lastAllele.subtractFromEnd(1, seq, cig, quals);
                if (lastAllele.length == 0) {
                    alleles.pop_back(); // remove 0-length alleles
                }
                newAllele.addToStart(seq, cig, quals);
                assert(newAllele.position == p);
                alleles.push_back(newAllele);
                assert(newAllele.alternateSequence.size() == newAllele.baseQualities.size());
            }
        } else if (newAllele.isNull()) {
            if (lastAllele.isComplex()) {
                // split apart the last allele if it's 'complex' but followed by a null allele
                vector<pair<int, string> > cigar = splitCigar(lastAllele.cigar);
                if (cigar.back().second == "M") {
                    int matchlen = cigar.back().first;
                    alleles.push_back(lastAllele);
                    Allele& pAllele = alleles.at(alleles.size() - 2);
                    string seq; vector<pair<int, string> > cig; vector<short> quals;
                    pAllele.subtractFromEnd(matchlen, seq, cig, quals);
                    alleles.back().subtractFromStart(pAllele.referenceLength, seq, cig, quals);
                }
            }
            alleles.push_back(newAllele);
        } else {
            // -> complex event, null allele
            if (mergeComplex && lastAllele.position + lastAllele.referenceLength == newAllele.position && !lastAllele.isNull()) {
                lastAllele.mergeAllele(newAllele, ALLELE_COMPLEX);
                assert(lastAllele.alternateSequence.size() == lastAllele.baseQualities.size());
            } else {
                alleles.push_back(newAllele);
            }
        }

    }

}

RegisteredAlignment& AlleleParser::registerAlignment(BamAlignment& alignment, RegisteredAlignment& ra, string& sampleName, string& sequencingTech) {

    string rDna = alignment.QueryBases;
    string rQual = alignment.Qualities;
    int rp = 0;  // read position, 0-based relative to read
    int csp = currentSequencePosition(alignment); // current sequence position, 0-based relative to currentSequence
    int sp = alignment.Position;  // sequence position

#ifdef VERBOSE_DEBUG
    if (parameters.debug2) {
        DEBUG2("registering alignment " << rp << " " << csp << " " << sp << endl <<
                "alignment readName " << alignment.Name << endl <<
                "alignment isPaired " << alignment.IsPaired() << endl <<
                "alignment isMateMapped " << alignment.IsMateMapped() << endl <<
                "alignment isProperPair " << alignment.IsProperPair() << endl <<
                "alignment mapQual " << alignment.MapQuality << endl <<
                "alignment sampleID " << sampleName << endl << 
                "alignment position " << alignment.Position << endl <<
                "alignment length " << alignment.Length << endl <<
                "alignment AlignedBases.size() " << alignment.AlignedBases.size() << endl <<
                "alignment GetEndPosition() " << alignment.GetEndPosition() << endl <<
                "alignment end position " << alignment.Position + alignment.AlignedBases.size());

        stringstream cigarss;
        int alignedLength = 0;
        for (vector<CigarOp>::const_iterator c = alignment.CigarData.begin(); c != alignment.CigarData.end(); ++c) {
            cigarss << c->Type << c->Length;
            if (c->Type == 'D')
                alignedLength += c->Length;
            if (c->Type == 'M')
                alignedLength += c->Length;
        }

        DEBUG2("alignment cigar " << cigarss.str());

        DEBUG2("current sequence pointer: " << csp);

        DEBUG2("read:          " << rDna);
        DEBUG2("aligned bases: " << alignment.AlignedBases);
        DEBUG2("reference seq: " << currentSequence.substr(csp, alignment.AlignedBases.size()));
    }
#endif

    /*
     * The cigar only records matches for sequences that have embedded
     * mismatches.
     *
     * Also, we don't store the entire undelying sequence; just the subsequence
     * that matches our current target region.
     * 
     * As we step through a match sequence, we look for mismatches.  When we
     * see one we set a positional flag indicating the location, and we emit a
     * 'Reference' allele that stretches from the the base after the last
     * mismatch to the base before the current one.
     *
     * An example follows:
     *
     * NNNNNNNNNNNMNNNNNNNNNNNNNNNN
     * reference  ^\-snp  reference
     *
     */

    vector<bool> indelMask (alignment.AlignedBases.size(), false);

    vector<CigarOp>::const_iterator cigarIter = alignment.CigarData.begin();
    vector<CigarOp>::const_iterator cigarEnd  = alignment.CigarData.end();
    for ( ; cigarIter != cigarEnd; ++cigarIter ) {
        unsigned int l = cigarIter->Length;
        char t = cigarIter->Type;
        DEBUG2("cigar item: " << t << l);

        if (t == 'M') { // match or mismatch
            int firstMatch = csp; // track the first match after a mismatch, for recording 'reference' alleles
            int mismatchStart = -1;
            bool inMismatch = false;

            // for each base in the match region
            // increment the csp, sp, and rp
            // if there is a mismatch, record the last matching stretch as a reference allele
            // presently just record one snp per mismatched position, whether or not they are in a series

            for (int i=0; i<l; i++) {

                // extract aligned base
                string b;
                try {
                    b = rDna.at(rp);
                } catch (std::out_of_range outOfRange) {
                    cerr << "Exception: Cannot read past the end of the alignment's sequence." << endl
                         << alignment.AlignedBases << endl
                         << currentSequence.substr(csp, alignment.AlignedBases.size()) << endl
                         << currentSequenceName << ":" << (long unsigned int) currentPosition + 1 << endl;
                    abort();
                }

                // convert base quality value into short int
                long double qual = qualityChar2LongDouble(rQual.at(rp));

                // get reference allele
                string sb;
                try {
                    sb = currentSequence.at(csp);
                } catch (std::out_of_range outOfRange) {
                    cerr << "Exception: Unable to read reference sequence base past end of current cached sequence." << endl
                         << alignment.AlignedBases << endl
                         << currentSequence.substr(csp, alignment.AlignedBases.size()) << endl
                         << currentSequenceName << ":" << (long unsigned int) currentPosition + 1 << endl;
                    abort();
                }

                // record mismatch if we have a mismatch here
                if (b != sb || sb == "N") {  // when the reference is N, we should always call a mismatch
                    if (firstMatch < csp) {
                        // TODO ; verify that the read and reference sequences *do* match
                        int length = csp - firstMatch;
                        //string matchingSequence = currentSequence.substr(csp - length, length);
                        //cerr << "matchingSequence " << matchingSequence << endl;
                        string readSequence = rDna.substr(rp - length, length);
                        //cerr << "readSequence " << readSequence << endl;
                        string qualstr = rQual.substr(rp - length, length);
                        // record 'reference' allele for last matching region
                        if (allATGC(readSequence)) {
                            ra.addAllele(Allele(ALLELE_REFERENCE,
                                        currentSequenceName, sp - length, &currentPosition, &currentReferenceBase, length, 
                                        rp, // bases left (for first base in ref allele)
                                        alignment.QueryBases.size() - rp, // bases right (for first base in ref allele)
                                        readSequence, sampleName, alignment.Name, sequencingTech,
                                        !alignment.IsReverseStrand(), alignment.MapQuality, qualstr,
                                        alignment.MapQuality, alignment.IsPaired(), alignment.IsMateMapped(), alignment.IsProperPair(),
                                        convert(length) + "M",
                                        &ra.alleles),
                                    parameters.allowComplex, parameters.maxComplexGap);
                        }
                    }
                    // register mismatch
                    if (qual >= parameters.BQL2) {
                        ++ra.mismatches;  // increment our mismatch counter if we're over BQL2
                        ++ra.snpCount; // always increment snp counter
                    }

                    // always emit a snp, if we have too many mismatches over
                    // BQL2 then we will discard the registered allele in the
                    // calling context
                    if (!inMismatch) {
                        mismatchStart = csp;
                        inMismatch = true;
                    }
                    firstMatch = csp + 1;

                } else if (inMismatch) {
                    inMismatch = false;
                    int length = csp - mismatchStart;
                    //string matchingSequence = currentSequence.substr(csp - length, length);
                    string readSequence = rDna.substr(rp - length, length);
                    string qualstr = rQual.substr(rp - length, length);
                    AlleleType mismatchtype = (length == 1) ? ALLELE_SNP : ALLELE_MNP;
                    long double lqual = minQuality(qualstr);
                    if (allATGC(readSequence)) {
                        ra.addAllele(Allele(mismatchtype, currentSequenceName, sp - length, &currentPosition,
                                    &currentReferenceBase, length,
                                    rp - length, // bases left
                                    alignment.QueryBases.size() - rp, // bases right
                                    readSequence,
                                    sampleName, alignment.Name, sequencingTech,
                                    !alignment.IsReverseStrand(), lqual,
                                    qualstr, alignment.MapQuality,
                                    alignment.IsPaired(),
                                    alignment.IsMateMapped(),
                                    alignment.IsProperPair(),
                                    convert(length) + "X",
                                    &ra.alleles),
                                parameters.allowComplex, parameters.maxComplexGap);
                    } else {
                        ra.addAllele(Allele(ALLELE_NULL, currentSequenceName, sp - length, &currentPosition,
                                    &currentReferenceBase, length,
                                    rp - length, // bases left
                                    alignment.QueryBases.size() - rp, // bases right
                                    readSequence,
                                    sampleName, alignment.Name, sequencingTech,
                                    !alignment.IsReverseStrand(), lqual, qualstr,
                                    alignment.MapQuality, alignment.IsPaired(),
                                    alignment.IsMateMapped(),
                                    alignment.IsProperPair(),
                                    convert(length) + "N",
                                    &ra.alleles),
                                parameters.allowComplex, parameters.maxComplexGap);
                    }
                }

                // update positions
                ++sp;
                ++csp;
                ++rp;
            }
            // catch mismatches at the end of the match
            if (inMismatch) {
                inMismatch = false;
                int length = csp - mismatchStart;
                //string matchingSequence = currentSequence.substr(csp - length, length);
                string readSequence = rDna.substr(rp - length, length);
                string qualstr = rQual.substr(rp - length, length);
                AlleleType mismatchtype = (length == 1) ? ALLELE_SNP : ALLELE_MNP;
                long double lqual = minQuality(qualstr);
                if (allATGC(readSequence)) {
                    ra.addAllele(Allele(mismatchtype, currentSequenceName, sp - length, &currentPosition,
                                &currentReferenceBase, length,
                                rp - length, // bases left
                                alignment.QueryBases.size() - rp, // bases right
                                readSequence,
                                sampleName, alignment.Name, sequencingTech,
                                !alignment.IsReverseStrand(), lqual, qualstr,
                                alignment.MapQuality, alignment.IsPaired(),
                                alignment.IsMateMapped(),
                                alignment.IsProperPair(),
                                convert(length) + "X",
                                &ra.alleles),
                            parameters.allowComplex, parameters.maxComplexGap);
                } else {
                    ra.addAllele(Allele(ALLELE_NULL, currentSequenceName, sp - length, &currentPosition,
                                &currentReferenceBase, length,
                                rp - length, // bases left
                                alignment.QueryBases.size() - rp, // bases right
                                readSequence,
                                sampleName, alignment.Name, sequencingTech,
                                !alignment.IsReverseStrand(), lqual, qualstr,
                                alignment.MapQuality, alignment.IsPaired(),
                                alignment.IsMateMapped(),
                                alignment.IsProperPair(),
                                convert(length) + "N",
                                &ra.alleles),
                            parameters.allowComplex, parameters.maxComplexGap);
                }
            // or, if we are not in a mismatch, construct the last reference allele of the match
            } else if (firstMatch < csp) {
                int length = csp - firstMatch;
                //string matchingSequence = currentSequence.substr(csp - length, length);
                string readSequence = rDna.substr(rp - length, length);
                string qualstr = rQual.substr(rp - length, length);
                if (allATGC(readSequence)) {
                    ra.addAllele(Allele(ALLELE_REFERENCE,
                            currentSequenceName, sp - length, &currentPosition, &currentReferenceBase, length,
                            rp, // bases left (for first base in ref allele)
                            alignment.QueryBases.size() - rp, // bases right (for first base in ref allele)
                            readSequence, sampleName, alignment.Name, sequencingTech,
                            !alignment.IsReverseStrand(), alignment.MapQuality, qualstr,
                            alignment.MapQuality, alignment.IsPaired(),
                            alignment.IsMateMapped(),
                            alignment.IsProperPair(),
                            convert(length) + "M",
                            &ra.alleles),
                        parameters.allowComplex, parameters.maxComplexGap);
                }
            }
        } else if (t == 'D') { // deletion

            // because deletions have no quality information,
            // use the surrounding sequence quality as a proxy
            // to provide quality scores of equivalent magnitude to insertions,
            // take N bp, right-centered on the position of the deletion
            // this logic prevents overflow of the read
            int spanstart;

            // this is used to calculate the quality string adding 2bp grounds
            // the indel in the surrounding sequence, which it is dependent
            // upon
            int L = l + 2;

            if (L > rQual.size()) {
                L = rQual.size();
                spanstart = 0;
            } else {
                // set lower bound to 0
                if (rp < (L / 2)) {
                    spanstart = 0;
                } else {
                    spanstart = rp - (L / 2);
                }
                // set upper bound to the string length
                if (spanstart + L > rQual.size()) {
                    spanstart = rQual.size() - L;
                }
            }

            string qualstr = rQual.substr(spanstart, L);

            // quality, scaled inversely by the ratio between the quality
            // string length and the length of the event
            long double qual = sumQuality(qualstr);

            // quality adjustment:
            // scale the quality by the inverse harmonic sum of the length of
            // the quality string X a scaling constant derived from the ratio
            // between the length of the quality string and the length of the
            // allele
            //qual += ln2phred(log((long double) l / (long double) L));
            qual += ln2phred(log((long double) L / (long double) l));
            qual /= harmonicSum(l);

            if (qual >= parameters.BQL2) {
                ra.mismatches += l;
                for (int i=0; i<l; i++) {
                    indelMask[sp - alignment.Position + i] = true;
                }
            }

            string refseq = currentSequence.substr(csp, l);
            if (allATGC(refseq)) {
                ra.addAllele(Allele(ALLELE_DELETION,
                        currentSequenceName, sp, &currentPosition, &currentReferenceBase, l,
                        rp, // bases left
                        alignment.QueryBases.size() - rp, // bases right
                        "", sampleName, alignment.Name, sequencingTech,
                        !alignment.IsReverseStrand(), qual, "",
                        alignment.MapQuality, alignment.IsPaired(),
                        alignment.IsMateMapped(), alignment.IsProperPair(),
                        convert(l) + "D",
                        &ra.alleles),
                    parameters.allowComplex, parameters.maxComplexGap);
            }
            ++ra.indelCount;

            sp += l;  // update sample position
            csp += l;

        } else if (t == 'I') { // insertion

            //string qualstr = rQual.substr(rp, l);
            int spanstart;

            // this is used to calculate the quality string adding 2bp grounds
            // the indel in the surrounding sequence, which it is dependent
            // upon
            int L = l + 2;

            if (L > rQual.size()) {
                L = rQual.size();
                spanstart = 0;
            } else {
                // set lower bound to 0
                if (rp < 1) {
                    spanstart = 0;
                } else {
                    spanstart = rp - 1;
                }
                // set upper bound to the string length
                if (spanstart + L > rQual.size()) {
                    spanstart = rQual.size() - L;
                }
            }

            string qualstr = rQual.substr(spanstart, L);

            // quality, scaled inversely by the ratio between the quality
            // string length and the length of the event
            long double qual = sumQuality(qualstr);

            // quality adjustment:
            // scale the quality by the inverse harmonic sum of the length of
            // the quality string X a scaling constant derived from the ratio
            // between the length of the quality string and the length of the
            // allele
            //qual += ln2phred(log((long double) l / (long double) L));
            qual += ln2phred(log((long double) L / (long double) l));
            qual /= harmonicSum(l);

            if (qual >= parameters.BQL2) {
                ra.mismatches += l;
                indelMask[sp - alignment.Position] = true;
            }

            string readseq = rDna.substr(rp, l);
            if (allATGC(readseq)) {
                ra.addAllele(Allele(ALLELE_INSERTION,
                        currentSequenceName, sp, &currentPosition, &currentReferenceBase, l,
                        rp - l, // bases left
                        alignment.QueryBases.size() - rp, // bases right
                        readseq,
                        sampleName, alignment.Name, sequencingTech,
                        !alignment.IsReverseStrand(), qual,
                        rQual.substr(rp, l), alignment.MapQuality, alignment.IsPaired(),
                        alignment.IsMateMapped(), alignment.IsProperPair(),
                        convert(l) + "I",
                        &ra.alleles),
                    parameters.allowComplex, parameters.maxComplexGap);
            }
            ++ra.indelCount;

            rp += l;

        // handle other cigar element types
        } else if (t == 'S') { // soft clip, clipped sequence present in the read not matching the reference
            // skip these bases in the read
            rp += l;// sp += l; csp += l;
        } else if (t == 'H') { // hard clip on the read, clipped sequence is not present in the read
            // the alignment position is the first non-clipped base.
            // thus, hard clipping seems to just be an indicator that we clipped something
            // here we do nothing
            //sp += l; csp += l;
        } else if (t == 'N') { // skipped region in the reference not present in read, aka splice
            sp += l; csp += l;
        }
        // ignore padding
        //} else if (t == 'P') { // padding, silent deletion from the padded reference sequence
        //    sp += l; csp += l;
        //}
    } // end cigar iter loop

    if (ra.alleles.empty()) {
        DEBUG2("generated no alleles from read");
        return ra;
    }

    // this deals with the case in which we have embedded Ns in the read
    // often this happens at the start or end of reads, thus affecting our RegisteredAlignment::start and ::end
    ra.start = ra.alleles.front().position;
    ra.end = ra.alleles.back().position + ra.alleles.back().referenceLength;

    double alignedBases = 0;
    double mismatchCount = 0;
    double matchCount = 0;
    double indelCount = 0;

    // tally mismatches in two categories, gaps and mismatched bases
    for (vector<Allele>::iterator a = ra.alleles.begin(); a != ra.alleles.end(); ++a) {
        Allele& allele = *a;
        switch (allele.type) {
            case ALLELE_REFERENCE:
                alignedBases += allele.length;
                matchCount += allele.length;
                break;
            case ALLELE_SNP:
            case ALLELE_MNP:
                alignedBases += allele.length;
                mismatchCount += allele.length;
                break;
            case ALLELE_INSERTION:
            case ALLELE_DELETION:
            case ALLELE_COMPLEX:
                ++indelCount;
                break;
            default:
                break;
        }
    }

    double mismatchRate = ( indelCount + mismatchCount ) / alignedBases;
    double snpRate = mismatchCount / alignedBases;
    double indelRate = indelCount / alignedBases;

    // store mismatch information about the alignment in the alleles
    // for each allele, normalize the mismatch rates by ignoring that allele,
    // this allows us to relate the mismatch rate without reference to called alleles
    for (vector<Allele>::iterator a = ra.alleles.begin(); a != ra.alleles.end(); ++a) {
        Allele& allele = *a;
        allele.readMismatchRate = mismatchRate;
        allele.readSNPRate = snpRate;
        allele.readIndelRate = indelRate;

        switch (allele.type) {
            case ALLELE_REFERENCE:
                allele.readMismatchRate = mismatchRate;
                allele.readSNPRate = snpRate;
                allele.readIndelRate = indelRate;
                break;
            case ALLELE_SNP:
            case ALLELE_MNP:
                allele.readSNPRate = ( mismatchCount - allele.length ) / alignedBases;
                allele.readIndelRate = indelRate;
                allele.readMismatchRate = indelRate + allele.readSNPRate;
                break;
            case ALLELE_INSERTION:
            case ALLELE_DELETION:
            case ALLELE_COMPLEX:
                allele.readSNPRate = snpRate;
                allele.readIndelRate = ( indelCount - 1 ) / alignedBases;
                allele.readMismatchRate = allele.readIndelRate + snpRate;
                break;
            default:
                break;
        }
    }

    // mark positions in each alignment which are within IDW bases of an indel
    // these are then filtered at each call to getAlleles()
    if (parameters.IDW > -1) { // -1 is the default value and means 'no indel exclusion'
        for (vector<bool>::iterator m = indelMask.begin(); m < indelMask.end(); ++m) {
            if (*m) {
                vector<bool>::iterator q = m - parameters.IDW;
                if (q < indelMask.begin()) q = indelMask.begin();
                for (; q <= m + parameters.IDW && q != indelMask.end(); ++q) {
                    *q = true;
                }
                m += parameters.IDW + 1;
            }
        }
        for (vector<Allele>::iterator a = ra.alleles.begin(); a != ra.alleles.end(); ++a) {
            Allele& allele = *a;
            int startpos = (allele.position - alignment.Position);
            int endpos = startpos + allele.length;
            vector<bool>::iterator im = indelMask.begin();
            // if there is anything masked, store it, otherwise just leave the
            // indelMask on this alignment empty, which means, "no masking" in
            // Allele::masked()
            for (vector<bool>::iterator q = im + startpos; q != im + endpos; ++q) {
                if (*q) { // there is a masked element
                    allele.indelMask.resize(allele.length);
                    copy(im + startpos, im + endpos, allele.indelMask.begin());
                    break;
                }
            }
        }
    }

    // ignore insertions, deletions, and N's which occur at the end of the read with
    // no reference-matching bases before the end of the read
    if (ra.alleles.back().isInsertion() || ra.alleles.back().isDeletion() || ra.alleles.back().isNull()) {
        ra.alleles.pop_back();
    }

#ifdef VERBOSE_DEBUG
    if (parameters.debug2) {
        cerr << "alleles:\n" << join(ra.alleles, "\n");
    }
#endif

    return ra;

}


void AlleleParser::updateAlignmentQueue(void) {

    DEBUG2("updating alignment queue");
    DEBUG2("currentPosition = " << currentPosition 
            << "; currentSequenceStart = " << currentSequenceStart 
            << "; currentSequence end = " << currentSequence.size() + currentSequenceStart);

    // make sure we have sequence for the *first* alignment
    extendReferenceSequence(currentAlignment);

    // push to the front until we get to an alignment that doesn't overlap our
    // current position or we reach the end of available alignments
    // filter input reads; only allow mapped reads with a certain quality
    DEBUG2("currentAlignment.Position == " << currentAlignment.Position 
            << ", currentAlignment.AlignedBases.size() == " << currentAlignment.AlignedBases.size()
            << ", currentPosition == " << currentPosition
            << ", currentSequenceStart == " << currentSequenceStart
            << " .. + currentSequence.size() == " << currentSequenceStart + currentSequence.size()
            );
    if (hasMoreAlignments && currentAlignment.Position <= currentPosition && currentAlignment.RefID == currentRefID) {
        do {
            DEBUG2("currentAlignment.Name == " << currentAlignment.Name);
            // get read group, and map back to a sample name
            string readGroup;
            if (!currentAlignment.GetTag("RG", readGroup)) {
                if (!oneSampleAnalysis) {
                    ERROR("Couldn't find read group id (@RG tag) for BAM Alignment " <<
                            currentAlignment.Name << " at position " << currentPosition
                            << " in sequence " << currentSequence << " EXITING!");
                    exit(1);
                } else {
                    readGroup = "unknown";
                }
            } else {
                if (oneSampleAnalysis) {
                    ERROR("No read groups specified in BAM header, but alignment " <<
                            currentAlignment.Name << " at position " << currentPosition
                            << " in sequence " << currentSequence << " has a read group.");
                    exit(1);
                }
            }

            // skip this alignment if we are not analyzing the sample it is drawn from
            if (readGroupToSampleNames.find(readGroup) == readGroupToSampleNames.end())
                continue;

            // skip this alignment if we are not using duplicate reads (we remove them by default)
            if (currentAlignment.IsDuplicate() && !parameters.useDuplicateReads)
                continue;

            // skip unmapped alignments, as they cannot be used in the algorithm
            if (!currentAlignment.IsMapped())
                continue;

            // skip alignments which are non-primary
            if (!currentAlignment.IsPrimaryAlignment()) // TODO add flag to optionally allow this
                continue;

            // otherwise, get the sample name and register the alignment to generate a sequence of alleles
            // we have to register the alignment to acquire some information required by filters
            // such as mismatches

            // initially skip reads with low mapping quality (what happens if MapQuality is not in the file)
            if (currentAlignment.MapQuality >= parameters.MQL0) {
                // extend our cached reference sequence to allow processing of this alignment
                extendReferenceSequence(currentAlignment);
                // left realign indels
                if (parameters.leftAlignIndels) {
                    int length = currentAlignment.GetEndPosition() - currentAlignment.Position + 1;
                    stablyLeftAlign(currentAlignment,
                        currentSequence.substr(currentSequencePosition(currentAlignment), length));
                }
                // get sample name
                string sampleName = readGroupToSampleNames[readGroup];
                string sequencingTech = readGroupToTechnology[readGroup];
                // decomposes alignment into a set of alleles
                // here we get the deque of alignments ending at this alignment's end position
                deque<RegisteredAlignment>& rq = registeredAlignments[currentAlignment.GetEndPosition()];
                // and insert the registered alignment into that deque
                rq.push_front(RegisteredAlignment(currentAlignment));
                RegisteredAlignment& ra = rq.front();
                registerAlignment(currentAlignment, ra, sampleName, sequencingTech);
                // backtracking if we have too many mismatches
                // or if there are no recorded alleles
                if (ra.alleles.empty()
                        || ((float) ra.mismatches / (float) currentAlignment.QueryBases.size()) > parameters.readMaxMismatchFraction
                        || ra.mismatches > parameters.RMU
                        || ra.snpCount > parameters.readSnpLimit
                        || ra.indelCount > parameters.readIndelLimit) {
                    rq.pop_front(); // backtrack
                } else {
                    // push the alleles into our registered alleles vector
                    for (vector<Allele>::iterator allele = ra.alleles.begin(); allele != ra.alleles.end(); ++allele) {
                        registeredAlleles.push_back(&*allele);
                    }
                }
            }
        } while ((hasMoreAlignments = bamMultiReader.GetNextAlignment(currentAlignment)) && currentAlignment.Position <= currentPosition
                && currentAlignment.RefID == currentRefID);
    }

    DEBUG2("... finished pushing new alignments");

}

// updates registered alleles and erases the unused portion of our cached reference sequence
void AlleleParser::updateRegisteredAlleles(void) {

    long int lowestPosition = currentSequenceStart + currentSequence.size();

    // remove reference alleles which are no longer overlapping the current position
    // http://stackoverflow.com/questions/347441/erasing-elements-from-a-vector
    vector<Allele*>& alleles = registeredAlleles;

    for (vector<Allele*>::iterator allele = alleles.begin(); allele != alleles.end(); ++allele) {
        long unsigned int position = (*allele)->position;
        if (currentPosition >= position + (*allele)->referenceLength) {
            *allele = NULL;
        }
        else {
            if (position < lowestPosition)
                lowestPosition = position;
        }
    }

    alleles.erase(remove(alleles.begin(), alleles.end(), (Allele*)NULL), alleles.end());

    if (lowestPosition <= currentPosition) {
        int diff = lowestPosition - currentSequenceStart;
        if (diff > 0 && currentSequenceStart + diff < currentPosition - CACHED_REFERENCE_WINDOW) {
            eraseReferenceSequence(diff);
        }
    }
}

void AlleleParser::updateInputVariants(void) {

    if (variantCallInputFile.is_open()) {

        if (hasMoreVariants && currentVariant->position - 1 <= currentPosition && currentVariant->sequenceName == currentSequenceName) {
            do {
                DEBUG2("getting input alleles from input VCF at position " << currentVariant->sequenceName << ":" << currentVariant->position);
                long int pos = currentVariant->position - 1;
                // get alternate alleles
                map<string, vector<vcf::VariantAllele> > variantAlleles = currentVariant->parsedAlternates();
                vector< vector<vcf::VariantAllele> > orderedVariantAlleles;
                for (vector<string>::iterator a = currentVariant->alt.begin(); a != currentVariant->alt.end(); ++a) {
                    orderedVariantAlleles.push_back(variantAlleles[*a]);
                }

                vector<Allele> genotypeAlleles;
                set<long int> alternatePositions;

                for (vector< vector<vcf::VariantAllele> >::iterator g = orderedVariantAlleles.begin(); g != orderedVariantAlleles.end(); ++g) {

                    vector<vcf::VariantAllele>& altAllele = *g;

                    vector<Allele> alleles;

                    for (vector<vcf::VariantAllele>::iterator v = altAllele.begin(); v != altAllele.end(); ++v) {

                        vcf::VariantAllele& variant = *v;
                        long int allelePos = variant.position - 1;
                        AlleleType type;
                        string alleleSequence = variant.alt;

                        int len = 0;
                        int reflen = 0;
                        string cigar;

                        // XXX
                        // FAIL
                        // you need to add in the reference bases between the non-reference ones!
                        // to allow for complex events!

                        if (variant.ref == variant.alt) {
                            // XXX note that for reference alleles, we only use the first base internally
                            // but this is technically incorrect, so this hack should be noted
                            len = variant.ref.size();
                            reflen = len;
                            //alleleSequence = alleleSequence.at(0); // take only the first base
                            type = ALLELE_REFERENCE;
                            cigar = convert(len) + "M";
                        } else if (variant.ref.size() == variant.alt.size()) {
                            len = variant.ref.size();
                            reflen = len;
                            if (variant.ref.size() == 1) {
                                type = ALLELE_SNP;
                            } else {
                                type = ALLELE_MNP;
                            }
                            cigar = convert(len) + "X";
                        } else if (variant.ref.size() > variant.alt.size()) {
                            len = variant.ref.size() - variant.alt.size();
                            reflen = len;
                            type = ALLELE_DELETION;
                            cigar = convert(len) + "D";
                        } else {
                            len = variant.alt.size() - variant.ref.size();
                            reflen = 0;
                            type = ALLELE_INSERTION;
                            cigar = convert(len) + "I";
                        }

                        Allele allele = genotypeAllele(type, alleleSequence, (unsigned int) len, cigar, (unsigned int) reflen, allelePos);
                        DEBUG2("input allele: " << allele);
                        alleles.push_back(allele);

                    }

                    // Variant::parsedAlternates() only gives us alternate alleles
                    // for now, add reference sequences back between the alternates here
                    if (alleles.size() > 1) {
                        vector<Allele> newAlleles;
                        vector<Allele>::iterator p = alleles.begin();
                        for (vector<Allele>::iterator a = alleles.begin(); a != alleles.end(); ++a) {
                            if (p != a) {
                                if (p->position + p->referenceLength < a->position) {
                                    // insert a reference allele
                                    long int pend = p->position + p->referenceLength;
                                    string refsequence = reference.getSubSequence(currentVariant->sequenceName, pend, a->position - pend);
                                    string cigar = convert(refsequence.size()) + "M";
                                    Allele refAllele = genotypeAllele(ALLELE_REFERENCE, refsequence, refsequence.size(), cigar, refsequence.size(), pend);
                                    newAlleles.push_back(refAllele);
                                }
                            }
                            newAlleles.push_back(*a);
                            p = a;
                        }
                        alleles = newAlleles;
                    }

                    // for any deletion alleles, grap the previous base (per standards in VCF and the rest of the parsing)
                    vector<Allele>::iterator p = alleles.begin();
                    for (vector<Allele>::iterator a = alleles.begin(); a != alleles.end(); ++a) {
                        if (a->isDeletion()) {
                            if (p != a) {
                                if (p->isReference()) {
                                    string seq; vector<pair<int, string> > cig; vector<short> quals;
                                    p->subtractFromEnd(1, seq, cig, quals);
                                    a->addToStart(seq, cig, quals);
                                }
                                // they will be merged otherwise in the complex merge step below
                            } else {
                                // add the previous reference base
                                vector<short> quals;
                                quals.assign(1, 0);
                                vector<pair<int, string> > cig;
                                cig.push_back(make_pair(1, "M"));
                                string seq = reference.getSubSequence(currentVariant->sequenceName, a->position - 1, 1);
                                a->addToStart(seq, cig, quals);
                            }
                        }
                        p = a;
                    }

                    Allele& allele = alleles.front();
                    if (alleles.size() > 1) {
                        vector<Allele>::iterator a = alleles.begin(); ++a;
                        for (; a != alleles.end(); ++a) {
                            if (a->referenceLength > 0) {
                                allele.mergeAllele(*a, ALLELE_COMPLEX);
                            }
                        }
                    }

                    inputVariantAlleles[allele.position].push_back(allele);
                    genotypeAlleles.push_back(allele);

                    if (allele.position + 1 != currentVariant->position) {
                        cerr << "parsed allele position is not the same as the variant position!" << endl;
                        cerr << *currentVariant << endl;
                        cerr << allele << endl;
                        exit(1);
                    }

                    if (allele.type != ALLELE_REFERENCE) {
                        alternatePositions.insert(allele.position);
                    }

                }

                // store the allele counts, if they are provided
                //
                if (currentVariant->info.find("AC") != currentVariant->info.end()
                    && currentVariant->info.find("AN") != currentVariant->info.end()) {
                    vector<string>& afstrs = currentVariant->info["AC"];
                    // XXX is the reference allele included?
                    vector<Allele>& altAlleles = inputVariantAlleles[currentVariant->position - 1];
                    vector<Allele>::iterator al = altAlleles.begin();
                    assert(altAlleles.size() == afstrs.size());
                    // XXX TODO include the reference allele--- how???
                    // do we add another tag to the VCF, use NS, or some other data?
                    map<Allele, int>& afs = inputAlleleCounts[currentVariant->position - 1];
                    int altcountsum = 0;
                    for (vector<string>::iterator afstr = afstrs.begin(); afstr != afstrs.end(); ++afstr, ++al) {
                        int c; convert(*afstr, c);
                        altcountsum += c;
                        afs[*al] = c;
                    }
                    int numalleles; convert(currentVariant->info["AN"].front(), numalleles);
                    int refcount = numalleles - altcountsum;
                    assert(refcount >= 0);
                    if (refcount > 0) {
                        // make ref allele
                        string ref = currentVariant->ref.substr(0, 1);
                        Allele refAllele = genotypeAllele(ALLELE_REFERENCE, ref, ref.size(), "1M", ref.size(), currentVariant->position - 1);
                        // and also cache the observed count of the reference allele
                        afs[refAllele] = refcount;
                    }
                }

                if (currentVariant->samples.empty())
                    continue;

                if (alternatePositions.size() == 1) {

                    // there can be only one alternate allele position, per immediately previous check
                    long int alternatePosition = *alternatePositions.begin();

                    map<int, vector<Genotype> > genotypesByPloidy;

                    // XXX hard-coded to a specific set of ploidys
                    for (int p = 1; p <= 2; ++p) {
                        vector<Genotype> genotypes = allPossibleGenotypes(p, genotypeAlleles);
                        genotypesByPloidy[p] = genotypes;
                    }

                    map<int, map<Genotype*, int> > vcfGenotypeOrder;
                    for (map<int, vector<Genotype> >::iterator gtg = genotypesByPloidy.begin(); gtg != genotypesByPloidy.end(); ++gtg) {
                        int groupPloidy = gtg->first;
                        vector<Genotype>& genotypes = gtg->second;
                        for (vector<Genotype>::iterator g = genotypes.begin(); g != genotypes.end(); ++g) {
                            Genotype* genotypePtr = &*g;
                            Genotype& genotype = *g;
                            vector<int> gtspec;
                            genotype.relativeGenotype(gtspec, genotypeAlleles);
                            // XXX TODO ... EVIL HACKS
                            if (groupPloidy == 2) {
                                int j = gtspec.front();
                                int k = gtspec.back();
                                vcfGenotypeOrder[groupPloidy][genotypePtr] = (k * (k + 1) / 2) + j;
                            } else if (groupPloidy == 1) {
                                vcfGenotypeOrder[groupPloidy][genotypePtr] = gtspec.front();
                            } else {
                                // XXX TODO ...
                            }
                        }
                    }

                    // get the GLs for each sample, and store for use in later computation
                    for (map<string, map<string, vector<string> > >::iterator s = currentVariant->samples.begin();
                            s != currentVariant->samples.end(); ++s) {

                        string sampleName = s->first;
                        map<string, vector<string> >& sample = s->second;
                        string& gt = sample["GT"].front();
                        map<int, int> genotype = vcf::decomposeGenotype(gt);

                        // default case, give a fixed count for each observed allele
                        int ploidy = 0;
                        for (map<int, int>::iterator g = genotype.begin(); g != genotype.end(); ++g) {
                            ploidy += g->second;
                        }

                        if (ploidy > 2) {
                            cerr << "warning, cannot handle ploidy > 2 for in the input VCF due to limitations "
                                 << "of the VCF specification's definition of Genotype ordering" << endl;
                            continue;
                        }

                        // in the case that we have genotype likelihoods in the VCF
                        if (sample.find("GL") != sample.end()) {
                            vector<string>& gls = sample["GL"];
                            vector<long double> genotypeLikelihoods;
                            genotypeLikelihoods.resize(gls.size());
                            transform(gls.begin(), gls.end(), genotypeLikelihoods.begin(), log10string2ln);

                            // now map the gls into genotype space
                            map<Genotype*, int>& genotypeOrder = vcfGenotypeOrder[ploidy];
                            for (map<Genotype*, int>::iterator gto = genotypeOrder.begin(); gto != genotypeOrder.end(); ++gto) {
                                Genotype& genotype = *gto->first;
                                int order = gto->second;
                                map<string, long double>& sampleGenotypeLikelihoods = inputGenotypeLikelihoods[alternatePosition][sampleName];
                                //cerr << sampleName << ":" << convert(genotype) << ":" << genotypeLikelihoods[order] << endl;
                                sampleGenotypeLikelihoods[convert(genotype)] = genotypeLikelihoods[order];
                            }
                        }
                    }

                } else {
                    /*
                    cerr << "warning, ambiguous VCF record (multiple mixed variant classes), unable to use for genotype likelihood input:" << endl
                        << *currentVariant << endl;
                    for (set<long int>::iterator i = alternatePositions.begin(); i != alternatePositions.end(); ++i) {
                        cerr << "has position: " << *i << endl;
                    }
                         //<< currentVariant->sequenceName << ":" << currentVariant->position << endl;
                    */
                }

            } while ((hasMoreVariants = variantCallInputFile.getNextVariant(*currentVariant))
                    && currentVariant->position - 1 <= currentPosition
                    && currentVariant->sequenceName == currentSequenceName);
        }
    }

}

void AlleleParser::addCurrentGenotypeLikelihoods(map<int, vector<Genotype> >& genotypesByPloidy,
    vector<vector<SampleDataLikelihood> >& sampleDataLikelihoods) {

    // check if there are any genotype likelihoods at the current position
    if (inputGenotypeLikelihoods.find(currentPosition) != inputGenotypeLikelihoods.end()) {

        map<string, map<string, long double> >& inputLikelihoodsBySample = inputGenotypeLikelihoods[currentPosition];

        vector<Genotype*> genotypePtrs;
        for (map<int, vector<Genotype> >::iterator gp = genotypesByPloidy.begin(); gp != genotypesByPloidy.end(); ++gp) {
            vector<Genotype>& genotypes = gp->second;
            for (vector<Genotype>::iterator g = genotypes.begin(); g != genotypes.end(); ++g) {
                genotypePtrs.push_back(&*g);
            }
        }
        // if there are, add them to the sample data likelihoods
        for (map<string, map<string, long double> >::iterator gls = inputLikelihoodsBySample.begin();
                gls != inputLikelihoodsBySample.end(); ++gls) {
            const string& sampleName = gls->first;
            map<string, long double>& likelihoods = gls->second;
            map<Genotype*, long double> likelihoodsPtr;
            for (map<string, long double>::iterator gl = likelihoods.begin(); gl != likelihoods.end(); ++gl) {
                const string& genotype = gl->first;
                long double l = gl->second;
                for (vector<Genotype*>::iterator g = genotypePtrs.begin(); g != genotypePtrs.end(); ++g) {
                    if (convert(**g) == genotype) {
                        likelihoodsPtr[*g] = l;
                    }
                }
            }

            Result sampleData;
            sampleData.name = sampleName;
            // TODO add null sample object to sampleData
            // do you need to????
            for (map<Genotype*, long double>::iterator p = likelihoodsPtr.begin(); p != likelihoodsPtr.end(); ++p) {
                sampleData.push_back(SampleDataLikelihood(sampleName, nullSample, p->first, p->second, 0));
            }
            sortSampleDataLikelihoods(sampleData);
            if (!sampleData.empty()) {
                sampleDataLikelihoods.push_back(sampleData);
            }
        }
    }
}


void AlleleParser::getInputAlleleCounts(vector<Allele>& genotypeAlleles, map<string, int>& inputACs) {
    // are there input ACs?
    //
    // if so, match them to the genotype alleles
    if (inputAlleleCounts.find(currentPosition) != inputAlleleCounts.end()) {
        map<Allele, int>& inputCounts = inputAlleleCounts[currentPosition];
        // XXX NB. We only use ACs for alleles in genotypeAlleles
        for (vector<Allele>::iterator a = genotypeAlleles.begin(); a != genotypeAlleles.end(); ++a) {
            if (inputCounts.find(*a) != inputCounts.end()) {
                inputACs[a->currentBase] = inputCounts[*a];
            }
        }
    }
}

void AlleleParser::removeNonOverlappingAlleles(vector<Allele*>& alleles, int haplotypeLength, bool getAllAllelesInHaplotype) {
    for (vector<Allele*>::iterator a = alleles.begin(); a != alleles.end(); ++a) {
        Allele* allele = *a;
        if (haplotypeLength > 1) {
            if (allele->type == ALLELE_REFERENCE) {
                // does the reference allele overlap the haplotype
                if (getAllAllelesInHaplotype
                        && !(currentPosition <= allele->position && allele->position < currentPosition + haplotypeLength)) {
                    *a = NULL;
                } else if (!(allele->position <= currentPosition
                        && allele->position + allele->referenceLength >= currentPosition + haplotypeLength)) {
                    *a = NULL;
                } else if (currentPosition < allele->position) { // not there yet
                    allele->processed = false;
                    *a = NULL;
                }
            } else { // snps, insertions, deletions
                if (getAllAllelesInHaplotype
                        && !(currentPosition <= allele->position && allele->position < currentPosition + haplotypeLength)) {
                    *a = NULL;
                } else if (!(currentPosition == allele->position && allele->referenceLength == haplotypeLength)) {
                    *a = NULL;
                } else if (currentPosition + haplotypeLength <= allele->position) {
                    allele->processed = false;
                    *a = NULL;
                }
            }
        } else {
            if (allele->type == ALLELE_REFERENCE) {
                if (!(allele->position + allele->referenceLength > currentPosition)) {
                    *a = NULL;
                } else if (currentPosition < allele->position) {
                    allele->processed = false;
                    *a = NULL;
                }
            } else { // snps, insertions, deletions
                if (currentPosition >= allele->position) {
                    *a = NULL;
                } else if (currentPosition < allele->position) {
                    allele->processed = false;
                    *a = NULL;
                }
            }
        }
    }
    alleles.erase(remove(alleles.begin(), alleles.end(), (Allele*)NULL), alleles.end());
}

// removes alleles which are filtered at the current position, and unsets their 'processed' flag so they are later evaluated
void AlleleParser::removeFilteredAlleles(vector<Allele*>& alleles) {
    for (vector<Allele*>::iterator allele = alleles.begin(); allele != alleles.end(); ++allele) {
        if ((*allele)->quality < parameters.BQL0 || (*allele)->masked() || (*allele)->currentBase == "N") {
            (*allele)->processed = false; // force re-processing later
            *allele = NULL;
        }
    }
    alleles.erase(remove(alleles.begin(), alleles.end(), (Allele*)NULL), alleles.end());
}

// steps our position/beddata/reference pointers through all positions in all
// targets, returns false when we are finished
//
// pushes and pulls alignments out of our queue of overlapping alignments via
// updateAlignmentQueue() as we progress
//
// returns true if we still have more targets to process
// false otherwise
bool AlleleParser::toNextTarget(void) {

    DEBUG2("seeking to next target with alignments...");

    // XXX
    // XXX
    // TODO cleanup the logic in this section
    // XXX
    // XXX

    // load first target if we have targets and have not loaded the first
    if (!parameters.useStdin && !targets.empty()) {

        bool ok = false;

        // try to load the first target if we need to
        if (!currentTarget)
            ok = loadTarget(&targets.front()) && getFirstAlignment();

        // step through targets until we get to one with alignments
        while (!ok && currentTarget != &targets.back()) {
            if (!loadTarget(++currentTarget)) {
                continue;
            }
            if (ok = getFirstAlignment()) {
                break;
            }
        }

        if (!ok) return false; // last target and couldn't get alignment

        // XXX hack
        clearRegisteredAlignments();
        currentSequenceStart = currentAlignment.Position;
        currentSequenceName = referenceIDToName[currentAlignment.RefID];
        currentRefID = currentAlignment.RefID;
        currentPosition = (currentPosition < currentAlignment.Position) ? currentAlignment.Position : currentPosition;
        currentSequence = uppercase(reference.getSubSequence(currentSequenceName, currentSequenceStart, currentAlignment.Length));

    // stdin, no targets cases
    } else if (!currentTarget && (parameters.useStdin || targets.empty())) {
        // if we have a target for limiting the analysis, use it
        // this happens when you specify stdin + a region string
        if (!targets.empty()) {
            currentTarget = &targets.front();
            loadTarget(currentTarget);
        }
        if (!getFirstAlignment()) {
            ERROR("Could not get first alignment from target");
            return false;
        }
        clearRegisteredAlignments();
        loadReferenceSequence(currentAlignment); // this seeds us with new reference sequence

    // we've reached the end of file, or stdin
    } else if (parameters.useStdin || targets.empty()) {
        return false;
    }

    DEBUG("getting first variant");
    getFirstVariant();
    DEBUG("got first variant");

    justSwitchedTargets = true;
    return true;

}

// TODO refactor this to allow reading from stdin or reading the whole file
// without loading each sequence as a target
bool AlleleParser::loadTarget(BedTarget* target) {

    currentTarget = target;

    DEBUG("processing target " << currentTarget->desc << " " <<
            currentTarget->seq << " " << currentTarget->left << " " <<
            currentTarget->right);
    DEBUG2("loading target reference subsequence");

    currentSequenceName = currentTarget->seq;

    int refSeqID = bamMultiReader.GetReferenceID(currentSequenceName);

    DEBUG2("reference sequence id " << refSeqID);

    DEBUG2("setting new position " << currentTarget->left);
    currentPosition = currentTarget->left;

    if (!bamMultiReader.SetRegion(refSeqID, currentTarget->left, refSeqID, currentTarget->right - 1)) {
        ERROR("Could not SetRegion to " << currentTarget->seq << ":" << currentTarget->left << ".." << currentTarget->right);
        return false;
    }

    if (variantCallInputFile.is_open()) {
        stringstream r;
        // TODO check that coordinates must be 1-bsaed
        r << currentTarget->seq << ":" << max(0, currentTarget->left - 1) << "-" << currentTarget->right - 1;
        if (!variantCallInputFile.setRegion(r.str())) {
            ERROR("Could not set the region of the variants input file to " <<
                    currentTarget->seq << ":" << currentTarget->left << ".." <<
                    currentTarget->right);
            return false;
        } else {
            DEBUG("set region of variant input file to " << 
                    currentTarget->seq << ":" << currentTarget->left << ".." <<
                    currentTarget->right);
        }
    }

    // now that we've jumped, reset the hasMoreAlignments counter
    hasMoreAlignments = true;

    // same for the variants record
    hasMoreVariants = true;

    DEBUG2("set region");

    return true;

}

bool AlleleParser::getFirstAlignment(void) {

    bool hasAlignments = true;
    if (!bamMultiReader.GetNextAlignment(currentAlignment)) {
        hasAlignments = false;
    } else {
        while (!currentAlignment.IsMapped()) {
            if (!bamMultiReader.GetNextAlignment(currentAlignment)) {
                hasAlignments = false;
                break;
            }
        }
    }

    if (hasAlignments) {
        DEBUG2("got first alignment in target region");
    } else {
        if (currentTarget) {
            ERROR("Could not find any mapped reads in target region " << currentSequenceName << ":" << currentTarget->left << ".." << currentTarget->right);
        } else {
            ERROR("Could not find any mapped reads in target region " << currentSequenceName);
        }
        return false;
    }

    return true;

}

bool AlleleParser::getFirstVariant(void) {

    hasMoreVariants = true;
    if (variantCallInputFile.is_open()) {
        if (!variantCallInputFile.getNextVariant(*currentVariant)) {
            hasMoreVariants = false;
        }

        if (hasMoreVariants) {
            DEBUG2("got first variant in target region");
        } else {
            return false;
        }
    }

    return true;

}

void AlleleParser::clearRegisteredAlignments(void) {
    DEBUG2("clearing registered alignments and alleles");
    registeredAlignments.clear();
    registeredAlleles.clear();
}

// TODO
// this should be simplified
// there are two modes of operation
// that in which we have targets
// and that without
//
// if we have targets, we need to keep track of which we're in
// and if we're outside of it, try to get to the next one
// and, if we have targets, we will try to jump around the bam file
//
// if we don't have targets we will just GetNextAlignment until we can't
// anymore.  all positionality of the parser will respond to input alignments.
//
// rewrite things so that we aren't strung out between 8 functions
//

// stepping
//
// if the next position is outside of target region
// seek to next target which is in-bounds for its sequence
// if none exist, return false
//
bool AlleleParser::toNextPosition(void) {

    if (currentSequenceName.empty()) {
        DEBUG2("loading first target");
        if (!toNextTarget()) {
            return false;
        }
    } 
    else {
        ++currentPosition;
    }

    if (!targets.empty() && (
                (!parameters.allowIndels && currentPosition >= currentTarget->right)
                || currentPosition > currentTarget->right - 1)) { // time to move to a new target
        DEBUG("next position " << (long int) currentPosition + 1 <<  " outside of current target right bound " << currentTarget->right + 1);
        if (!toNextTarget()) {
            DEBUG("no more targets, finishing");
            return false;
        }
    }

    // stdin, no targets case
    if (parameters.useStdin || targets.empty()) {
        // TODO rectify this with the other copies of this stanza...
        // implicit step of target sequence
        // XXX this must wait for us to clean out all of our alignments at the end of the target
        while (hasMoreAlignments && !currentAlignment.IsMapped()) {
            hasMoreAlignments = bamMultiReader.GetNextAlignment(currentAlignment);
        }
        if (hasMoreAlignments && registeredAlignments.empty() && currentRefID != currentAlignment.RefID) {
            clearRegisteredAlignments();
            loadReferenceSequence(currentAlignment);
            justSwitchedTargets = true;
        } else if (!hasMoreAlignments && registeredAlignments.empty()) {
            DEBUG("no more alignments in input");
            return false;
        } else if (!hasMoreAlignments && currentPosition > currentSequence.size() + currentSequenceStart) {
            DEBUG("at end of sequence");
            return false;
        }
    }


    currentReferenceBase = currentReferenceBaseChar();

    // handle the case in which we don't have targets but in which we've switched reference sequence

    DEBUG2("processing position " << (long unsigned int) currentPosition + 1 << " in sequence " << currentSequenceName);
    updateAlignmentQueue();
    DEBUG2("updating variants");
    updateInputVariants();
    DEBUG2("updating registered alleles");
    updateRegisteredAlleles(); // this removes unused left-flanking sequence
    //DEBUG2("updating prior variant alleles");
    //updatePriorAlleles();

    // if we have alignments which ended at the previous base, erase them and their alleles
    // TODO check that this doesn't leak...
    DEBUG2("erasing old registered alignments");
    map<long unsigned int, deque<RegisteredAlignment> >::iterator f = registeredAlignments.find(currentPosition - 2);
    if (f != registeredAlignments.end()) {
        registeredAlignments.erase(f);
    }

    // and do the same for the variants from the input VCF
    DEBUG2("erasing old input variant alleles");
    map<long int, vector<Allele> >::iterator v = inputVariantAlleles.find(currentPosition - 3);
    if (v != inputVariantAlleles.end()) {
        inputVariantAlleles.erase(v);
    }

    DEBUG2("erasing old genotype likelihoods");
    map<long int, map<string, map<string, long double> > >::iterator l = inputGenotypeLikelihoods.find(currentPosition - 3);
    if (l != inputGenotypeLikelihoods.end()) {
        inputGenotypeLikelihoods.erase(l);
    }

    DEBUG2("erasing old allele frequencies");
    map<long int, map<Allele, int> >::iterator af = inputAlleleCounts.find(currentPosition - 3);
    if (af != inputAlleleCounts.end()) {
        inputAlleleCounts.erase(af);
    }

    // so we have to make sure it's still there (this matters in low-coverage)
    DEBUG2("updating reference sequence cache");
    preserveReferenceSequenceWindow(CACHED_REFERENCE_WINDOW);

    return true;

}

// XXX for testing only, steps targets but does nothing
bool AlleleParser::dummyProcessNextTarget(void) {

    if (!toNextTarget()) {
        DEBUG("no more targets, finishing");
        return false;
    }

    while (bamMultiReader.GetNextAlignment(currentAlignment)) {
    }

    return true;
}

// adjusts the registered alignment and contained alleles so that one allele
// covers the entire haplotype window
// returns a vector of pointers to alleles generated in this process
// alleles which are discarded are not explicitly removed, but 'squashed',
// which triggers their collection later
bool RegisteredAlignment::fitHaplotype(int haplotypeStart, int haplotypeLength, Allele*& aptr) {

    // if the read overlaps the haplotype window,
    // generate one Allele to describe the read in that region
    // and "squash" the unused ones
    vector<Allele*> newAllelesPtr;
    vector<Allele> newAlleles;

    int haplotypeEnd = haplotypeStart + haplotypeLength;
    
    //if (containedAlleleTypes == ALLELE_REFERENCE) {
    //    return false;
    //}
    //cerr << "registered alignment alleles," << endl << alleles << endl;
    //cerr << "start: " << start << " end: " << end << endl;
    //cerr << "haplotypestart: " << haplotypeStart << " haplotypeend: " << haplotypeEnd << endl;

    if (start <= haplotypeStart && end >= haplotypeEnd) {
        vector<Allele>::iterator a = alleles.begin();
        //cerr << "trying to find overlapping haplotype alleles for the range " << haplotypeStart << " to " << haplotypeEnd << endl;
        while (a + 1 != alleles.end() && a->position + a->referenceLength <= haplotypeStart) {
            ++a;
        }
        vector<Allele>::iterator b = a;
        while (b + 1 != alleles.end() && b->position + b->referenceLength < haplotypeEnd) {
            ++b;
        }

        // do not attempt to build haplotype alleles where there are non-contiguous reads
        for (vector<Allele>::iterator p = alleles.begin(); p != alleles.end(); ++p) {
            if (p != alleles.begin()) {
                if (p->position != (p - 1)->position + (p - 1)->referenceLength) {
                    return true;
                }
            }
        }

        if (a == b && a->isReference()) {
            // no change needed when we just have a reference allele
            return true;
        }

        string seq;
        vector<pair<int, string> > cigar;
        vector<short> quals;

        // now "a" should overlap the start of the haplotype block, and "b" the end
        //cerr << "block start overlaps: " << *a << endl;
        //cerr << "block end overlaps: " << *a << endl;

        // adjust a to match the start of the haplotype block
        if (a->position == haplotypeStart) {
            // nothing to do!
        } else if (a->position < haplotypeStart) {
            // squeeze bases off the front of this allele onto the last allele
            // generating a new allele if there isn't one
            Allele newAllele = *a;
            newAllele.subtractFromEnd(a->position + a->referenceLength - haplotypeStart, seq, cigar, quals);
            a->subtractFromStart(haplotypeStart - a->position, seq, cigar, quals);
            newAlleles.push_back(newAllele);
        }

        if (b->position + b->referenceLength == haplotypeEnd) {
            // nothing to do!!!!
        } else if (b->position + b->referenceLength > haplotypeEnd) {
            Allele newAllele = *b;
            newAllele.subtractFromStart(haplotypeEnd - b->position, seq, cigar, quals);
            b->subtractFromEnd(b->position + b->referenceLength - haplotypeEnd, seq, cigar, quals);
            newAlleles.push_back(newAllele);
        }

        // now, for everything between a and b, merge them into one allele
        while (a != b) {
            vector<pair<int, string> > cigarV = splitCigar(a->cigar);
            vector<Allele>::iterator p = a + 1;
            // update the quality of the merged allele in the same way as we do
            // for complex events
            if (!a->isReference())  {
                p->quality = min(a->quality, p->quality);  // note that phred and log are inverted
                p->lnquality = max(a->lnquality, p->lnquality);
            }
            p->addToStart(a->alternateSequence, cigarV, a->baseQualities);
            a->squash();
            ++a;
        }

        // remove any 0-length alleles, these are useless
        alleles.erase(remove_if(alleles.begin(), alleles.end(), isEmptyAllele), alleles.end());

        for (vector<Allele>::iterator p = newAlleles.begin(); p != newAlleles.end(); ++p) {
            alleles.push_back(*p);
        }

        AllelePositionCompare apcomp;
        sort(alleles.begin(), alleles.end(), apcomp);

        // now the pointers have changed, so find the allele we want... again!!!!!!
        //cerr << "registered alignment alleles, after haplotype construction," << endl << alleles << endl;
        bool hasHaplotypeAllele = false;
        for (vector<Allele>::iterator p = alleles.begin(); p != alleles.end(); ++p) {
            //cerr << *p << endl;
            if (p->position == haplotypeStart && p->position + p->referenceLength == haplotypeEnd) {
                aptr = &*p;
		if (isDividedIndel(*p)) {
		    hasHaplotypeAllele = false;
		} else {
		    hasHaplotypeAllele = true;
		}
                break;
            }
        }

        if (hasHaplotypeAllele) {
            return true;
        } else {
	    return false;
            //assert(hasHaplotypeAllele);
        }

    } else {
        return true;
    }

}

void AlleleParser::buildHaplotypeAlleles(vector<Allele>& alleles, Samples& samples, map<string, vector<Allele*> >& alleleGroups, int allowedAlleleTypes) {

    int haplotypeLength = 1;
    for (vector<Allele>::iterator a = alleles.begin(); a != alleles.end(); ++a) {
        Allele& allele = *a;
        if (!allele.isReference() && allele.referenceLength > haplotypeLength) {
            haplotypeLength = allele.referenceLength;
        }
    }

    if (haplotypeLength > 1) {

        int oldHaplotypeLength = haplotypeLength;
        do {
            oldHaplotypeLength = haplotypeLength;
            //cerr << "haplotype length is " << haplotypeLength << endl;
            getAlleles(samples, allowedAlleleTypes, haplotypeLength, true);
            alleleGroups.clear();
            groupAlleles(samples, alleleGroups);
            alleles = genotypeAlleles(alleleGroups, samples, parameters.onlyUseInputAlleles);
            for (vector<Allele>::iterator a = alleles.begin(); a != alleles.end(); ++a) {
                Allele& allele = *a;
                if (!allele.isReference() && allele.position + allele.referenceLength > currentPosition + haplotypeLength) {
                    haplotypeLength = (allele.position + allele.referenceLength) - currentPosition;
                }
            }
        } while (haplotypeLength != oldHaplotypeLength);

        // for each non-reference allele within the haplotype length of this
        // position, adjust the length and reference sequences of the adjacent
        // alleles 
        DEBUG2("fitting haplotype block " << currentPosition << " to " << currentPosition + haplotypeLength << ", " << haplotypeLength << "bp");

        registeredAlleles.clear();
        samples.clear();

        vector<Allele*> newAlleles;
        for (map<long unsigned int, deque<RegisteredAlignment> >::iterator ras = registeredAlignments.begin(); ras != registeredAlignments.end(); ++ras) {
            deque<RegisteredAlignment>& rq = ras->second;
            for (deque<RegisteredAlignment>::iterator rai = rq.begin(); rai != rq.end(); ++rai) {
                RegisteredAlignment& ra = *rai;
                vector<Allele*> oldptrs;
                for (vector<Allele>::iterator a = ra.alleles.begin(); a != ra.alleles.end(); ++a) {
                    oldptrs.push_back(&*a);
                }
                Allele* aptr;
                //cerr << "before fitting haplotype, we have these alleles:" << endl << ra.alleles << endl;
                if (ra.fitHaplotype(currentPosition, haplotypeLength, aptr)) {
                    //registeredAlleles.insert(registeredAlleles.end(), generatedAlleles.begin(), generatedAlleles.end());
		    // for debugging:
		    //assert(aptr->position == currentPosition);
                    //assert(aptr->referenceLength == haplotypeLength);
                    //cerr << "generated haplotype-matching allele: " << *aptr << endl;
                    //cerr << "and these alleles, " << ra.alleles << endl;

		    for (vector<Allele>::iterator a = ra.alleles.begin(); a != ra.alleles.end(); ++a) {
			a->processed = false; // re-trigger use of all alleles
			registeredAlleles.push_back(&*a);
		    }
                }
            }
        }

        updateRegisteredAlleles();

        // now re-get the alleles
        getAlleles(samples, allowedAlleleTypes, haplotypeLength);

        // re-group the alleles using groupAlleles()
        alleleGroups.clear();
        groupAlleles(samples, alleleGroups);

        // are there two alleles with the same alt sequence?
        // if so, homogenize them, and then re-sort the alleles
        map<string, int> altseqCounts;
        bool multipleAllelesWithIdenticalAlts = false;
        for (map<string, vector<Allele*> >::iterator a = alleleGroups.begin(); a != alleleGroups.end(); ++a) {
            Allele& allele = *a->second.front();
            if (allele.isReference()) {
                continue;
            }
            if (altseqCounts[allele.alternateSequence] > 0) {
                multipleAllelesWithIdenticalAlts = true;
                break;
            } else {
                ++altseqCounts[allele.alternateSequence];
            }
        }

        if (multipleAllelesWithIdenticalAlts) {
            homogenizeAlleles(alleleGroups);
            getAlleles(samples, allowedAlleleTypes, haplotypeLength);
            alleleGroups.clear();
            groupAlleles(samples, alleleGroups);
            alleles = genotypeAlleles(alleleGroups, samples, parameters.onlyUseInputAlleles);
        } else {
            // re-run genotypeAlleles() to update the set of genotype alleles,
            // which is passed by reference
            alleles = genotypeAlleles(alleleGroups, samples, parameters.onlyUseInputAlleles);
        }

    }

    lastHaplotypeLength = haplotypeLength;

}

bool AlleleParser::getNextAlleles(Samples& samples, int allowedAlleleTypes) {
    long int nextPosition = currentPosition + lastHaplotypeLength;
    lastHaplotypeLength = 1;
    while (currentPosition < nextPosition) {
        if (!toNextPosition()) {
            return false;
        } else {
            getAlleles(samples, allowedAlleleTypes);
        }
    }
    return true;
}

void AlleleParser::getAlleles(Samples& samples, int allowedAlleleTypes, int haplotypeLength, bool getAllAllelesInHaplotype) {

    DEBUG2("getting alleles");

    // if we just switched targets, clean up everything in our input vector
    if (justSwitchedTargets) {
        DEBUG2("just switched targets, cleaning up sample alleles");
        for (Samples::iterator s = samples.begin(); s != samples.end(); ++s)
            s->second.clear();
        justSwitchedTargets = false; // TODO this whole flagged stanza is hacky;
                                     // to clean up, store the samples map in
                                     // the AlleleParser and clear it when we jump
    } else {
        // otherwise, update and remove non-overlapping and filtered alleles
        for (Samples::iterator s = samples.begin(); s != samples.end(); ++s) {
            Sample& sample = s->second;
            for (Sample::iterator g = sample.begin(); g != sample.end(); ++g) {
                vector<Allele*>& alleles = g->second;
                removeNonOverlappingAlleles(alleles, haplotypeLength, getAllAllelesInHaplotype); // removes alleles which no longer overlap our current position
                updateAllelesCachedData(alleles);  // calls allele.update() on each Allele*
                removeFilteredAlleles(alleles); // removes alleles which are filtered at this position, 
                                                  // and requeues them for processing by unsetting their 'processed' flag
            }
            sample.sortReferenceAlleles();  // put reference alleles into the correct group, according to their current base
        }
    }

    // if we have targets and are outside of the current target, don't return anything

    // add the reference allele to the analysis
    if (parameters.useRefAllele) {
        if (currentReferenceAllele) delete currentReferenceAllele; // clean up after last position
        currentReferenceAllele = referenceAllele(parameters.MQR, parameters.BQR);
        samples[currentSequenceName].clear();
        samples[currentSequenceName][currentReferenceAllele->currentBase].push_back(currentReferenceAllele);
        //alleles.push_back(currentReferenceAllele);
    }

    // get the variant alleles *at* the current position
    // and the reference alleles *overlapping* the current position
    for (vector<Allele*>::const_iterator a = registeredAlleles.begin(); a != registeredAlleles.end(); ++a) {
        Allele& allele = **a;
        if (!allele.processed
                && allowedAlleleTypes & allele.type
                && ((haplotypeLength > 1 &&
                    ((allele.type == ALLELE_REFERENCE
                      && allele.position <= currentPosition 
                      && allele.position + allele.referenceLength >= currentPosition + haplotypeLength)
                      || 
                     (allele.position == currentPosition
                      && allele.referenceLength == haplotypeLength)
                      ||
                     (getAllAllelesInHaplotype
                      && allele.type != ALLELE_REFERENCE
                      && allele.position >= currentPosition
                      && allele.position < currentPosition + haplotypeLength)))
                  ||
                    (haplotypeLength == 1 &&
                     ((allele.type == ALLELE_REFERENCE
                      && allele.position <= currentPosition
                      && allele.position + allele.referenceLength > currentPosition)
                      || 
                     (allele.position == currentPosition)))
                   ) ) {
            allele.update();
            if (allele.quality >= parameters.BQL0 && !allele.masked() && allele.currentBase != "N"
                    && (allele.isReference() || !allele.alternateSequence.empty())) { // filters haplotype construction chaff
                samples[allele.sampleID][allele.currentBase].push_back(*a);
                // XXX testing
                if (!getAllAllelesInHaplotype) {
                    allele.processed = true;
                    if (haplotypeLength > 1) {
                        if (!allele.isReference() && !(allele.position == currentPosition && allele.referenceLength == haplotypeLength)) {
                            cerr << "non-reference allele should not be added to result alleles because it does not match the haplotype!:" << endl;
                            cerr << "haplotype is from " << currentPosition << " to " << currentPosition + haplotypeLength << ", " << haplotypeLength << "bp" << endl;
                            cerr << allele << endl;
                            assert(false);
                        }
                    }
                }
            }
        }
    }

    vector<string> samplesToErase;
    // now remove empty alleles from our return so as to not confuse processing
    for (Samples::iterator s = samples.begin(); s != samples.end(); ++s) {

        const string& name = s->first;
        Sample& sample = s->second;

        // move updated reference alleles to the right bin
        // everything else will get axed
        //sample.sortReferenceAlleles();

        bool empty = true;
        vector<string> genotypesToErase;
        // and remove any empty groups which remain
        for (Sample::iterator g = sample.begin(); g != sample.end(); ++g) {
            if (g->second.empty()) {
                //cerr << "sample " << name << " has an empty " << g->first << endl;
                //sample.erase(g);
                genotypesToErase.push_back(g->first);
            } else {
                // accumulate bitmap of unique types
                empty = false;
            }
        }

        for (vector<string>::iterator gt = genotypesToErase.begin(); gt != genotypesToErase.end(); ++gt) {
            sample.erase(*gt);
        }

        // and remove the entire sample if it has no alleles
        if (empty || currentSamplePloidy(name) == 0) {
            samplesToErase.push_back(name);
        }
    }

    for (vector<string>::iterator name = samplesToErase.begin(); name != samplesToErase.end(); ++name) {
        samples.erase(*name);
    }

    DEBUG2("done getting alleles");

}

Allele* AlleleParser::referenceAllele(int mapQ, int baseQ) {
    string base = currentReferenceBaseString();
    //string name = reference.filename;
    string name = currentSequenceName; // this behavior matches old bambayes
    string sequencingTech = "reference";
    string baseQstr = "";
    //baseQstr += qualityInt2Char(baseQ);
    Allele* allele = new Allele(ALLELE_REFERENCE, 
            currentSequenceName,
            currentPosition,
            &currentPosition, 
            &currentReferenceBase,
            1, 0, 0, base, name, name, sequencingTech,
            true, baseQ,
            baseQstr,
            mapQ,
            false, false, false, "1M", NULL); // pair information
    allele->genotypeAllele = true;
    allele->baseQualities.push_back(baseQ);
    allele->update();
    return allele;
}

vector<Allele> AlleleParser::genotypeAlleles(
        map<string, vector<Allele*> >& alleleGroups, // alleles grouped by equivalence
        Samples& samples, // alleles grouped by sample
        bool useOnlyInputAlleles
        ) {

    vector<pair<Allele, int> > unfilteredAlleles;

    DEBUG2("getting genotype alleles");

    for (map<string, vector<Allele*> >::iterator group = alleleGroups.begin(); group != alleleGroups.end(); ++group) {
        // for each allele that we're going to evaluate, we have to have at least one supporting read with
        // map quality >= MQL1 and the specific quality of the allele has to be >= BQL1
        DEBUG("allele group " << group->first);
        vector<Allele*>& alleles = group->second;
        if (alleles.size() < parameters.minAltTotal) {
            DEBUG2("allele group lacks sufficient observations in the whole population (min-alternate-total)");
            continue;
        }
        bool passesFilters = false;
        int qSum = 0;
        for (vector<Allele*>::iterator a = alleles.begin(); a != alleles.end(); ++a) {
            DEBUG2("allele " << **a);
            Allele& allele = **a;
            qSum += allele.quality;
            if (!passesFilters && allele.quality >= parameters.BQL1 && allele.mapQuality >= parameters.MQL1) {
                passesFilters = true;
            }
        }
        if (passesFilters) {
            Allele& allele = *(alleles.front());
            int length = allele.length;
            int reflength = allele.referenceLength;
            string altseq = allele.alternateSequence;
            if (allele.type == ALLELE_REFERENCE) {
                length = 1;
                reflength = 1;
                altseq = currentReferenceBase;
            }
            unfilteredAlleles.push_back(make_pair(genotypeAllele(allele.type, altseq, length, allele.cigar, reflength, allele.position), qSum));
        }
    }
    DEBUG2("found genotype alleles");

    map<Allele, int> filteredAlleles;

    DEBUG2("filtering genotype alleles which are not supported by at least " << parameters.minAltCount 
            << " observations comprising at least " << parameters.minAltFraction << " of the observations in a single individual");
    for (vector<pair<Allele, int> >::iterator p = unfilteredAlleles.begin();
            p != unfilteredAlleles.end(); ++p) {

        Allele& genotypeAllele = p->first;
        int qSum = p->second;
        DEBUG2("genotype allele: " << genotypeAllele << " qsum " << qSum);

        for (Samples::iterator s = samples.begin(); s != samples.end(); ++s) {
            Sample& sample = s->second; 
            int alleleCount = 0;
            int qsum = 0;
            Sample::iterator c = sample.find(genotypeAllele.currentBase);
            if (c != sample.end()) {
                vector<Allele*>& obs = c->second;
                alleleCount = obs.size();
                for (vector<Allele*>::iterator a = obs.begin(); a != obs.end(); ++a) {
                    Allele& allele = **a;
                    qsum += allele.quality;
                }
            }
            int observationCount = sample.observationCount();
            if (qsum >= parameters.minAltQSum
                    && alleleCount >= parameters.minAltCount 
                    && ((float) alleleCount / (float) observationCount) >= parameters.minAltFraction) {
                DEBUG(genotypeAllele << " has support of " << alleleCount 
                    << " in individual " << s->first << " and fraction " 
                    << (float) alleleCount / (float) observationCount);
                filteredAlleles[genotypeAllele] = qSum;
                break;
                //out << *genotypeAllele << endl;
            }
        }
    }
    DEBUG2("filtered genotype alleles");


    vector<Allele> resultAlleles;
    vector<Allele> resultIndelAndMNPAlleles;

    string refBase = currentReferenceBaseString();

    if (parameters.useBestNAlleles == 0) {
        // this means "use everything"
        bool hasRefAllele = false;
        for (map<Allele, int>::iterator p = filteredAlleles.begin();
                p != filteredAlleles.end(); ++p) {
            if (p->first.currentBase == refBase)
                hasRefAllele = true;
            DEBUG2("adding allele to result alleles " << p->first.currentBase);
            resultAlleles.push_back(p->first);
        }
        // and add the reference allele if we need it
        if (parameters.forceRefAllele && !hasRefAllele) {
            DEBUG2("including reference allele");
            resultAlleles.insert(resultAlleles.begin(), genotypeAllele(ALLELE_REFERENCE, refBase, 1, "1M", 1, currentPosition));
        }
    } else {
        // this means, use the N best
        vector<pair<Allele, int> > sortedAlleles;
        for (map<Allele, int>::iterator p = filteredAlleles.begin();
                p != filteredAlleles.end(); ++p) {
            sortedAlleles.push_back(make_pair(p->first, p->second));
        }
        DEBUG2("sorting alleles to get best alleles");
        AllelePairIntCompare alleleQualityCompare;
        sort(sortedAlleles.begin(), sortedAlleles.end(), alleleQualityCompare);

        DEBUG2("getting " << parameters.useBestNAlleles << " best SNP alleles, and all other alleles");
        bool hasRefAllele = false;
        for (vector<pair<Allele, int> >::iterator a = sortedAlleles.begin(); a != sortedAlleles.end(); ++a) {
            Allele& allele = a->first;
            if (allele.currentBase == refBase) {
                hasRefAllele = true;
            }
            if (allele.type & (ALLELE_DELETION | ALLELE_INSERTION | ALLELE_MNP | ALLELE_COMPLEX)) {
                DEBUG2("adding allele to result alleles " << allele.currentBase);
                resultIndelAndMNPAlleles.push_back(allele);
            } else {
                DEBUG2("adding allele to SNP alleles " << allele.currentBase);
                resultAlleles.push_back(allele);
            }
            DEBUG2("allele quality sum " << a->second);
        }
        DEBUG2("found " << sortedAlleles.size() << " SNP/ref alleles of which we now have " << resultAlleles.size() << endl
               << "and " << resultIndelAndMNPAlleles.size() << " INDEL and MNP alleles");

        // if we have reached the limit of allowable alleles, and still
        // haven't included the reference allele, include it
        if (parameters.forceRefAllele && !hasRefAllele) {
            DEBUG2("including reference allele in analysis");
            resultAlleles.insert(resultAlleles.begin(), genotypeAllele(ALLELE_REFERENCE, refBase, 1, "1M", 1, currentPosition));
        }

        // if we now have too many alleles (most likely one too many), get rid of some
        while (resultAlleles.size() > parameters.useBestNAlleles) {
            resultAlleles.pop_back();
        }

        // drop the SNPs back into the set of alleles
        for (vector<Allele>::iterator a = resultIndelAndMNPAlleles.begin(); a != resultIndelAndMNPAlleles.end(); ++a) {
            resultAlleles.push_back(*a);
        }

    }

    // now add in the alleles from the input variant set

    if (useOnlyInputAlleles)
        resultAlleles.clear();

    map<long int, vector<Allele> >::iterator v = inputVariantAlleles.find(currentPosition);
    if (v != inputVariantAlleles.end()) {
        vector<Allele>& inputalleles = v->second;
        for (vector<Allele>::iterator a = inputalleles.begin(); a != inputalleles.end(); ++a) {
            Allele& allele = *a;
            // check if the allele is already present
            bool alreadyPresent = false;
            for (vector<Allele>::iterator r = resultAlleles.begin(); r != resultAlleles.end(); ++r) {
                if (r->equivalent(allele)) {
                    alreadyPresent = true;
                    break;
                }
            }
            if (!alreadyPresent) {
                resultAlleles.push_back(allele);
            }
        }
    }

    DEBUG2("found " << resultAlleles.size() << " result alleles");
    return resultAlleles;

}

// homopolymer run length.  number of consecutive nucleotides (prior to this
// position) in the genome reference sequence matching the alternate allele,
// after substituting the alternate in place of the reference sequence allele
int AlleleParser::homopolymerRunLeft(string altbase) {

    int position = currentPosition - 1;
    int sequenceposition = position - currentSequenceStart;
    int runlength = 0;
    while (sequenceposition >= 0 && currentSequence.substr(sequenceposition, 1) == altbase) {
        ++runlength;
        --position;
        sequenceposition = position - currentSequenceStart;
    }
    return runlength;

}

int AlleleParser::homopolymerRunRight(string altbase) {

    int position = currentPosition + 1;
    int sequenceposition = position - currentSequenceStart;
    int runlength = 0;
    while (sequenceposition >= 0 && currentSequence.substr(sequenceposition, 1) == altbase) {
        ++runlength;
        ++position;
        sequenceposition = position - currentSequenceStart;
    }
    return runlength;

}

// returns the number of repeats for each subsequence at the current position
// of size up to maxsize.  filters out repeat sequences which are redundant
// (e.g. TATA -> TA)
map<string, int> AlleleParser::repeatCounts(int maxsize) {
    map<string, int> counts;
    int position = currentPosition;
    int sequenceposition = position - currentSequenceStart;
    for (int i = 1; i <= maxsize; ++i) {
        // subseq here i bases
        string seq = currentSequence.substr(sequenceposition, i);
        // go left.
        int j = sequenceposition - i;
        int left = 0;
        while (j - i >= 0 && seq == currentSequence.substr(j, i)) {
            j -= i;
            ++left;
        }
        // go right.
        j = sequenceposition + i;
        int right = 0;
        while (j + i < currentSequence.size() && seq == currentSequence.substr(j, i)) {
            j += i;
            ++right;
        }
        // if we went left and right a non-zero number of times, 
        if (right > 0 || left > 0) {
            counts[seq] = right + left + 1;
        }
    }

    // filter out redundant repeat information
    if (counts.size() > 1) {
        map<string, int> filteredcounts;
        map<string, int>::iterator c = counts.begin();
        string prev = c->first;
        filteredcounts[prev] = c->second;  // shortest sequence
        ++c;
        for (; c != counts.end(); ++c) {
            int i = 0;
            string seq = c->first;
            while (i + prev.length() <= seq.length() && seq.substr(i, prev.length()) == prev) {
                i += prev.length();
            }
            if (i < seq.length()) {
                filteredcounts[seq] = c->second;
                prev = seq;
            }
        }
        return filteredcounts;
    } else {
        return counts;
    }
}

bool AlleleParser::hasInputVariantAllelesAtCurrentPosition(void) {
    map<long int, vector<Allele> >::iterator v = inputVariantAlleles.find(currentPosition);
    if (v != inputVariantAlleles.end()) {
        return true;
    } else {
        return false;
    }
}

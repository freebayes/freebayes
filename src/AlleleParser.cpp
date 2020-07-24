#include "AlleleParser.h"
#include "multichoose.h" // includes generic functions, so it must be included here
                         // otherwise we will get a linker error
                         // see: http://stackoverflow.com/questions/36039/templates-spread-across-multiple-files
                         // http://www.cplusplus.com/doc/tutorial/templates/ "Templates and Multi-file projects"
#include "multipermute.h"
#include "Logging.h"

using namespace std;

namespace {  // anonymous namespace

// Convert a std::string into an integer, ignoring any commas.
int stringToInt(string str) {
    str.erase(remove(str.begin(), str.end(), ','), str.end());
    return atoi(str.c_str());
}

}  // anonymous namespace


// open BAM input file
void AlleleParser::openBams(void) {

    // report differently if we have one or many bam files
    if (parameters.bams.size() == 1) {
        DEBUG("Opening BAM format alignment input file: " << parameters.bams.front() << " ...");
    } else if (parameters.bams.size() > 1) {
        DEBUG("Opening " << parameters.bams.size() << " BAM format alignment input files");
        for (vector<string>::const_iterator b = parameters.bams.begin();
                b != parameters.bams.end(); ++b) {
            DEBUG2(*b);
        }
    }

#ifdef HAVE_BAMTOOLS
    if (parameters.useStdin) {
        if (!bamMultiReader.Open(parameters.bams)) {
            ERROR("Could not read BAM data from stdin");
            cerr << bamMultiReader.GetErrorString() << endl;
            exit(1);
        }
    } else {
        if (!bamMultiReader.Open(parameters.bams)) {
            ERROR("Could not open input BAM files");
            cerr << bamMultiReader.GetErrorString() << endl;
            exit(1);
        } else {
            if (!bamMultiReader.LocateIndexes()) {
                ERROR("Opened BAM reader without index file, jumping is disabled.");
                cerr << bamMultiReader.GetErrorString() << endl;
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
        if (!bamMultiReader.SetExplicitMergeOrder(bamMultiReader.MergeByCoordinate)) {
            ERROR("could not set sort order to coordinate");
            cerr << bamMultiReader.GetErrorString() << endl;
            exit(1);
        }
    }

#else

    if (parameters.useStdin) {
        if (!bamMultiReader.Open("-")) {
            ERROR("Could not read BAM data from stdin");
            exit(1);
        }
    } else {
      for (std::vector<std::string>::const_iterator i = parameters.bams.begin(); i != parameters.bams.end(); ++i){
            string b = *i;
            // set reference to supplied fasta if the alignment file ends with cram
            if(b.substr(b.size()-4).compare("cram") == 0){
                DEBUG("Setting cram reference for bam reader")
                bamMultiReader.SetCramReference(parameters.fasta);
            }else{
                // reset the reference if this alignment file is no cram
                DEBUG("Unsetting cram reference for bam reader")
                bamMultiReader.SetCramReference("");
            }

            if (!bamMultiReader.Open(*i)) {
                ERROR("Could not open input BAM file: " + *i);
                exit(1);
            }else {
            /*if (!bamMultiReader.LocateIndexes()) {
                    ERROR("Opened BAM reader without index file, jumping is disabled.");
                    cerr << bamMultiReader.GetErrorString() << endl;
                    if (!targets.empty()) {
                        ERROR("Targets specified but no BAM index file provided.");
                        ERROR("FreeBayes cannot jump through targets in BAM files without BAM index files, exiting.");
                        ERROR("Please generate a BAM index file eithe, e.g.:");
                        ERROR("    \% bamtools index -in <bam_file>");
                        ERROR("    \% samtools index <bam_file>");
                        exit(1);
                    }
            }*/
            }
        }
        /*if (!bamMultiReader.SetExplicitMergeOrder(bamMultiReader.MergeByCoordinate)) {
            ERROR("could not set sort order to coordinate");
            cerr << bamMultiReader.GetErrorString() << endl;
            exit(1);
	    }*/

    }
#endif

    // from PR 319 below
#ifdef HAVE_BAMTOOLS
    if (!parameters.useStdin) {
        BamReader reader;
        for (vector<string>::const_iterator b = parameters.bams.begin();
             b != parameters.bams.end(); ++b) {
            reader.Open(*b);
            string bamHeader = reader.GetHeaderText();
            vector<string> headerLines = split(bamHeader, '\n');
            bamHeaderLines.insert(bamHeaderLines.end(), headerLines.begin(), headerLines.end());
            reader.Close();
        }
    } else {
        bamHeaderLines = split(bamMultiReader.GetHeaderText(), '\n');
    }
#else

    // retrieve header information
    string bamHeader = bamMultiReader.GETHEADERTEXT;
    bamHeaderLines = split(bamHeader, '\n');

#endif

    DEBUG(" done");
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
                size_t colpos = r->find(":");
                if (colpos != string::npos) {
                    string fieldname = r->substr(0, colpos);
                    if (fieldname == "PL") {
                        tech = r->substr(colpos+1);
                    } else if (fieldname == "ID") {
                        readGroupID = r->substr(colpos+1);
                    }
                }
            }

            if (tech.empty()) {
                if (!sequencingTechnologies.empty()) {
                    cerr << "no sequencing technology specified in @RG tag (no PL: in @RG tag) " << endl << headerLine << endl;
                }
            } else {
                map<string, string>::iterator s = readGroupToTechnology.find(readGroupID);
                if (s != readGroupToTechnology.end()) {
                    if (s->second != tech) {
                        ERROR("multiple technologies (PL) map to the same read group (RG)" << endl
                              << endl
                              << "technologies " << tech << " and " << s->second << " map to " << readGroupID << endl
                              << endl
                              << "As freebayes operates on a virtually merged stream of its input files," << endl
                              << "it will not be possible to determine what technology an alignment belongs to" << endl
                              << "at runtime." << endl
                              << endl
                              << "To resolve the issue, ensure that RG ids are unique to one technology" << endl
                              << "across all the input files to freebayes." << endl
                              << endl
                              << "See bamaddrg (https://github.com/ekg/bamaddrg) for a method which can" << endl
                              << "add RG tags to alignments." << endl);
                        exit(1);
                    }
                    // if it's the same technology and RG combo, no worries
                }
                readGroupToTechnology[readGroupID] = tech;
                technologies[tech] = true;
            }

            if (readGroupID.empty()) {
                cerr << "could not find ID: in @RG tag " << endl << headerLine << endl;
                continue;
            }
            //string name = nameParts.back();
            //mergedHeader.append(1, '\n');
            //cerr << "found read group id " << readGroupID << " containing sample " << name << endl;
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
                size_t colpos = r->find(":");
                if (colpos != string::npos) {
                    string fieldname = r->substr(0, colpos);
                    if (fieldname == "SM") {
                        name = r->substr(colpos+1);
                    } else if (fieldname == "ID") {
                        readGroupID = r->substr(colpos+1);
                    }
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
                    ERROR("multiple samples (SM) map to the same read group (RG)" << endl
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
    headerss
        << "##fileformat=VCFv4.2" << endl
        << "##fileDate=" << dateStr() << endl
        << "##source=freeBayes " << VERSION_GIT << endl
        << "##reference=" << reference.filename << endl;

    for (REFVEC::const_iterator it = referenceSequences.begin();
	 it != referenceSequences.end(); ++it)
      headerss << "##contig=<ID=" << it->REFNAME << ",length=" << it->REFLEN << ">" << endl;

    headerss
        << "##phasing=none" << endl
        << "##commandline=\"" << parameters.commandline << "\"" << endl
        << "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of samples with data\">" << endl
        << "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total read depth at the locus\">" << endl
        << "##INFO=<ID=DPB,Number=1,Type=Float,Description=\"Total read depth per bp at the locus; bases in reads overlapping / bases in haplotype\">" << endl

        // allele frequency metrics
        << "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Total number of alternate alleles in called genotypes\">" << endl
        << "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">" << endl
        << "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Estimated allele frequency in the range (0,1]\">" << endl

        // observation counts
        << "##INFO=<ID=RO,Number=1,Type=Integer,Description=\"Count of full observations of the reference haplotype.\">" << endl
        << "##INFO=<ID=AO,Number=A,Type=Integer,Description=\"Count of full observations of this alternate haplotype.\">" << endl
        << "##INFO=<ID=PRO,Number=1,Type=Float,Description=\"Reference allele observation count, with partial observations recorded fractionally\">" << endl
        << "##INFO=<ID=PAO,Number=A,Type=Float,Description=\"Alternate allele observations, with partial observations recorded fractionally\">" << endl

        // qualities
        << "##INFO=<ID=QR,Number=1,Type=Integer,Description=\"Reference allele quality sum in phred\">" << endl
        << "##INFO=<ID=QA,Number=A,Type=Integer,Description=\"Alternate allele quality sum in phred\">" << endl
        << "##INFO=<ID=PQR,Number=1,Type=Float,Description=\"Reference allele quality sum in phred for partial observations\">" << endl
        << "##INFO=<ID=PQA,Number=A,Type=Float,Description=\"Alternate allele quality sum in phred for partial observations\">" << endl


        // binomial balance metrics
        << "##INFO=<ID=SRF,Number=1,Type=Integer,Description=\"Number of reference observations on the forward strand\">" << endl
        << "##INFO=<ID=SRR,Number=1,Type=Integer,Description=\"Number of reference observations on the reverse strand\">" << endl
        << "##INFO=<ID=SAF,Number=A,Type=Integer,Description=\"Number of alternate observations on the forward strand\">" << endl
        << "##INFO=<ID=SAR,Number=A,Type=Integer,Description=\"Number of alternate observations on the reverse strand\">" << endl
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
        << "##INFO=<ID=RPL,Number=A,Type=Float,Description=\"Reads Placed Left: number of reads supporting the alternate balanced to the left (5') of the alternate allele\">" << endl
        //<< "##INFO=<ID=RPLR,Number=A,Type=Float,Description=\"Reads Placed Left: number of reads supporting the alternate balanced to the left (5') of the alternate allele\">" << endl
        << "##INFO=<ID=RPR,Number=A,Type=Float,Description=\"Reads Placed Right: number of reads supporting the alternate balanced to the right (3') of the alternate allele\">" << endl
        //<< "##INFO=<ID=RPRR,Number=A,Type=Float,Description=\"Reads Placed Right: number of reads supporting the alternate balanced to the right (3') of the alternate allele\">" << endl
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
        /*
        << "##INFO=<ID=XRM,Number=1,Type=Float,Description=\"Reference allele read mismatch rate: The rate of SNPs + MNPs + INDELs in reads supporting the reference allele.\">" << endl
        << "##INFO=<ID=XRS,Number=1,Type=Float,Description=\"Reference allele read SNP rate: The rate of per-base mismatches (SNPs + MNPs) in reads supporting the reference allele.\">" << endl
        << "##INFO=<ID=XRI,Number=1,Type=Float,Description=\"Reference allele read INDEL rate: The rate of INDELs (gaps) in reads supporting the reference allele.\">" << endl
        << "##INFO=<ID=XAM,Number=A,Type=Float,Description=\"Alternate allele read mismatch rate: The rate of SNPs + MNPs + INDELs in reads supporting the alternate allele, excluding the called variant.\">" << endl
        << "##INFO=<ID=XAS,Number=A,Type=Float,Description=\"Alternate allele read SNP rate: The rate of per-base mismatches (SNPs + MNPs) in reads supporting the alternate allele, excluding the called variant.\">" << endl
        << "##INFO=<ID=XAI,Number=A,Type=Float,Description=\"Alternate allele read INDEL rate: The rate of INDELs (gaps) in reads supporting the alternate allele, excluding the called variant.\">" << endl
        */
        // error rate ratios
        //<< "##INFO=<ID=ARM,Number=A,Type=Float,Description=\"Alternate allele / reference allele read mismatch ratio: The rate of SNPs + MNPs + INDELs in reads supporting the alternate allele versus reads supporting the reference allele, excluding the called variant.\">" << endl
        //<< "##INFO=<ID=ARS,Number=A,Type=Float,Description=\"Alternate allele / reference allele read SNP ratio: The rate of per-base mismatches (SNPs + MNPs) in reads supporting the alternate allele versus reads supporting the reference allele, excluding the called variant.\">" << endl
        //<< "##INFO=<ID=ARI,Number=A,Type=Float,Description=\"Alternate allele / reference allele read INDEL ratio: The ratio in rate rate of INDELs (gaps) in reads supporting the alternate allele versus reads supporting the reference allele, excluding the called variant.\">" << endl

        // supplementary information about the site
        << "##INFO=<ID=ODDS,Number=1,Type=Float,Description=\"The log odds ratio of the best genotype combination to the second-best.\">" << endl
        << "##INFO=<ID=GTI,Number=1,Type=Integer,Description=\"Number of genotyping iterations required to reach convergence or bailout.\">" << endl
        //<< "##INFO=<ID=TS,Number=0,Type=Flag,Description=\"site has transition SNP\">" << endl
        //<< "##INFO=<ID=TV,Number=0,Type=Flag,Description=\"site has transversion SNP\">" << endl
        //<< "##INFO=<ID=CpG,Number=0,Type=Flag,Description=\"CpG site (either CpG, TpG or CpA)\">" << endl
        << "##INFO=<ID=TYPE,Number=A,Type=String,Description=\"The type of allele, either snp, mnp, ins, del, or complex.\">" << endl
        << "##INFO=<ID=CIGAR,Number=A,Type=String,Description=\"The extended CIGAR representation of each alternate allele, with the exception that '=' is replaced by 'M' to ease VCF parsing.  Note that INDEL alleles do not have the first matched base (which is provided by default, per the spec) referred to by the CIGAR.\">" << endl
        //<< "##INFO=<ID=SNP,Number=0,Type=Flag,Description=\"SNP allele at site\">" << endl
        //<< "##INFO=<ID=MNP,Number=0,Type=Flag,Description=\"MNP allele at site\">" << endl
        //<< "##INFO=<ID=INS,Number=0,Type=Flag,Description=\"insertion allele at site\">" << endl
        //<< "##INFO=<ID=DEL,Number=0,Type=Flag,Description=\"deletion allele at site\">" << endl
        //<< "##INFO=<ID=COMPLEX,Number=0,Type=Flag,Description=\"complex allele (insertion/deletion/substitution composite) at site\">" << endl
        << "##INFO=<ID=NUMALT,Number=1,Type=Integer,Description=\"Number of unique non-reference alleles in called genotypes at this position.\">" << endl
        << "##INFO=<ID=MEANALT,Number=A,Type=Float,Description=\"Mean number of unique non-reference allele observations per sample with the corresponding alternate alleles.\">" << endl
        //<< "##INFO=<ID=HWE,Number=1,Type=Float,Description=\"Phred-scaled discrete HWE prior probability of the genotyping across all samples.\">" << endl
        << "##INFO=<ID=LEN,Number=A,Type=Integer,Description=\"allele length\">" << endl
        << "##INFO=<ID=MQM,Number=A,Type=Float,Description=\"Mean mapping quality of observed alternate alleles\">" << endl
        << "##INFO=<ID=MQMR,Number=1,Type=Float,Description=\"Mean mapping quality of observed reference alleles\">" << endl
        << "##INFO=<ID=PAIRED,Number=A,Type=Float,Description=\"Proportion of observed alternate alleles which are supported by properly paired read fragments\">" << endl
        << "##INFO=<ID=PAIREDR,Number=1,Type=Float,Description=\"Proportion of observed reference alleles which are supported by properly paired read fragments\">" << endl
        << "##INFO=<ID=MIN_DP,Number=1,Type=Integer,Description=\"Minimum depth in gVCF output block.\">" << endl
        << "##INFO=<ID=END,Number=1,Type=Integer,Description=\"Last position (inclusive) in gVCF output record.\">" << endl;

    // sequencing technology tags, which vary according to input data
    for (vector<string>::iterator st = sequencingTechnologies.begin(); st != sequencingTechnologies.end(); ++st) {
        string& tech = *st;
        headerss << "##INFO=<ID=technology." << tech << ",Number=A,Type=Float,Description=\"Fraction of observations supporting the alternate observed in reads from " << tech << "\">" << endl;
    }

    if (parameters.showReferenceRepeats) {
        headerss << "##INFO=<ID=REPEAT,Number=1,Type=String,Description=\"Description of the local repeat structures flanking the current position\">" << endl;
    }

    string gqType = "Float";
    if (parameters.strictVCF)
        gqType = "Integer";

        // format fields for genotypes
    headerss << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << endl
        << "##FORMAT=<ID=GQ,Number=1,Type=" << gqType << ",Description=\"Genotype Quality, the Phred-scaled marginal (or unconditional) probability of the called genotype\">" << endl
        // this can be regenerated with RA, AA, QR, QA
        << "##FORMAT=<ID=GL,Number=G,Type=Float,Description=\"Genotype Likelihood, log10-scaled likelihoods of the data given the called genotype for each possible genotype generated from the reference and alternate alleles given the sample ploidy\">" << endl
	//<< "##FORMAT=<ID=GLE,Number=1,Type=String,Description=\"Genotype Likelihood Explicit, same as GL, but with tags to indicate the specific genotype.  For instance, 0^-75.22|1^-223.42|0/0^-323.03|1/0^-99.29|1/1^-802.53 represents both haploid and diploid genotype likilehoods in a biallelic context\">" << endl
        << "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">" << endl
        << "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Number of observation for each allele\">" << endl
        << "##FORMAT=<ID=RO,Number=1,Type=Integer,Description=\"Reference allele observation count\">" << endl
        << "##FORMAT=<ID=QR,Number=1,Type=Integer,Description=\"Sum of quality of the reference observations\">" << endl
        << "##FORMAT=<ID=AO,Number=A,Type=Integer,Description=\"Alternate allele observation count\">" << endl
        << "##FORMAT=<ID=QA,Number=A,Type=Integer,Description=\"Sum of quality of the alternate observations\">" << endl
        << "##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description=\"Minimum depth in gVCF output block.\">" << endl
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
    // variant input for analysis and targeting
    if (!parameters.variantPriorsFile.empty()) {
        variantCallInputFile.open(parameters.variantPriorsFile);
        currentVariant = new vcflib::Variant(variantCallInputFile);
        usingVariantInputAlleles = true;

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

    // haplotype alleles for constructing haplotype alleles
    if (!parameters.haplotypeVariantFile.empty()) {
        haplotypeVariantInputFile.open(parameters.haplotypeVariantFile);
        usingHaplotypeBasisAlleles = true;
    }
}

void AlleleParser::loadBamReferenceSequenceNames(void) {

    //--------------------------------------------------------------------------
    // read reference sequences from input file
    //--------------------------------------------------------------------------

    // store the names of all the reference sequences in the BAM file
    referenceSequences = bamMultiReader.GETREFDATA;
    int i = 0;
    for (REFVEC::iterator r = referenceSequences.begin(); r != referenceSequences.end(); ++r) {
        referenceIDToName[i] = r->REFNAME;
        ++i;
    }

    DEBUG("Number of ref seqs: " << bamMultiReader.GETREFNUM);
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

bool AlleleParser::hasMoreInputVariants(void) {
    pair<int, long> next = nextInputVariantPosition();
    return next.first != -1;
}

bool AlleleParser::loadNextPositionWithAlignmentOrInputVariant(BAMALIGN& alignment) {

    pair<int, long> next = nextInputVariantPosition();
    if (next.first != -1) {
        int varRefID = next.first;
        if (!hasMoreAlignments || varRefID < alignment.REFID || (varRefID == alignment.REFID && next.second < alignment.POSITION)) {
	  return loadNextPositionWithInputVariant();
        } else {
	  loadReferenceSequence(alignment);
        }
    } else {
      loadReferenceSequence(alignment);
    }
    return true;
}

bool AlleleParser::loadNextPositionWithInputVariant(void) {
  pair<int, long> next = nextInputVariantPosition();
  if (next.first != -1) {
    //cerr << "Next is " << next.first << ":" << next.second << endl;
    loadReferenceSequence(referenceIDToName[next.first]);
    currentPosition = next.second;
    rightmostHaplotypeBasisAllelePosition = currentPosition;
    return true;
  } else {
    return false;
  }
}

// alignment-based method for loading the first bit of our reference sequence
void AlleleParser::loadReferenceSequence(BAMALIGN& alignment) {
  loadReferenceSequence(referenceIDToName[alignment.REFID]);
  currentPosition = alignment.POSITION;
}

void AlleleParser::loadReferenceSequence(string& seqname) {
    if (currentSequenceName != seqname) {
        currentSequenceName = seqname;
        currentSequenceStart = 0;
        currentRefID = bamMultiReader.GETREFID(currentSequenceName);
        currentSequence = uppercase(reference.getRawSequence(currentSequenceName));
        // check the first few characters and verify they are not garbage
        string validBases = "ACGTURYKMSWBDHVN-";
        size_t found = currentSequence.substr(0, 100).find_first_not_of(validBases);
        if (found != string::npos) {
            ERROR("Found non-DNA character " << currentSequence.at(found)
                  << " at position " << found << " in " << seqname << endl
                  << "Is your reference compressed or corrupted? "
                  << "freebayes requires an uncompressed reference sequence.");
            exit(1);
        }
        currentSequence = reference.getSequence(currentSequenceName);
    }
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

        size_t foundLastColon = region.rfind(":");

        // we only have a single string, use the whole sequence as the target
        if (foundLastColon == string::npos) {
            startSeq = region;
            startPos = 0;
            stopPos = -1;
        } else {
            startSeq = region.substr(0, foundLastColon);
            string sep = "..";
            size_t foundRangeSep = region.find(sep, foundLastColon);
            if (foundRangeSep == string::npos) {
                sep = "-";
                foundRangeSep = region.find(sep, foundLastColon);
            }
            if (foundRangeSep == string::npos) {
                startPos = stringToInt(region.substr(foundLastColon + 1));
                // differ from bamtools in this regard, in that we process only
                // the specified position if a range isn't given
                stopPos = startPos + 1;
            } else {
                startPos = stringToInt(region.substr(foundLastColon + 1, foundRangeSep - foundLastColon).c_str());
                // if we have range sep specified, but no second number, read to the end of sequence
                if (foundRangeSep + sep.size() != region.size()) {
                    stopPos = stringToInt(region.substr(foundRangeSep + sep.size()).c_str()); // end-exclusive, bed-format
                } else {
                    stopPos = -1;
                }
            }
        }

        //DEBUG("startPos == " << startPos);
        //DEBUG("stopPos == " << stopPos);

        // REAL BED format is 0 based, half open (end base not included)
        BedTarget bd(startSeq,
                    (startPos == 0) ? 0 : startPos,
                    ((stopPos == -1) ? reference.sequenceLength(startSeq) : stopPos) - 1); // internally, we use 0-base inclusive end
        DEBUG("will process reference sequence " << startSeq << ":" << bd.left << ".." << bd.right + 1);
        targets.push_back(bd);
        bedReader.targets.push_back(bd);

    }

    // check validity of targets wrt. reference
    for (vector<BedTarget>::iterator e = targets.begin(); e != targets.end(); ++e) {
        BedTarget& bd = *e;
        // internally, we use 0-base inclusive end
        if (bd.left < 0 || bd.right + 1 > reference.sequenceLength(bd.seq)) {
            ERROR("Target region coordinates (" << bd.seq << " "
                    << bd.left << " " << bd.right + 1
                    << ") outside of reference sequence bounds ("
                    << bd.seq << " " << reference.sequenceLength(bd.seq) << ") terminating.");
            exit(1);
        }
        if (bd.right < bd.left) {
            ERROR("Invalid target region coordinates (" << bd.seq << " " << bd.left << " " << bd.right + 1 << ")"
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
    REFVEC::iterator refIter = referenceSequences.begin();
    REFVEC::iterator refEnd  = referenceSequences.end();
    for( ; refIter != refEnd; ++refIter) {
        REFDATA refData = *refIter;
        string refName = refData.REFNAME;
        BedTarget bd(refName, 0, refData.REFLEN); // 0-based inclusive internally
        DEBUG2("will process reference sequence " << refName << ":" << bd.left << ".." << bd.right + 1);
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
        for (REFVEC::iterator r = referenceSequences.begin(); r != referenceSequences.end(); ++r) {
            sampleCNV.setPloidy(referenceSampleName, r->REFNAME, 0, r->REFLEN, 1);
        }
    }

}

int AlleleParser::currentSamplePloidy(string const& sample) {
    return sampleCNV.ploidy(sample, currentSequenceName, currentPosition);
}

int AlleleParser::copiesOfLocus(Samples& samples) {
    int copies = 0;
    for (Samples::iterator s = samples.begin(); s != samples.end(); ++s) {
        string const& name = s->first;
        copies += currentSamplePloidy(name);
    }
    return copies;
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
        // expects 0-based, fully-closed, and we're only checking a single
        // base, so start == end.
        if (bedReader.targetsOverlap(currentSequenceName, currentPosition, currentPosition)) {
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
    currentPosition = 0;
    currentTarget = NULL; // to be initialized on first call to getNextAlleles
    currentReferenceAllele = NULL; // same, NULL is brazenly used as an initialization flag
    justSwitchedTargets = false;  // flag to trigger cleanup of Allele*'s and objects after jumping targets
    hasMoreAlignments = true; // flag to track when we run out of alignments in the current target or BAM files
    currentSequenceStart = 0;
    lastHaplotypeLength = 0;
    usingHaplotypeBasisAlleles = false;
    usingVariantInputAlleles = false;
    rightmostHaplotypeBasisAllelePosition = 0;
    rightmostInputAllelePosition = 0;
    nullSample = new Sample();
    referenceSampleName = "reference_sample";

    // initialization
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
int AlleleParser::currentSequencePosition(const BAMALIGN& alignment) {
  return alignment.POSITION - currentSequenceStart;
}

// relative current position within the cached currentSequence
int AlleleParser::currentSequencePosition() {
    return currentPosition - currentSequenceStart;
}

char AlleleParser::currentReferenceBaseChar(void) {
    return toupper(*currentReferenceBaseIterator());
}

string AlleleParser::currentReferenceBaseString(void) {
    return currentSequence.substr(floor(currentPosition) - currentSequenceStart, 1);
}

string::iterator AlleleParser::currentReferenceBaseIterator(void) {
    return currentSequence.begin() + (floor(currentPosition) - currentSequenceStart);
}

string AlleleParser::currentReferenceHaplotype(void) {
    return currentSequence.substr(floor(currentPosition) - currentSequenceStart, lastHaplotypeLength);
}

string AlleleParser::referenceSubstr(long int pos, unsigned int len) {
    return uppercase(reference.getSubSequence(currentSequenceName, floor(pos), len));
}

bool AlleleParser::isCpG(string& altbase) {
    // bounds check
    if (floor(currentPosition) - currentSequenceStart - 1 < 0
            || floor(currentPosition) - currentSequenceStart + 1 >= currentSequence.size()) {
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

void capBaseQuality(BAMALIGN& alignment, int baseQualityCap) {
  string rQual = alignment.QUALITIES;
  char qualcap = qualityInt2Char(baseQualityCap);
  for (string::iterator c = rQual.begin(); c != rQual.end(); ++c) {
    if (qualityChar2ShortInt(*c) > baseQualityCap) {
      *c = qualcap;
    }
  }
}

void RegisteredAlignment::addAllele(Allele newAllele, bool mergeComplex, int maxComplexGap, bool boundIndels) {

    if (newAllele.alternateSequence.size() != newAllele.baseQualities.size()) {
        cerr << "new allele qualities not == in length to sequence: " << newAllele << endl;
        assert(false);
    }

    //cerr << "adding allele " << newAllele << " to " << alleles.size() << " alleles" << endl;

    alleleTypes |= newAllele.type;
    alleles.push_back(newAllele);
}

void RegisteredAlignment::clumpAlleles(bool mergeComplex, int maxComplexGap, bool boundIndels) {
    // remove any empty alleles, and skip if we go totally empty
    alleles.erase(remove_if(alleles.begin(), alleles.end(), isEmptyAllele), alleles.end());
    if (!alleles.size()) return;
    vector<bool> toMerge(alleles.size());
    if (maxComplexGap >= 0) {
        for (int i = 1; i < alleles.size()-1; ++i) {
            const Allele& lastAllele = alleles[i-1];
            const Allele& currAllele = alleles[i];
            const Allele& nextAllele = alleles[i+1];
            if (lastAllele.isNull() || currAllele.isNull() || nextAllele.isNull()) continue;
            if (!lastAllele.isReference() && !nextAllele.isReference()
                && (currAllele.isReference() && currAllele.referenceLength <= maxComplexGap
                    || !currAllele.isReference())) {
                toMerge[i-1] = true;
                toMerge[i] = true;
                toMerge[i+1] = true;
            } else if (!lastAllele.isReference() && !currAllele.isReference() && !nextAllele.isReference()) {
                toMerge[i-1] = true;
                toMerge[i] = true;
                toMerge[i+1] = true;
            } else if (!lastAllele.isReference() && !currAllele.isReference()) {
                toMerge[i-1] = true;
                toMerge[i] = true;
            } else if (!nextAllele.isReference() && !currAllele.isReference()) {
                toMerge[i] = true;
                toMerge[i+1] = true;
            }
        }
    }
    // find clumps by combining all reference alleles <= maxComplexGap bases with their neighbors
    vector<Allele> newAlleles;
    for (int i = 0; i < toMerge.size(); ++i) {
        bool merge = toMerge[i];
        if (merge) {
            Allele merged = alleles[i];
            while (++i < toMerge.size() && toMerge[i]) {
                merged.mergeAllele(alleles[i], ALLELE_COMPLEX);
            }
            newAlleles.push_back(merged);
            if (i < toMerge.size() && !toMerge[i]) { // we broke on this allele
                newAlleles.push_back(alleles[i]);
            }
        } else {
            newAlleles.push_back(alleles[i]);
        }
    }
    alleles = newAlleles;
    newAlleles.clear();
    alleles.erase(remove_if(alleles.begin(), alleles.end(), isEmptyAllele), alleles.end());
    // maintain flanking bases
    for (int i = 1; i < alleles.size()-1; ++i) {
        Allele& lastAllele = alleles[i-1];
        Allele& currAllele = alleles[i];
        Allele& nextAllele = alleles[i+1];
        if (!lastAllele.length || !currAllele.length || !nextAllele.length) continue;
        vector<pair<int, string> > lastCigar = splitCigar(lastAllele.cigar);
        vector<pair<int, string> > currCigar = splitCigar(currAllele.cigar);
        vector<pair<int, string> > nextCigar = splitCigar(nextAllele.cigar);
        string currFirstOp = currCigar.front().second;
        string currLastOp = currCigar.back().second;
        string lastLastOp = lastCigar.back().second;
        string nextFirstOp = nextCigar.front().second;
        if ((currFirstOp == "I" || currFirstOp == "D")
            && (lastLastOp == "M" || lastLastOp == "X")) {
            // split from the last onto curr
            string seq; vector<pair<int, string> > cig; vector<short> quals;
            lastAllele.subtractFromEnd(1, seq, cig, quals);
            currAllele.addToStart(seq, cig, quals);
        }
        if ((currLastOp == "I" || currLastOp == "D")
            && (nextFirstOp == "M" || nextFirstOp == "X")) {
            string seq; vector<pair<int, string> > cig; vector<short> quals;
            nextAllele.subtractFromStart(1, seq, cig, quals);
            currAllele.addToEnd(seq, cig, quals);
            // split from the next onto curr
            if (nextAllele.length == 0) i+=3;// skip past this allele if we've axed it
        }
    }
    // now we want to remove any null alleles, etc.
    // also indels at the start and end of reads
    alleles.erase(remove_if(alleles.begin(), alleles.end(), isEmptyAllele), alleles.end());
    //cerr << "merged to alleles " << alleles << endl;
}

// TODO erase alleles which are beyond N bp before the current position on position step
void AlleleParser::updateHaplotypeBasisAlleles(long int pos, int referenceLength) {
    if (pos + referenceLength > rightmostHaplotypeBasisAllelePosition) {
        stringstream r;
        //r << currentSequenceName << ":" << rightmostHaplotypeBasisAllelePosition << "-" << pos + referenceLength + CACHED_BASIS_HAPLOTYPE_WINDOW;
        //cerr << "getting variants in " << r.str() << endl;

        // tabix expects 1-based, fully closed regions for ti_parse_region()
        // (which is what setRegion() calls eventually)
        if (haplotypeVariantInputFile.setRegion(currentSequenceName,
                                                rightmostHaplotypeBasisAllelePosition + 1,
                                                pos + referenceLength + CACHED_BASIS_HAPLOTYPE_WINDOW + 1)) {
            //cerr << "the vcf line " << haplotypeVariantInputFile.line << endl;
            // get the variants in the target region
            vcflib::Variant var(haplotypeVariantInputFile);
            while (haplotypeVariantInputFile.getNextVariant(var)) {
                //cerr << "input variant: " << var << endl;

                // the following stanza is for parsed
                // alternates. instead use whole haplotype calls, as
                // alternates can be parsed prior to providing the
                // file as input.
                /*
                  for (vector<string>::iterator a = var.alt.begin(); a != var.alt.end(); ++a) {
                  haplotypeBasisAlleles[var.position].insert(AllelicPrimitive(var.ref.size(), *a));
                  }
                */

                map<string, vector<vcflib::VariantAllele> > variants = var.parsedAlternates();
                for (map<string, vector<vcflib::VariantAllele> >::iterator a = variants.begin(); a != variants.end(); ++a) {
                    for (vector<vcflib::VariantAllele>::iterator v = a->second.begin(); v != a->second.end(); ++v) {
                        //cerr << v->ref << "/" << v->alt << endl;
                        if (v->ref != v->alt) {
                            //cerr << "basis allele " << v->position << " " << v->ref << "/" << v->alt << endl;
                            haplotypeBasisAlleles[v->position].push_back(AllelicPrimitive(v->ref, v->alt));
                            //cerr << "number of alleles at position " <<  haplotypeBasisAlleles[v->position].size() << endl;
                        }
                    }
                }

            }
        } else {
            // indicates empty region
            //ERROR("Could not set haplotype-basis VCF file to target region");
            //exit(1);
        }
        // set the rightmost haplotype position to trigger the next update
        rightmostHaplotypeBasisAllelePosition = pos + referenceLength + CACHED_BASIS_HAPLOTYPE_WINDOW;
    }
}


bool AlleleParser::allowedHaplotypeBasisAllele(long int pos, string& ref, string& alt) {
    // check the haplotypeBasisAllele map for membership of the allele in question in the current sequence
    //cerr << "is allowed: " << pos << " " << ref << "/" << alt << " ?" << endl;
    if (!usingHaplotypeBasisAlleles) {
        return true; // always true if we aren't using the haplotype basis allele system
    } else {
        map<long int, vector<AllelicPrimitive> >::iterator p = haplotypeBasisAlleles.find(pos);
        if (p != haplotypeBasisAlleles.end()) {
            vector<AllelicPrimitive>& alleles = p->second;
            for (vector<AllelicPrimitive>::iterator z = alleles.begin(); z != alleles.end(); ++z) {
                //cerr << "overlapping allele " << z->ref << ":" << z->alt << endl;
                if (z->ref == ref && z->alt == alt) {
                    //cerr << "yess" << endl;
                    return true;
                }
            }
        }
        return false;
    }

}

Allele AlleleParser::makeAllele(RegisteredAlignment& ra,
                                AlleleType type,
                                long int pos,
                                int length,
                                int basesLeft,
                                int basesRight,
                                string& readSequence,
                                string& sampleName,
                                BAMALIGN& alignment,
                                string& sequencingTech,
                                long double qual,
                                string& qualstr
    ) {

    string cigar;
    int reflen = length;

    if (type == ALLELE_REFERENCE) {
        cigar = convert(length) + "M";
    } else if (type == ALLELE_SNP || type == ALLELE_MNP) {
        cigar = convert(length) + "X";
    } else if (type == ALLELE_INSERTION) {
        reflen = 0;
        cigar = convert(length) + "I";
    } else if (type == ALLELE_DELETION) {
        cigar = convert(length) + "D";
    } else if (type == ALLELE_NULL) {
        cigar = convert(length) + "N";
    }

    string refSequence;
    if (type != ALLELE_NULL) { // only used for non null allele, avoid soft clipping edge cases
        refSequence = currentSequence.substr(pos - currentSequenceStart, reflen);
    }

    long int repeatRightBoundary = pos;

    // check if it's allowed

    // if it isn't allowed
    // and referenceLength > 0, make a reference allele with reference quality
    // if referenceLength == 0 (insertion), make a reference allele with 0 length (it will be filtered out in another context)

    // if it is allowed, make a normal allele

    // if not, adjust the allele so that it's a reference allele with preset BQ and length
    // in effect, this means creating a reference allele of the reference length of the allele with 0 BQ

    // NB, if we are using haplotype basis alleles the algorithm forces
    // alleles that aren't in the haplotype basis set into the reference space
    if (type != ALLELE_REFERENCE
        && type != ALLELE_NULL
        && !allowedHaplotypeBasisAllele(pos + 1,
                                        refSequence,
                                        readSequence)) {
        type = ALLELE_REFERENCE;
        length = referenceLengthFromCigar(cigar);
        cigar = convert(length) + "M";
        // by adjusting the cigar, we implicitly adjust
        // allele.referenceLength, which is calculated when the allele is made
        qualstr = string(length, qualityInt2Char(0));
        readSequence = currentSequence.substr(pos - currentSequenceStart, length);
    }

    // cache information about repeat structure in the alleles, to
    // allow haplotype construction to be forced to extend across
    // tandem repeats and homopolymers when indels are present
    if (type == ALLELE_INSERTION || type == ALLELE_DELETION) {
        string alleleseq;
        if (type == ALLELE_INSERTION) {
            alleleseq = readSequence;
        } else if (type == ALLELE_DELETION) {
            alleleseq = refSequence;
        }
        map<long int, map<string, int> >::iterator rc = cachedRepeatCounts.find(pos);
        if (rc == cachedRepeatCounts.end()) {
            cachedRepeatCounts[pos] = repeatCounts(pos - currentSequenceStart, currentSequence, 12);
            rc = cachedRepeatCounts.find(pos);
        }
        map<string, int>& matchedRepeatCounts = rc->second;
        for (map<string, int>::iterator r = matchedRepeatCounts.begin(); r != matchedRepeatCounts.end(); ++r) {
            const string& repeatunit = r->first;
            int rptcount = r->second;
            string repeatstr = repeatunit * rptcount;
            // assumption of left-alignment may be problematic... so this should be updated
            if (repeatstr.size() >= parameters.minRepeatSize && isRepeatUnit(alleleseq, repeatunit)) {
                // determine the boundaries of the repeat
                long int p = pos - currentSequenceStart;
                // adjust to ensure we hit the first of the repeatstr
                size_t startpos = currentSequence.find(repeatstr, max((long int) 0, p - (long int) repeatstr.size() - 1));
                long int leftbound = startpos + currentSequenceStart;
                if (startpos == string::npos) {
                    cerr << "could not find repeat sequence?" << endl;
                    cerr << "repeat sequence: " << repeatstr << endl;
                    cerr << "currentsequence start: " << currentSequenceStart << endl;
                    cerr << currentSequence << endl;
                    cerr << "matched repeats:" << endl;
                    for (map<string, int>::iterator q = matchedRepeatCounts.begin(); q != matchedRepeatCounts.end(); ++q) {
                        cerr << q->first  << " : " << q->second << endl;
                        cerr << "... at position " << pos << endl;
                    }
                    break; // ignore right-repeat boundary in this case
                }
                repeatRightBoundary = leftbound + repeatstr.size() + 1; // 1 past edge of repeat
            }
        }

        // a dangerous game
        int start = pos - currentSequenceStart;
        double minEntropy = parameters.minRepeatEntropy;
        while (minEntropy > 0 && // ignore if turned off
               // don't run off the end of the current sequence
               repeatRightBoundary - currentSequenceStart < currentSequence.size() &&
               // there is no point in going past the alignment end
               // because we won't make a haplotype call unless we have a covering observation from a read
               repeatRightBoundary < alignment.ENDPOSITION &&
               entropy(currentSequence.substr(start, repeatRightBoundary - pos)) < minEntropy) {
            ++repeatRightBoundary;
        }

        // edge case, the indel is an insertion and matches the reference to the right
        // this means there is a repeat structure in the read, but not the ref
        if (currentSequence.substr(pos - currentSequenceStart, length) == readSequence) {
            repeatRightBoundary = max(repeatRightBoundary, pos + length + 1);
        }
    }

    string qnamer = alignment.QNAME;

    return Allele(type,
                  currentSequenceName,
                  pos,
                  &currentPosition,
                  &currentReferenceBase,
                  length,
                  repeatRightBoundary,
                  basesLeft,
                  basesRight,
                  readSequence,
                  sampleName,
                  qnamer,
                  ra.readgroup,
                  sequencingTech,
                  !alignment.ISREVERSESTRAND,
                  max(qual, (long double) 0), // ensure qual is at least 0
                  qualstr,
                  alignment.MAPPINGQUALITY,
                  alignment.ISPAIRED,
                  alignment.ISMATEMAPPED,
                  alignment.ISPROPERPAIR,
                  cigar,
                  &ra.alleles,
                  alignment.POSITION,
                  alignment.ENDPOSITION);

}

RegisteredAlignment& AlleleParser::registerAlignment(BAMALIGN& alignment, RegisteredAlignment& ra, string& sampleName, string& sequencingTech) {

    string rDna = alignment.QUERYBASES;
    string rQual = alignment.QUALITIES;
    if (qualityChar2LongDouble(rQual.at(0)) == -1) {
        // force rQual to be 0
        char q0 = qualityInt2Char(0);
        for (size_t i = 0; i < rQual.size(); ++i) {
            rQual[i] = q0;
        }
    }
    int rp = 0;  // read position, 0-based relative to read
    int csp = currentSequencePosition(alignment); // current sequence position, 0-based relative to currentSequence
    int sp = alignment.POSITION;  // sequence position
    if (usingHaplotypeBasisAlleles) {
        updateHaplotypeBasisAlleles(sp, alignment.ALIGNEDBASES);
    }

#ifdef VERBOSE_DEBUG
    if (parameters.debug2) {
        DEBUG2("registering alignment " << rp << " " << csp << " " << sp << endl <<
               "alignment readName " << alignment.QNAME << endl <<
               "alignment isPaired " << alignment.ISPAIRED << endl <<
               "alignment isMateMapped " << alignment.ISMATEMAPPED << endl <<
               "alignment isProperPair " << alignment.ISPROPERPAIR << endl <<
               "alignment mapQual " << alignment.MAPPINGQUALITY << endl <<
               "alignment sampleID " << sampleName << endl <<
               "alignment position " << alignment.POSITION << endl <<
               "alignment length " << alignment.ALIGNMENTLENGTH << endl <<
               "alignment AlignedBases.size() " << alignment.ALIGNEDBASES << endl <<
               "alignment GetEndPosition() " << alignment.ENDPOSITION << endl <<
               "alignment end position " << alignment.POSITION + alignment.ALIGNEDBASES);

        stringstream cigarss;
        int alignedLength = 0;
	CIGAR cig = alignment.GETCIGAR;
        for (CIGAR::const_iterator c = cig.begin(); c != cig.end(); ++c) {
            cigarss << c->CIGTYPE << c->CIGLEN;
            if (c->CIGTYPE == 'D')
                alignedLength += c->CIGLEN;
            if (c->CIGTYPE == 'M')
                alignedLength += c->CIGLEN;
        }

        DEBUG2("alignment cigar " << cigarss.str());

        DEBUG2("current sequence pointer: " << csp);

        DEBUG2("read:          " << rDna);
        DEBUG2("aligned bases: " << alignment.QUERYBASES);
        DEBUG2("qualities:     " << alignment.QUALITIES);
        DEBUG2("reference seq: " << currentSequence.substr(csp, alignment.ALIGNEDBASES));
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

    /*    std::cerr << "********" << std::endl
	      << alignment.QueryBases << std::endl
	      << alignment.AlignedBases << std::endl;
    vector<CigarOp>::const_iterator cigarIter2 = alignment.CigarData.begin();
    vector<CigarOp>::const_iterator cigarEnd2  = alignment.CigarData.end();
    for (; cigarIter2 != cigarEnd2; ++cigarIter2)
      std::cerr << cigarIter2->Length << cigarIter2->Type;
    std::cerr << std::endl;
    */

    CIGAR cigar = alignment.GETCIGAR;
    CIGAR::const_iterator cigarIter = cigar.begin();
    CIGAR::const_iterator cigarEnd  = cigar.end();
    for ( ; cigarIter != cigarEnd; ++cigarIter ) {
        int l = cigarIter->CIGLEN;
        char t = cigarIter->CIGTYPE;
	  DEBUG2("cigar item: " << t << l);

        if (t == 'M' || t == 'X' || t == '=') { // match or mismatch
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
                         << alignment.QNAME << endl
                         << currentSequenceName << ":" << (long unsigned int) currentPosition + 1 << endl
		      //<< alignment.AlignedBases << endl
                         << currentSequence.substr(csp, alignment.ALIGNEDBASES) << endl;
		    cerr << " RP " << rp << " " << rDna  <<" len " << rDna.length() <<  std::endl;
                    abort();
                }

                // convert base quality value into short int
                long double qual = qualityChar2LongDouble(rQual.at(rp));

                // get reference allele
                string sb;
                try {
                    sb = currentSequence.at(csp);
                } catch (std::out_of_range outOfRange) {
                    cerr << "Exception: Alignment reports a match past the end of the current reference sequence." << endl
                         << "This suggests alignment corruption or a mismatch between this reference and the alignments." << endl
                         << "Are you sure that you are calling against the same reference you aligned to?" << endl
                         << "Sequence: " << currentSequenceName << ":" << (long unsigned int) currentPosition + 1 << endl
                         << "Alignment: " << alignment.QNAME << " @ " << alignment.POSITION << "-" << alignment.ENDPOSITION << endl;
                    break;
                }

                // record mismatch if we have a mismatch here
                if (b != sb || sb == "N") {  // when the reference is N, we should always call a mismatch
                    if (firstMatch < csp) {
                        int length = csp - firstMatch;
                        string readSequence = rDna.substr(rp - length, length);
                        string qualstr = rQual.substr(rp - length, length);
                        // record 'reference' allele for last matching region
                        if (allATGC(readSequence)) {
                            ra.addAllele(
                                makeAllele(ra,
                                           ALLELE_REFERENCE,
                                           sp - length,
                                           length,
                                           rp, // bases left (for first base in ref allele)
                                           alignment.SEQLEN - rp, // bases right (for first base in ref allele)
                                           readSequence,
                                           sampleName,
                                           alignment,
                                           sequencingTech,
                                           alignment.MAPPINGQUALITY, // reference allele quality == mapquality
                                           qualstr),
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
                    string readSequence = rDna.substr(rp - length, length);
                    string qualstr = rQual.substr(rp - length, length);
                    for (int j = 0; j < length; ++j) {
                        long double lqual = qualityChar2LongDouble(qualstr.at(j));
                        string qualp = qualstr.substr(j, 1);
                        string rs = readSequence.substr(j, 1);
                        if (allATGC(rs)) {
                            ra.addAllele(
                                makeAllele(ra,
                                           ALLELE_SNP,
                                           sp - length + j,
                                           1,
                                           rp - length - j, // bases left
                                           alignment.SEQLEN - rp + j, // bases right
                                           rs,
                                           sampleName,
                                           alignment,
                                           sequencingTech,
                                           lqual,
                                           qualp),
                                parameters.allowComplex, parameters.maxComplexGap);

                        } else {
                            ra.addAllele(
                                makeAllele(ra,
                                           ALLELE_NULL,
                                           sp - length + j,
                                           1,
                                           rp - length - j, // bases left
                                           alignment.SEQLEN - rp + j, // bases right
                                           rs,
                                           sampleName,
                                           alignment,
                                           sequencingTech,
                                           lqual,
                                           qualp),
                                parameters.allowComplex, parameters.maxComplexGap);
                        }
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
                string readSequence = rDna.substr(rp - length, length);
                string qualstr = rQual.substr(rp - length, length);
                for (int j = 0; j < length; ++j) {
                    long double lqual = qualityChar2LongDouble(qualstr.at(j));
                    string qualp = qualstr.substr(j, 1);
                    string rs = readSequence.substr(j, 1);
                    if (allATGC(rs)) {
                        ra.addAllele(
                            makeAllele(ra,
                                       ALLELE_SNP,
                                       sp - length + j,
                                       1,
                                       rp - length - j, // bases left
                                       alignment.SEQLEN - rp + j, // bases right
                                       rs,
                                       sampleName,
                                       alignment,
                                       sequencingTech,
                                       lqual,
                                       qualp),
                            parameters.allowComplex, parameters.maxComplexGap);

                    } else {
                        ra.addAllele(
                            makeAllele(ra,
                                       ALLELE_NULL,
                                       sp - length + j,
                                       1,
                                       rp - length - j, // bases left
                                       alignment.SEQLEN - rp + j, // bases right
                                       rs,
                                       sampleName,
                                       alignment,
                                       sequencingTech,
                                       lqual,
                                       qualp),
                            parameters.allowComplex, parameters.maxComplexGap);
                    }
                }
                // or, if we are not in a mismatch, construct the last reference allele of the match
            } else if (firstMatch < csp) {
                int length = csp - firstMatch;
                //string matchingSequence = currentSequence.substr(csp - length, length);
                string readSequence = rDna.substr(rp - length, length);
                string qualstr = rQual.substr(rp - length, length);
                if (allATGC(readSequence)) {
                    ra.addAllele(
                        makeAllele(ra,
                                   ALLELE_REFERENCE,
                                   sp - length,
                                   length,
                                   rp, // bases left (for first base in ref allele)
                                   alignment.SEQLEN - rp, // bases right (for first base in ref allele)
                                   readSequence,
                                   sampleName,
                                   alignment,
                                   sequencingTech,
                                   alignment.MAPPINGQUALITY, // ... hmm
                                   qualstr),
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

            long double qual;
            if (parameters.useMinIndelQuality) {
                qual = minQuality(qualstr);
                //qual = averageQuality(qualstr);
            } else {
                // quality, scaled inversely by the ratio between the quality
                // string length and the length of the event
                qual = sumQuality(qualstr);
                // quality adjustment:
                // scale the quality by the inverse harmonic sum of the length of
                // the quality string X a scaling constant derived from the ratio
                // between the length of the quality string and the length of the
                // allele
                //qual += ln2phred(log((long double) l / (long double) L));
                qual += ln2phred(log((long double) L / (long double) l));
                qual /= harmonicSum(l);
            }

            string refseq = currentSequence.substr(csp, l);
            // some aligners like to report deletions at the beginnings and ends of reads.
            // without any sequence in the read to support this, it is hard to believe
            // that these deletions are real, so we ignore them here.
	      CIGAR cigar = alignment.GETCIGAR;
	      if (cigarIter != cigar.begin()      // guard against deletion at beginning
		  && (cigarIter+1) != cigar.end() // and against deletion at end
                && allATGC(refseq)) {
                string nullstr;
                ra.addAllele(
                    makeAllele(ra,
                               ALLELE_DELETION,
                               sp,
                               l,
                               rp, // bases left (for first base in ref allele)
                               alignment.SEQLEN - rp, // bases right (for first base in ref allele)
                               nullstr, // no read sequence for deletions
                               sampleName,
                               alignment,
                               sequencingTech,
                               qual,
                               nullstr), // no qualstr for deletions
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

            long double qual;
            if (parameters.useMinIndelQuality) {
                qual = minQuality(qualstr);
                //qual = averageQuality(qualstr); // does not work as well as the min
            } else {
                // quality, scaled inversely by the ratio between the quality
                // string length and the length of the event
                qual = sumQuality(qualstr);
                // quality adjustment:
                // scale the quality by the inverse harmonic sum of the length of
                // the quality string X a scaling constant derived from the ratio
                // between the length of the quality string and the length of the
                // allele
                //qual += ln2phred(log((long double) l / (long double) L));
                qual += ln2phred(log((long double) L / (long double) l));
                qual /= harmonicSum(l);
            }

            string readseq = rDna.substr(rp, l);
            if (allATGC(readseq)) {
                string qualstr = rQual.substr(rp, l);
                ra.addAllele(
                    makeAllele(ra,
                               ALLELE_INSERTION,
                               sp,
                               l,
                               rp - l, // bases left (for first base in ref allele)
                               alignment.SEQLEN - rp, // bases right (for first base in ref allele)
                               readseq,
                               sampleName,
                               alignment,
                               sequencingTech,
                               qual,
                               qualstr),
                    parameters.allowComplex, parameters.maxComplexGap);
            }
            ++ra.indelCount;

            rp += l;

            // handle other cigar element types
        } else if (t == 'S') { // soft clip, clipped sequence present in the read not matching the reference
            if (sp - l < 0) {
                // nothing to do, soft clip is beyond the beginning of the reference
            } else {
                string qualstr = rQual.substr(rp, l);
                string readseq = alignment.QUERYBASES.substr(rp, l);
                // skip these bases in the read
                ra.addAllele(
                    makeAllele(ra,
                               ALLELE_NULL,
                               sp,
                               l,
                               rp, // bases left (for first base in ref allele)
                               alignment.SEQLEN - rp, // bases right
                               readseq,
                               sampleName,
                               alignment,
                               sequencingTech,
                               alignment.MAPPINGQUALITY,
                               qualstr),
                    parameters.allowComplex, parameters.maxComplexGap);
            }
            rp += l;// sp += l; csp += l;
        } else if (t == 'H') { // hard clip on the read, clipped sequence is not present in the read
            // the alignment position is the first non-clipped base.
            // thus, hard clipping seems to just be an indicator that we clipped something
            // here we do nothing
            //sp += l; csp += l;
        } else if (t == 'N') { // skipped region in the reference not present in read, aka splice
            // skip these bases in the read
            // the following block could be enabled to process them, if they are desired
            /*
            string nullstr;
            ra.addAllele(
                makeAllele(ra,
                           ALLELE_NULL,
                           sp - l,
                           l,
                           rp - l, // bases left
                           alignment.SEQLEN - rp, // bases right
                           nullstr,
                           sampleName,
                           alignment,
                           sequencingTech,
                           alignment.MAPPINGQUALITY,
                           nullstr),
                parameters.allowComplex, parameters.maxComplexGap);
            */
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

    // ignore insertions, deletions, and N's which occur at the end of the read with
    // no reference-matching bases before the end of the read
    /*
    if (parameters.boundIndels &&
        (ra.alleles.back().isInsertion()
         || ra.alleles.back().isDeletion()
         || ra.alleles.back().isNull())) {
        ra.alleles.pop_back();
    }
    */
    ra.clumpAlleles(parameters.allowComplex, parameters.maxComplexGap, parameters.boundIndels);

    DEBUG2("alleles:" << endl << join(ra.alleles, "\n"));

    /*
      cerr << "ra.alleles.size() = " << ra.alleles.size() << endl;
      for (vector<Allele>::iterator a = ra.alleles.begin(); a != ra.alleles.end(); ++a) {
      cerr << *a << endl;
      }
    */

    return ra;

}


void AlleleParser::updateAlignmentQueue(long int position,
                                        vector<Allele*>& newAlleles,
                                        bool gettingPartials) {

    DEBUG2("updating alignment queue");
    DEBUG2("currentPosition = " << position
           << "; currentSequenceStart = " << currentSequenceStart
           << "; currentSequence end = " << currentSequence.size() + currentSequenceStart);

    // make sure we have sequence for the *first* alignment
    //extendReferenceSequence(currentAlignment);

    // push to the front until we get to an alignment that doesn't overlap our
    // current position or we reach the end of available alignments
    // filter input reads; only allow mapped reads with a certain quality
    DEBUG2("currentAlignment.Position == " << currentAlignment.POSITION
           << ", currentAlignment.AlignedBases.size() == " << currentAlignment.ALIGNEDBASES
           << ", currentPosition == " << position
           << ", currentSequenceStart == " << currentSequenceStart
           << " .. + currentSequence.size() == " << currentSequenceStart + currentSequence.size()
        );

    if (hasMoreAlignments
        && currentAlignment.POSITION <= position
        && currentAlignment.REFID == currentRefID) {
        do {
            DEBUG2("top of alignment parsing loop");
            DEBUG("alignment: " << currentAlignment.QNAME);
            // get read group, and map back to a sample name
            string readGroup;
#ifdef HAVE_BAMTOOLS
            if (!currentAlignment.GetTag("RG", readGroup)) {
#else
	      currentAlignment.GetZTag("RG", readGroup);
	    if (readGroup.empty()) {
#endif
                if (!oneSampleAnalysis) {
                    ERROR("Couldn't find read group id (@RG tag) for BAM Alignment " <<
                          currentAlignment.QNAME << " at " << currentSequenceName << ":"
                          << position + 1 << " EXITING!");
                    exit(1);
                } else {
                    readGroup = "unknown";
                }
            } else {
                if (oneSampleAnalysis) {
                    ERROR("No read groups specified in BAM header, but alignment " <<
                          currentAlignment.QNAME << " at " << currentSequenceName << ":"
                          << position + 1 << " has a read group.");
                    exit(1);
                }
            }

            // skip this alignment if we are not analyzing the sample it is drawn from
            if (readGroupToSampleNames.find(readGroup) == readGroupToSampleNames.end()) {
                ERROR("could not find sample matching read group id " << readGroup);
                continue;
            }

            // skip this alignment if we are not using duplicate reads (we remove them by default)
            if (currentAlignment.ISDUPLICATE && !parameters.useDuplicateReads) {
                DEBUG("skipping alignment " << currentAlignment.QNAME << " because it is a duplicate read");
                continue;
            }

            // skip unmapped alignments, as they cannot be used in the algorithm
            if (!currentAlignment.ISMAPPED) {
                DEBUG("skipping alignment " << currentAlignment.QNAME << " because it is not mapped");
                continue;
            }

            // skip alignments which have no aligned bases
            if (currentAlignment.ALIGNEDBASES == 0) {
                DEBUG("skipping alignment " << currentAlignment.QNAME << " because it has no aligned bases");
                continue;
            }

            // skip alignments which are non-primary
            if (currentAlignment.SecondaryFlag()) {
                DEBUG("skipping alignment " << currentAlignment.QNAME << " because it is not marked primary");
                continue;
            }

            if (!gettingPartials && currentAlignment.ENDPOSITION < position) {
                cerr << currentAlignment.QNAME << " at " << currentSequenceName << ":" << currentAlignment.POSITION << " is out of order!"
                     << " expected after " << position << endl;
                continue;
            }

            // otherwise, get the sample name and register the alignment to generate a sequence of alleles
            // we have to register the alignment to acquire some information required by filters
            // such as mismatches

            // initially skip reads with low mapping quality (what happens if MapQuality is not in the file)
            if (currentAlignment.MAPPINGQUALITY >= parameters.MQL0) {
                // extend our cached reference sequence to allow processing of this alignment
                //extendReferenceSequence(currentAlignment);
                // left realign indels
                if (parameters.leftAlignIndels) {
                    int length = currentAlignment.ENDPOSITION - currentAlignment.POSITION + 1;
                    stablyLeftAlign(currentAlignment,
                                    currentSequence.substr(currentSequencePosition(currentAlignment), length));
                }
                // get sample name
                string sampleName = readGroupToSampleNames[readGroup];
                string sequencingTech;
                map<string, string>::iterator t = readGroupToTechnology.find(readGroup);
                if (t != readGroupToTechnology.end()) {
                    sequencingTech = t->second;
                }
                // limit base quality if cap set
                if (parameters.baseQualityCap != 0) {
                    capBaseQuality(currentAlignment, parameters.baseQualityCap);
                }
                // do we exceed coverage anywhere?
                // do we touch anything where we had exceeded coverage?
                // if so skip this read, and mark and remove processed alignments and registered alleles overlapping the coverage capped position
                bool considerAlignment = true;
                if (parameters.skipCoverage > 0) {
                    for (unsigned long int i =  currentAlignment.POSITION; i < currentAlignment.ENDPOSITION; ++i) {
                        unsigned long int x = ++coverage[i];
                        if (x > parameters.skipCoverage && !gettingPartials) {
                            considerAlignment = false;
                            // we're exceeding coverage at this position for the first time, so clean up
                            if (!coverageSkippedPositions.count(i)) {
                                // clean up reads overlapping this position
                                removeCoverageSkippedAlleles(registeredAlleles, i);
                                removeCoverageSkippedAlleles(newAlleles, i);
                                // remove the alignments overlapping this position
                                removeRegisteredAlignmentsOverlappingPosition(i);
                                // record that the position is capped
                                coverageSkippedPositions.insert(i);
                            }
                        }
                    }
                }
                // decomposes alignment into a set of alleles
                // here we get the deque of alignments ending at this alignment's end position
                deque<RegisteredAlignment>& rq = registeredAlignments[currentAlignment.ENDPOSITION];
                //cerr << "parameters capcoverage " << parameters.capCoverage << " " << rq.size() << endl;
                if (considerAlignment) {
                    // and insert the registered alignment into that deque
                    rq.push_front(RegisteredAlignment(currentAlignment));
                    RegisteredAlignment& ra = rq.front();
                    registerAlignment(currentAlignment, ra, sampleName, sequencingTech);
                    // backtracking if we have too many mismatches
                    // or if there are no recorded alleles
                    if (ra.alleles.empty()
                        || ((float) ra.mismatches / (float) currentAlignment.SEQLEN) > parameters.readMaxMismatchFraction
                        || ra.mismatches > parameters.RMU
                        || ra.snpCount > parameters.readSnpLimit
                        || ra.indelCount > parameters.readIndelLimit) {
                        rq.pop_front(); // backtrack
                    } else {
                        // push the alleles into our new alleles vector
                        for (vector<Allele>::iterator allele = ra.alleles.begin(); allele != ra.alleles.end(); ++allele) {
                            newAlleles.push_back(&*allele);
                        }
                    }
                }
            }
	    } while ((hasMoreAlignments = GETNEXT(bamMultiReader, currentAlignment))
                 && currentAlignment.POSITION <= position
                 && currentAlignment.REFID == currentRefID);
    }

    DEBUG2("... finished pushing new alignments");

}

void AlleleParser::removeRegisteredAlignmentsOverlappingPosition(long unsigned int pos) {
    map<long unsigned int, deque<RegisteredAlignment> >::iterator f = registeredAlignments.begin();
    map<long unsigned int, set<deque<RegisteredAlignment>::iterator> > alignmentsToErase;
    set<Allele*> allelesToErase;
    while (f != registeredAlignments.end()) {
        for (deque<RegisteredAlignment>::iterator d = f->second.begin(); d != f->second.end(); ++d) {
            if (d->start <= pos && d->end > pos) {
                alignmentsToErase[f->first].insert(d);
                for (vector<Allele>::iterator a = d->alleles.begin(); a != d->alleles.end(); ++a) {
                    allelesToErase.insert(&*a);
                }
            }
        }
        ++f;
    }
    // clean up registered alleles--- maybe this should be done externally?
    for (vector<Allele*>::iterator a = registeredAlleles.begin(); a != registeredAlleles.end(); ++a) {
        if (allelesToErase.count(*a)) {
            *a = NULL;
        }
    }
    registeredAlleles.erase(remove(registeredAlleles.begin(), registeredAlleles.end(), (Allele*)NULL), registeredAlleles.end());
    if (alignmentsToErase.size()) {
        for (map<long unsigned int, set<deque<RegisteredAlignment>::iterator> >::iterator e = alignmentsToErase.begin();
             e != alignmentsToErase.end(); ++e) {
            deque<RegisteredAlignment> updated;
            map<long unsigned int, deque<RegisteredAlignment> >::iterator f = registeredAlignments.find(e->first);
            assert(f != registeredAlignments.end());
            for (deque<RegisteredAlignment>::iterator d = f->second.begin(); d != f->second.end(); ++d) {
                if (!e->second.count(d)) {
                    updated.push_back(*d);
                }
            }
            f->second = updated;
        }
    }
}

void AlleleParser::addToRegisteredAlleles(vector<Allele*>& alleles) {
    registeredAlleles.insert(registeredAlleles.end(),
                             alleles.begin(),
                             alleles.end());
}

// updates registered alleles and erases the unused portion of our cached reference sequence
void AlleleParser::updateRegisteredAlleles(void) {

    // remove reference alleles which are no longer overlapping the current position
    vector<Allele*>& alleles = registeredAlleles;

    for (vector<Allele*>::iterator allele = alleles.begin(); allele != alleles.end(); ++allele) {
        long unsigned int position = (*allele)->position;
        if (position + (*allele)->referenceLength < currentPosition) {
            *allele = NULL;
        }
    }

    alleles.erase(remove(alleles.begin(), alleles.end(), (Allele*)NULL), alleles.end());

}

pair<int, long int> AlleleParser::nextInputVariantPosition(void) {
    // are we past the last one in the sequence?
    if (usingVariantInputAlleles &&
        ((inputVariantAlleles.find(currentRefID) != inputVariantAlleles.end()
          && inputVariantAlleles[currentRefID].upper_bound(currentPosition) != inputVariantAlleles[currentRefID].end())
         || inputVariantAlleles.upper_bound(currentRefID) != inputVariantAlleles.end())) {
        map<long, vector<Allele> >& inChrom = inputVariantAlleles[currentRefID];
        map<long, vector<Allele> >::iterator ic = inChrom.upper_bound(currentPosition);
        if (ic != inChrom.end()) {
            return make_pair(currentRefID, ic->first);
        } else {
            // find next chrom with input alleles
            map<int, map<long, vector<Allele> > >::iterator nc
              = inputVariantAlleles.upper_bound(currentRefID);
            if (nc != inputVariantAlleles.end()) {
                return make_pair(nc->first, nc->second.begin()->first);
            } else {
                return make_pair(-1, 0);
            }
        }
    }
    return make_pair(-1, 0);
}

void AlleleParser::getAllInputVariants(void) {
    string nullstr;
    getInputVariantsInRegion(nullstr);
}

void AlleleParser::getInputVariantsInRegion(string& seq, long start, long end) {

    if (!usingVariantInputAlleles) return;

    // get the variants in the target region
    vcflib::Variant var(variantCallInputFile);
    if (!seq.empty()) {
        variantCallInputFile.setRegion(seq, start, end);
    }
    bool ok;
    while ((ok = variantCallInputFile.getNextVariant(*currentVariant))) {

        long int pos = currentVariant->position - 1;
        // get alternate alleles
        bool includePreviousBaseForIndels = true;
        map<string, vector<vcflib::VariantAllele> > variantAlleles = currentVariant->parsedAlternates();
        // TODO this would be a nice option: why does it not work?
        //map<string, vector<vcflib::VariantAllele> > variantAlleles = currentVariant->flatAlternates();
        vector< vector<vcflib::VariantAllele> > orderedVariantAlleles;
        for (vector<string>::iterator a = currentVariant->alt.begin();
          a != currentVariant->alt.end(); ++a) {
            orderedVariantAlleles.push_back(variantAlleles[*a]);
        }

        vector<Allele> genotypeAlleles;
        set<long int> alternatePositions;

        for (vector< vector<vcflib::VariantAllele> >::iterator
          g = orderedVariantAlleles.begin();
          g != orderedVariantAlleles.end(); ++g) {

            vector<vcflib::VariantAllele>& altAllele = *g;

            vector<Allele> alleles;

            for (vector<vcflib::VariantAllele>::iterator v = altAllele.begin();
              v != altAllele.end(); ++v) {
                vcflib::VariantAllele& variant = *v;
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
                    type = ALLELE_DELETION;
                    len = variant.ref.size() - variant.alt.size();
                    allelePos -= 1;
                    reflen = len + 2;
                    alleleSequence =
                        reference.getSubSequence(currentVariant->sequenceName, allelePos, 1)
                        + alleleSequence
                        + reference.getSubSequence(currentVariant->sequenceName, allelePos+1+len, 1);
                    cigar = "1M" + convert(len) + "D" + "1M";
                } else {
                    // we always include the flanking bases for these elsewhere, so here too in order to be consistent and trigger use
                    type = ALLELE_INSERTION;
                    // add previous base and post base to match format typically used for calling
                    allelePos -= 1;
                    alleleSequence =
                        reference.getSubSequence(currentVariant->sequenceName, allelePos, 1)
                        + alleleSequence
                        + reference.getSubSequence(currentVariant->sequenceName, allelePos+1, 1);
                    len = variant.alt.size() - var.ref.size();
                    cigar = "1M" + convert(len) + "I" + "1M";
                    reflen = 2;
                }
                // TODO deal woth complex subs

                Allele allele = genotypeAllele(type, alleleSequence, (unsigned int) len, cigar, (unsigned int) reflen, allelePos);
                DEBUG("input allele: " << allele.referenceName << " " << allele);
                //cerr << "input allele: " << allele.referenceName << " " << allele << endl;

                //alleles.push_back(allele);
                genotypeAlleles.push_back(allele);

                if (allele.type != ALLELE_REFERENCE) {
                    inputVariantAlleles[bamMultiReader.GETREFID(currentVariant->sequenceName)][allele.position].push_back(allele);
                    alternatePositions.insert(allele.position);
                }
            }
        }
    }
}

void AlleleParser::updateInputVariants(long int pos, int referenceLength) {

    //cerr << "updating input variants (?) " << pos << " + " << referenceLength << " >? " << rightmostInputAllelePosition << endl;
    if (!usingVariantInputAlleles) return;

    if (pos + referenceLength > rightmostInputAllelePosition) {
        long int start = rightmostInputAllelePosition;
        if (start == 0) {
            start = rightmostHaplotypeBasisAllelePosition;
        }

        /*
        stringstream r;
        r << currentSequenceName << ":" << start
          << "-" << pos + referenceLength + CACHED_BASIS_HAPLOTYPE_WINDOW;
        cerr << "getting variants in " << r.str() << endl;
        */

        // tabix expects 1-based, fully closed regions for ti_parse_region()
        // (which is what setRegion() calls eventually)
        bool gotRegion = false;
        if (referenceLength > 0) {
            gotRegion = variantCallInputFile.setRegion(currentSequenceName,
                                                       start + 1,
                                                       pos + referenceLength + CACHED_BASIS_HAPLOTYPE_WINDOW + 1);
        } else {
            // whole chromosome
            gotRegion = variantCallInputFile.setRegion(currentSequenceName);
        }

        if (gotRegion) {

            // get the variants in the target region
            vcflib::Variant var(variantCallInputFile);
            bool ok;
            while ((ok = variantCallInputFile.getNextVariant(*currentVariant))) {

                DEBUG("getting input alleles from input VCF at position " << currentVariant->sequenceName << ":" << currentVariant->position);
                long int pos = currentVariant->position - 1;
                // get alternate alleles
                bool includePreviousBaseForIndels = true;
                map<string, vector<vcflib::VariantAllele> > variantAlleles = currentVariant->parsedAlternates();
                // TODO this would be a nice option: why does it not work?
                //map<string, vector<vcflib::VariantAllele> > variantAlleles = currentVariant->flatAlternates();
                vector< vector<vcflib::VariantAllele> > orderedVariantAlleles;
                for (vector<string>::iterator a = currentVariant->alt.begin(); a != currentVariant->alt.end(); ++a) {
                    orderedVariantAlleles.push_back(variantAlleles[*a]);
                }

                vector<Allele> genotypeAlleles;
                set<long int> alternatePositions;

                for (vector< vector<vcflib::VariantAllele> >::iterator g = orderedVariantAlleles.begin(); g != orderedVariantAlleles.end(); ++g) {

                    vector<vcflib::VariantAllele>& altAllele = *g;

                    vector<Allele> alleles;

                    for (vector<vcflib::VariantAllele>::iterator v = altAllele.begin(); v != altAllele.end(); ++v) {
                        vcflib::VariantAllele& variant = *v;
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
                            type = ALLELE_DELETION;
                            len = variant.ref.size() - variant.alt.size();
                            allelePos -= 1;
                            reflen = len + 2;
                            alleleSequence =
                                reference.getSubSequence(currentSequenceName, allelePos, 1)
                                + alleleSequence
                                + reference.getSubSequence(currentSequenceName, allelePos+1+len, 1);
                            cigar = "1M" + convert(len) + "D" + "1M";
                        } else {
                            // we always include the flanking bases for these elsewhere, so here too in order to be consistent and trigger use
                            type = ALLELE_INSERTION;
                            // add previous base and post base to match format typically used for calling
                            allelePos -= 1;
                            alleleSequence =
                                reference.getSubSequence(currentSequenceName, allelePos, 1)
                                + alleleSequence
                                + reference.getSubSequence(currentSequenceName, allelePos+1, 1);
                            len = variant.alt.size() - var.ref.size();
                            cigar = "1M" + convert(len) + "I" + "1M";
                            reflen = 2;
                        }
                        // TODO deal woth complex subs

                        Allele allele = genotypeAllele(type, alleleSequence, (unsigned int) len, cigar, (unsigned int) reflen, allelePos);
                        DEBUG("input allele: " << allele.referenceName << " " << allele);

                        //alleles.push_back(allele);
                        genotypeAlleles.push_back(allele);

                        if (allele.type != ALLELE_REFERENCE) {
                            inputVariantAlleles[bamMultiReader.GETREFID(allele.referenceName)][allele.position].push_back(allele);
                            alternatePositions.insert(allele.position);
                        }

                    }

                }

                // store the allele counts, if they are provided
                //
            }

            if (!ok) hasMoreVariants = false;
        }
        /*
        for (map<long int, vector<Allele> >::iterator v = inputVariantAlleles.begin(); v != inputVariantAlleles.end(); ++v) {
            vector<Allele>& iv = v->second;
            cerr << "input variants pos = " << v->first << endl;
            for (vector<Allele>::iterator a = iv.begin(); a != iv.end(); ++a) {
                cerr << *a << endl;
            }
        }
        */
        //rightmostHaplotypeBasisAllelePosition = pos + referenceLength + CACHED_BASIS_HAPLOTYPE_WINDOW;
        //rightmostInputAllelePosition = pos + referenceLength + CACHED_BASIS_HAPLOTYPE_WINDOW;
    }

}

/*
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
*/

void AlleleParser::removeAllelesWithoutReadSpan(vector<Allele*>& alleles, int probeLength, int haplotypeLength) {
    for (vector<Allele*>::iterator a = alleles.begin(); a != alleles.end(); ++a) {
        Allele* allele = *a;
        if (!(allele->position == currentPosition && allele->referenceLength == haplotypeLength))
            continue;
        // require additionally
        int additionalRequiredBases = probeLength - allele->alternateSequence.size();
        int requiredFlank = ceil((double) additionalRequiredBases / 2);
        DEBUG2(allele << " needs at least " << additionalRequiredBases
              << " bpleft " << allele->read5pNonNullBases() << " bpright " << allele->read3pNonNullBases());
        if (additionalRequiredBases > 0 &&
            (allele->read5pNonNullBases() < additionalRequiredBases
             || allele->read3pNonNullBases() < additionalRequiredBases)) {
            DEBUG("removing " << allele << " as it does not have the required probe length");
            *a = NULL;
        }
    }
    alleles.erase(remove(alleles.begin(), alleles.end(), (Allele*)NULL), alleles.end());
}

void AlleleParser::removeNonOverlappingAlleles(vector<Allele*>& alleles, int haplotypeLength, bool getAllAllelesInHaplotype) {
    for (vector<Allele*>::iterator a = alleles.begin(); a != alleles.end(); ++a) {
        Allele* allele = *a;
        if (allele->type == ALLELE_REFERENCE) {
            // does the reference allele overlap the haplotype
            if (getAllAllelesInHaplotype
                && !(currentPosition <= allele->position && allele->position < currentPosition + haplotypeLength)) {
                //cerr << *a << " is not in haplotype" << endl;
                *a = NULL;
            } else if (!(allele->position <= currentPosition
                         && allele->position + allele->referenceLength >= currentPosition + haplotypeLength)) {
                //cerr << *a << " is not fully overlapping haplotype from " << currentPosition << " to " << currentPosition + haplotypeLength << endl;
                *a = NULL;
            } else if (currentPosition < allele->position) { // not there yet
                //cerr << *a << " is not before current position" << endl;
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
    }
    alleles.erase(remove(alleles.begin(), alleles.end(), (Allele*)NULL), alleles.end());
}

// removes alleles which are filtered at the current position, and unsets their 'processed' flag so they are later evaluated
void AlleleParser::removeFilteredAlleles(vector<Allele*>& alleles) {
    for (vector<Allele*>::iterator allele = alleles.begin(); allele != alleles.end(); ++allele) {
        if ((*allele)->quality < parameters.BQL0 || (*allele)->currentBase == "N") {
            (*allele)->processed = false; // force re-processing later
            *allele = NULL;
        }
    }
    alleles.erase(remove(alleles.begin(), alleles.end(), (Allele*)NULL), alleles.end());
}

void AlleleParser::removePreviousAlleles(vector<Allele*>& alleles, long int position) {
    for (vector<Allele*>::iterator a = alleles.begin(); a != alleles.end(); ++a) {
        Allele* allele = *a;
        if (*a != NULL && allele->position + allele->referenceLength < position) {
            allele->processed = true;
            *a = NULL;
        }
    }
    alleles.erase(remove(alleles.begin(), alleles.end(), (Allele*)NULL), alleles.end());
}

void AlleleParser::removeCoverageSkippedAlleles(vector<Allele*>& alleles, long int position) {
    for (vector<Allele*>::iterator a = alleles.begin(); a != alleles.end(); ++a) {
        Allele* allele = *a;
        if (*a != NULL && allele->alignmentStart <= position && allele->alignmentEnd > position) {
            allele->processed = true;
            *a = NULL;
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

    DEBUG("to next target");

    clearRegisteredAlignments();
    coverageSkippedPositions.clear();
    cachedRepeatCounts.clear();
    coverage.clear();

    // reset haplotype length; there is no last call in this sequence; it isn't relevant
    lastHaplotypeLength = 0;

    if (targets.empty() && usingVariantInputAlleles) {
        // we are processing everything, so load the entire input variant allele set
        getAllInputVariants();
    }

    // load first target if we have targets and have not loaded the first
    if (!parameters.useStdin && !targets.empty()) {

        bool ok = false;

        // try to load the first target if we need to
        if (!currentTarget) {
            ok = loadTarget(&targets.front()) && getFirstAlignment();
        }

        // step through targets until we get to one with alignments
        while (!ok && currentTarget != &targets.back()) {
            if (!loadTarget(++currentTarget)) {
                continue;
            }
            if ((ok = getFirstAlignment())) {
                break;
            }
        }

        if (!ok) {
            return loadNextPositionWithInputVariant();
        }

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
        loadNextPositionWithAlignmentOrInputVariant(currentAlignment);
        //loadReferenceSequence(currentAlignment); // this seeds us with new reference sequence
        // however, if we have a target list of variants and we should also respect them
    // we've reached the end of file, or stdin
    } else if (parameters.useStdin || targets.empty()) {
        return false;
    }

    if (currentTarget && usingVariantInputAlleles) {
        getInputVariantsInRegion(currentTarget->seq, currentTarget->left, currentTarget->right);
    }

    loadReferenceSequence(currentSequenceName);

    justSwitchedTargets = true;
    return true;

}


// TODO refactor this to allow reading from stdin or reading the whole file
// without loading each sequence as a target
bool AlleleParser::loadTarget(BedTarget* target) {

    currentTarget = target;

    DEBUG("processing target " << currentTarget->desc << " " <<
          currentTarget->seq << " " << currentTarget->left << " " <<
          currentTarget->right + 1);
    DEBUG2("loading target reference subsequence");

    loadReferenceSequence(currentTarget->seq);

    DEBUG2("setting new position " << currentTarget->left);
    currentPosition = currentTarget->left;
    rightmostHaplotypeBasisAllelePosition = currentTarget->left;

#ifdef HAVE_BAMTOOLS
    if (!bamMultiReader.SetRegion(currentRefID, currentTarget->left, currentRefID, currentTarget->right + 1)) { // bamtools expects 0-based, half-open
        ERROR("Could not SetRegion to " << currentTarget->seq << ":" << currentTarget->left << ".." << currentTarget->right + 1);
        cerr << bamMultiReader.GetErrorString() << endl;
        return false;
    }
#else
    if (!bamMultiReader.SetRegion(SeqLib::GenomicRegion(currentRefID, currentTarget->left, currentTarget->right + 1))) { // bamtools expects 0-based, half-open
        ERROR("Could not SetRegion to " << currentTarget->seq << ":" << currentTarget->left << ".." << currentTarget->right + 1);
        return false;
    }
#endif

    if (variantCallInputFile.is_open()) {
        stringstream r;
        // tabix expects 1-based, fully closed regions for ti_parse_region()
        // (which is what setRegion() calls eventually)
        r << currentTarget->seq << ":" << currentTarget->left + 1 << "-" << currentTarget->right + 1;
        if (!variantCallInputFile.setRegion(r.str())) {
            WARNING("Could not set the region of the variants input file to " <<
                    currentTarget->seq << ":" << currentTarget->left << ".." <<
                    currentTarget->right + 1);
            //return false;
        } else {
            DEBUG("set region of variant input file to " <<
                    currentTarget->seq << ":" << currentTarget->left << ".." <<
                    currentTarget->right + 1);
        }
    }

    // now that we've jumped, reset the hasMoreAlignments counter
    hasMoreAlignments = true;

    DEBUG2("set region");

    return true;

}

bool AlleleParser::getFirstAlignment(void) {

    bool hasAlignments = true;
    if (!GETNEXT(bamMultiReader, currentAlignment)) {
      hasAlignments = false;
    } else {
      while (!currentAlignment.ISMAPPED) {
	if (!GETNEXT(bamMultiReader, currentAlignment)) {
	  hasAlignments = false;
	  break;
	}
      }
    }

    if (hasAlignments) {
        DEBUG2("got first alignment in target region");
    } else {
        if (currentTarget) {
            DEBUG("Could not find any mapped reads in target region " << currentSequenceName << ":" << currentTarget->left << ".." << currentTarget->right + 1);
        } else {
            DEBUG("Could not find any mapped reads in target region " << currentSequenceName);
        }
        return false;
    }

    return true;

}

bool AlleleParser::getFirstVariant(void) {

    hasMoreVariants = false;
    if (variantCallInputFile.is_open()) {
        if (!variantCallInputFile.getNextVariant(*currentVariant)) {
            hasMoreVariants = false;
        } else {
            hasMoreVariants = true;
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
    // is this our first position? (indicated by empty currentSequenceName)
    // if so, load it up
    bool first_pos = false;
    if (currentSequenceName.empty()) {
        DEBUG("loading first target");
        if (!toNextTarget()) {
            return false;
        }
        first_pos = true;
    }

    // here we assume we are processing an entire BAM or one contiguous region
    if (parameters.useStdin || targets.empty()) {
        // here we loop over unaligned reads at the beginning of a target
        // we need to get to a mapped read to figure out where we are
        while (hasMoreAlignments && !currentAlignment.ISMAPPED) {
            hasMoreAlignments = GETNEXT(bamMultiReader, currentAlignment);
        }
        // determine if we have more alignments or not
        if (!hasMoreAlignments) {
            if (hasMoreInputVariants()) {
                // continue as we have more variants
                DEBUG("continuing because we have more input variants");
                loadNextPositionWithInputVariant();
            } else if (registeredAlignments.empty()) {
                DEBUG("no more alignments in input");
                return false;
            } else if (currentPosition >= currentSequence.size() + currentSequenceStart) {
                DEBUG("no more alignments in input");
                DEBUG("at end of sequence");
                return false;
            } else {
                ++currentPosition;
            }
        } else {
            // step the position
            if (!first_pos) {
                ++currentPosition;
            }
            // if the current position of this alignment is outside of the reference sequence length
            // we need to switch references
            if (currentPosition >= reference.sequenceLength(currentSequenceName)
                || (registeredAlignments.empty() && currentRefID != currentAlignment.REFID)) {
                DEBUG("at end of sequence");
                clearRegisteredAlignments();
                coverageSkippedPositions.clear();
                cachedRepeatCounts.clear();
                coverage.clear();
                loadNextPositionWithAlignmentOrInputVariant(currentAlignment);
                justSwitchedTargets = true;
            }
        }
    } else {
        // or if it's not we should step to the next position
        if (!first_pos) {
            ++currentPosition;
        }
        // if we've run off the right edge of a target, jump
        if (currentPosition > currentTarget->right) {
            // time to move to a new target
            DEBUG("next position " << (long int) currentPosition
                  <<  " outside of current target right bound " << currentTarget->right + 1);
            // try to get to the next one, and if this fails, bail out
            if (!toNextTarget()) {
                DEBUG("no more targets, finishing");
                return false;
            }
            justSwitchedTargets = true;
        }
    }

    // so we have to make sure it's still there (this matters in low-coverage)
    currentReferenceBase = currentReferenceBaseChar();

    // handle the case in which we don't have targets but in which we've switched reference sequence

    DEBUG("processing position " << (long unsigned int) currentPosition + 1 << " in sequence " << currentSequenceName);
    vector<Allele*> newAlleles;
    updateAlignmentQueue(currentPosition, newAlleles);
    addToRegisteredAlleles(newAlleles);
    DEBUG2("updating variants");
    // done typically at each new read, but this handles the case where there is no data for a while
    //updateInputVariants(currentPosition, 1);

    // remove past registered alleles
    DEBUG2("marking previous alleles as processed and removing from registered alleles");
    removePreviousAlleles(registeredAlleles, currentPosition);

    // if we have alignments which ended at the previous base, erase them and their alleles
    DEBUG2("erasing old registered alignments");
    map<long unsigned int, deque<RegisteredAlignment> >::iterator f = registeredAlignments.begin();
    set<long unsigned int> positionsToErase;
    set<Allele*> allelesToErase;
    while (f != registeredAlignments.end()
           && f->first < currentPosition - lastHaplotypeLength) {
        for (deque<RegisteredAlignment>::iterator d = f->second.begin(); d != f->second.end(); ++d) {
            for (vector<Allele>::iterator a = d->alleles.begin(); a != d->alleles.end(); ++a) {
                allelesToErase.insert(&*a);
            }
        }
        positionsToErase.insert(f->first);
        ++f;
    }
    for (vector<Allele*>::iterator a = registeredAlleles.begin(); a != registeredAlleles.end(); ++a) {
        if (allelesToErase.count(*a)) {
            *a = NULL;
        }
    }
    registeredAlleles.erase(remove(registeredAlleles.begin(), registeredAlleles.end(), (Allele*)NULL), registeredAlleles.end());
    for (set<long unsigned int>::iterator p = positionsToErase.begin(); p != positionsToErase.end(); ++p) {
        registeredAlignments.erase(*p);
    }

    // and do the same for the variants from the input VCF
    DEBUG2("erasing old input variant alleles");
    int refid = bamMultiReader.GETREFID(currentSequenceName);
    if (inputVariantAlleles.find(refid) != inputVariantAlleles.end()) {
        map<long int, vector<Allele> >::iterator v = inputVariantAlleles[refid].begin();
        while (v != inputVariantAlleles[refid].end() && v->first < currentPosition) {
            inputVariantAlleles[refid].erase(v++);
        }
        for (map<int, map<long int, vector<Allele> > >::iterator v = inputVariantAlleles.begin();
             v != inputVariantAlleles.end(); ++v) {
            if (v->first != refid) inputVariantAlleles.erase(v);
        }
    }

    DEBUG2("erasing old input haplotype basis alleles");
    map<long int, vector<AllelicPrimitive> >::iterator z = haplotypeBasisAlleles.begin();
    while (z != haplotypeBasisAlleles.end() && z->first < currentPosition) {
        haplotypeBasisAlleles.erase(z++);
    }

    DEBUG2("erasing old cached repeat counts");
    map<long int, map<string, int> >::iterator rc = cachedRepeatCounts.begin();
    while (rc != cachedRepeatCounts.end() && rc->first < currentPosition) {
        cachedRepeatCounts.erase(rc++);
    }

    DEBUG2("erasing old coverage cap");
    while (coverageSkippedPositions.size() && *coverageSkippedPositions.begin() < currentPosition) {
        coverageSkippedPositions.erase(coverageSkippedPositions.begin());
    }

    DEBUG2("erasing old coverage counts");
    map<long unsigned int, long unsigned int>::iterator cov = coverage.begin();
    while (cov != coverage.end() && cov->first < currentPosition) {
        coverage.erase(cov++);
    }

    return true;

}

// XXX for testing only, steps targets but does nothing
bool AlleleParser::dummyProcessNextTarget(void) {

    if (!toNextTarget()) {
        DEBUG("no more targets, finishing");
        return false;
    }

    while (GETNEXT(bamMultiReader, currentAlignment)) { }

    return true;
}

void AlleleParser::removeDuplicateAlleles(Samples& samples, map<string, vector<Allele*> >& alleleGroups, int allowedAlleleTypes, int haplotypeLength, Allele& refallele) {

    map<string, int> seqCounts;
    bool multipleAllelesWithIdenticalAlts = false;
    string refseq = currentReferenceHaplotype();
    ++seqCounts[refseq];
    for (map<string, vector<Allele*> >::iterator a = alleleGroups.begin(); a != alleleGroups.end(); ++a) {
        Allele& allele = *a->second.front();
        if (seqCounts[allele.alternateSequence] > 0) {
            multipleAllelesWithIdenticalAlts = true;
            break;
        } else {
            ++seqCounts[allele.alternateSequence];
        }
    }

    if (multipleAllelesWithIdenticalAlts) {
        homogenizeAlleles(alleleGroups, refseq, refallele);
        getAlleles(samples, allowedAlleleTypes, haplotypeLength, false, true);
        alleleGroups.clear();
        groupAlleles(samples, alleleGroups);  // groups by alternate sequence
    }

}

// adjusts the registered alignment and contained alleles so that one allele
// covers the entire haplotype window
// returns a vector of pointers to alleles generated in this process
// alleles which are discarded are not explicitly removed, but 'squashed',
// which triggers their collection later
bool RegisteredAlignment::fitHaplotype(int haplotypeStart, int haplotypeLength, Allele*& aptr, bool allowPartials) {

    // if the read overlaps the haplotype window,
    // generate one Allele to describe the read in that region
    // and "squash" the unused ones
    vector<Allele*> newAllelesPtr;
    vector<Allele> newAlleles;

    int haplotypeEnd = haplotypeStart + haplotypeLength;

    //if (containedAlleleTypes == ALLELE_REFERENCE) {
    //    return false;
    //}
    /*
    cerr << "start: " << start << " end: " << end << endl;
    cerr << "haplotypestart: " << haplotypeStart << " haplotypeend: " << haplotypeEnd << endl;
    cerr << "registered alignment alleles," << endl << alleles << endl;
    */

    // save and bail out if we can't construct a haplotype allele
    vector<Allele> savedAlleles = alleles;

    if ((allowPartials && (start <= haplotypeEnd || end >= haplotypeStart))
        || (start <= haplotypeStart && end >= haplotypeEnd)) {
        vector<Allele>::iterator a = alleles.begin();
        //cerr << "trying to find overlapping haplotype alleles for the range " << haplotypeStart << " to " << haplotypeEnd << endl;
        //cerr << alleles << endl;
        while (a+1 != alleles.end()) {
            if (a->position <= haplotypeStart && a->position + a->referenceLength > haplotypeStart) {
                break;
            }
            ++a;
        }
        if (!(a->position <= haplotypeStart && a->position + a->referenceLength > haplotypeStart)) {
            return false;
        }
        vector<Allele>::iterator b = alleles.begin();
        while (b + 1 != alleles.end()) {
            if (b->position < haplotypeEnd && b->position + b->referenceLength >= haplotypeEnd) {
                break;
            }
            ++b;
        }
        if (!(b->position < haplotypeEnd && b->position + b->referenceLength >= haplotypeEnd)) {
            return false; // nothing to do here
        }

        // do not attempt to build haplotype alleles where there are non-contiguous reads
        /*
        for (vector<Allele>::iterator p = alleles.begin(); p != alleles.end(); ++p) {
            if (p != alleles.begin()) {
                if (p->position != (p - 1)->position + (p - 1)->referenceLength) {
                    cerr << "non-contiguous reads, cannot construct haplotype allele" << endl;
                    return true;
                }
            }
        }
        */

        // conceptually it will be easier to work on the haplotype obs if the reference alleles match the haplotype specification
        //if (a == b && a->isReference()) {
            // break the reference observation
            //cerr << "we just have a reference allele" << endl;
            //return true;
        //}

        string seq;
        vector<pair<int, string> > cigar;
        vector<short> quals;

        // now "a" should overlap the start of the haplotype block, and "b" the end
        //cerr << "block start overlaps: " << *a << endl;
        //cerr << "block end overlaps: " << *b << endl;
        //cerr << "haplotype start: " << haplotypeStart << endl;

        for (vector<Allele>::iterator p = a; p != (b+1); ++p) {
            if (p->isNull()) return false; // can't assemble across NULL alleles
        }

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
            //cerr << "subtracting " << haplotypeEnd - b->position << " from start " << newAllele << endl;
            newAllele.subtractFromStart(haplotypeEnd - b->position, seq, cigar, quals);
            if (isUnflankedIndel(newAllele)) {
                if (b + 1 != alleles.end()) {
                    ++b;
                }
            } else {
                b->subtractFromEnd(b->position + b->referenceLength - haplotypeEnd, seq, cigar, quals);
                newAlleles.push_back(newAllele);
            }
        }

        // now, for everything between a and b, merge them into one allele
        while (a != b) {
            vector<pair<int, string> > cigarV = splitCigar(a->cigar);
            vector<Allele>::iterator p = a + 1;
            // update the quality of the merged allele in the same way as we do
            // for complex events
            if (!a->isReference() && !a->isNull())  {
                p->quality = min(a->quality, p->quality);  // note that phred and log are inverted
                p->lnquality = max(a->lnquality, p->lnquality);
            }
            p->addToStart(a->alternateSequence, cigarV, a->baseQualities);
            a->squash();
            ++a;
        }

        // remove any 0-length alleles, these are useless
        // this operation requires independent removal of references to these alleles (e.g. registeredAlleles.clear())
        alleles.erase(remove_if(alleles.begin(), alleles.end(), isEmptyAllele), alleles.end());

        for (vector<Allele>::iterator p = newAlleles.begin(); p != newAlleles.end(); ++p) {
            alleles.push_back(*p);
        }

        AllelePositionCompare apcomp;
        sort(alleles.begin(), alleles.end(), apcomp);

        // now the pointers have changed, so find the allele we want... again!!!!!!
        //cerr << "registered alignment alleles, after haplotype construction," << endl << alleles << endl;
        bool hasHaplotypeAllele = false;
        bool dividedIndel = false;
        for (vector<Allele>::iterator p = alleles.begin(); p != alleles.end(); ++p) {
            // fix the "base"
            if (!p->isReference()) {
                p->update(haplotypeLength);
            }
            //cerr << *p << endl;
            if (p->position == haplotypeStart && p->position + p->referenceLength == haplotypeEnd) {
                aptr = &*p;
                if (isUnflankedIndel(*p)) {
                    hasHaplotypeAllele = false;
                    dividedIndel = true;
                } else {
                    hasHaplotypeAllele = true;
                }
                break;
            }
        }

        if (hasHaplotypeAllele) {
            //cerr << "registered alignment alleles after (pass)," << endl << alleles << endl;
            return true;
        } else {
            if (!allowPartials) {
                alleles = savedAlleles; // reset alleles
            }
            //cerr << "registered alignment alleles after (fail)," << endl << alleles << endl;
            return false;
            //assert(hasHaplotypeAllele);
        }

    } else {
        cerr << "registered alignment alleles after (pass)," << endl << alleles << endl;
        return true;
    }

}

void AlleleParser::buildHaplotypeAlleles(
    vector<Allele>& alleles,
    Samples& samples,
    map<string, vector<Allele*> >& alleleGroups,
    // provides observation group counts, counts of partial observations
    map<string, vector<Allele*> >& partialObservationGroups,
    map<Allele*, set<Allele*> >& partialObservationSupport,
    int allowedAlleleTypes) {

    int haplotypeLength = 1;
    for (vector<Allele>::iterator a = alleles.begin(); a != alleles.end(); ++a) {
        Allele& allele = *a;
        if (allele.isReference()) continue;
        // check if there are any complex alleles
        if (allele.referenceLength > haplotypeLength) {
            DEBUG("reference length of " << allele << " is " << allele.referenceLength
                  << " so extending haplotype");
            haplotypeLength = allele.referenceLength;
        }
        // check if we are embedded in a repeat structure
        if (allele.repeatRightBoundary > currentPosition + haplotypeLength) {
            DEBUG("right boundary " << allele.repeatRightBoundary << " for " << allele << " is past "
                  << currentPosition + haplotypeLength);
            haplotypeLength = allele.repeatRightBoundary - currentPosition;
        }
    }

    // return here if we have no registered alignments
    if (registeredAlignments.empty()) return;

    // always attempt to determine haplotype length in this fashion
    {

        DEBUG("haplotype length is " << haplotypeLength);

        // NB: for indels in tandem repeats, if the indel sequence is
        // derived from the repeat structure, build the haplotype
        // across the entire repeat pattern.  This ensures we actually
        // can discriminate between reference and indel/complex
        // alleles in the most common misalignment case.  For indels
        // that match the repeat structure, we have cached the right
        // boundary of the repeat.  We build the haplotype to the
        // maximal boundary indicated by the present alleles.

        int oldHaplotypeLength = haplotypeLength;
        do {
            oldHaplotypeLength = haplotypeLength;

            // rebuild samples
            samples.clear();

            long int maxAlignmentEnd = registeredAlignments.rbegin()->first;
            for (long int i = currentPosition+1; i < maxAlignmentEnd; ++i) {
                deque<RegisteredAlignment>& ras = registeredAlignments[i];
                for (deque<RegisteredAlignment>::iterator r = ras.begin(); r != ras.end(); ++r) {
                    RegisteredAlignment& ra = *r;
                    if ((ra.start > currentPosition && ra.start < currentPosition + haplotypeLength)
                        || (ra.end > currentPosition && ra.end < currentPosition + haplotypeLength)) {
                        Allele* aptr;
                        bool allowPartials = true;
                        ra.fitHaplotype(currentPosition, haplotypeLength, aptr, allowPartials);
                        for (vector<Allele>::iterator a = ra.alleles.begin(); a != ra.alleles.end(); ++a) {
                            registeredAlleles.push_back(&*a);
                        }
                    }
                }
            }

            getAlleles(samples, allowedAlleleTypes, haplotypeLength, true, true);
            alleleGroups.clear();
            groupAlleles(samples, alleleGroups);
            alleles = genotypeAlleles(alleleGroups, samples, parameters.onlyUseInputAlleles);
            for (vector<Allele>::iterator a = alleles.begin(); a != alleles.end(); ++a) {
                Allele& allele = *a;
                if (!allele.isReference()) {
                    long int alleleend = (allele.position + allele.referenceLength);
                    // this adjustment forces reference observations to overlap the ends of the indels
                    //if (allele.isInsertion() || allele.isDeletion()) {
                    //    alleleend += 1;
                    //}
                    long int hapend = max((long int) alleleend,
                                          allele.repeatRightBoundary);
                    /*
                    cerr << currentPosition + haplotypeLength << " vs " << alleleend
                         << " end " << hapend << " ? " << allele.position + allele.referenceLengthFromCigar()
                         << " hapend for " << allele << endl;
                    */
                    if (hapend > currentPosition + haplotypeLength) {
                        DEBUG("adjusting haplotype length to " << hapend - currentPosition
                              << " to overlap allele end " << alleleend
                              << " or right repeat boundary " << allele.repeatRightBoundary
                              << " " << allele);
                        haplotypeLength = hapend - currentPosition;
                    }
                }
            }
        } while (haplotypeLength != oldHaplotypeLength); // && haplotypeLength < parameters.maxHaplotypeLength);


        // TODO?
        //haplotypeLength = min(parameters.maxHaplotypeLength, haplotypeLength);

        // TODO adjust haplotypes over indels to include +1 bp on 3' end
        // this will force reference observations across the entire allele

        // for each non-reference allele within the haplotype length of this
        // position, adjust the length and reference sequences of the adjacent
        // alleles
        DEBUG("fitting haplotype block " << currentPosition << " to " << currentPosition + haplotypeLength << ", " << haplotypeLength << "bp");

        lastHaplotypeLength = haplotypeLength;

        registeredAlleles.clear();
        samples.clear();

        vector<Allele*> haplotypeObservations;
        getCompleteObservationsOfHaplotype(samples, haplotypeLength, haplotypeObservations);
        addToRegisteredAlleles(haplotypeObservations);
        DEBUG("added to registered alleles");

        // add partial observations
        // first get all the alleles up to the end of the haplotype window
        vector<Allele*> partialHaplotypeObservations;
        if (parameters.usePartialObservations && haplotypeLength > 1) {
            getPartialObservationsOfHaplotype(samples, haplotypeLength, partialHaplotypeObservations);
        }
        DEBUG("got partial observations of haplotype");
        //addToRegisteredAlleles(partialHaplotypeObservations);
        // now align the sequences of these alleles to the haplotype alleles
        // and put them into the partials bin in each sample

        // correct quality and alternate sequence for reference
        for (vector<Allele*>::iterator h = haplotypeObservations.begin(); h != haplotypeObservations.end(); ++h) {
            if ((*h)->position == currentPosition && (*h)->referenceLength == haplotypeLength) {
                (*h)->currentBase = (*h)->alternateSequence;
                (*h)->setQuality();
                (*h)->update(haplotypeLength);
                if ((*h)->isReference()) {  // HACK.. undoes damage of update() call
                    (*h)->currentBase = (*h)->alternateSequence;
                }
            }
        }
        for (vector<Allele*>::iterator p = partialHaplotypeObservations.begin(); p != partialHaplotypeObservations.end(); ++p) {
            (*p)->currentBase = (*p)->alternateSequence;
            (*p)->setQuality();
            (*p)->update(haplotypeLength);
        }
        DEBUG("done updating");

        if (parameters.debug) {
            cerr << "refr_seq\t" << currentPosition << "\t\t" << reference.getSubSequence(currentSequenceName, currentPosition, haplotypeLength) << endl;
            for (vector<Allele*>::iterator h = haplotypeObservations.begin(); h != haplotypeObservations.end(); ++h) {
                if ((*h)->position == currentPosition && (*h)->referenceLength == haplotypeLength) {
                    cerr << "haplo_obs\t" << (*h)->position << "\t" << (*h)->lnquality << "\t"
                        //<< (*h)->currentBase << "\t"
                         << string(max((long int)0,(*h)->position-currentPosition), ' ')
                         << (*h)->alternateSequence << "\t" << *h << endl;
                }
            }
            for (vector<Allele*>::iterator p = partialHaplotypeObservations.begin(); p != partialHaplotypeObservations.end(); ++p) {
                if ((*p)->position >= currentPosition && (*p)->position < currentPosition+haplotypeLength) {
                    cerr << "part_obs\t" << (*p)->position << "\t" << (*p)->lnquality << "\t"
                        //<< (*p)->currentBase << "\t"
                         << string(max((long int)0,(*p)->position-currentPosition), ' ')
                         << (*p)->alternateSequence << "\t" << *p << endl;
                }
            }
        }

        // now re-get the alleles
        getAlleles(samples, allowedAlleleTypes, haplotypeLength, false, true);

        // re-group the alleles using groupAlleles()
        alleleGroups.clear();
        groupAlleles(samples, alleleGroups);

        /*
        if (parameters.debug) {
            DEBUG("after re-grouping alleles");
            for (Samples::iterator s = samples.begin(); s != samples.end(); ++s) {
                cerr << s->first << endl;
                for (Sample::iterator t = s->second.begin(); t != s->second.end(); ++t) {
                    cerr << t->first << " " << t->second << endl << endl;
                }
            }
        }
        */

        Allele refAllele = genotypeAllele(ALLELE_REFERENCE,
                                          reference.getSubSequence(currentSequenceName, currentPosition, haplotypeLength),
                                          haplotypeLength,
                                          convert(haplotypeLength)+"M",
                                          haplotypeLength,
                                          currentPosition);

        // are there two alleles with the same alt sequence?
        // if so, homogenize them, and then re-sort the alleles

        // ensure uniqueness of registered alleles
        sort(registeredAlleles.begin(), registeredAlleles.end());
        registeredAlleles.erase(unique(registeredAlleles.begin(), registeredAlleles.end()), registeredAlleles.end());

        removeDuplicateAlleles(samples, alleleGroups, allowedAlleleTypes, haplotypeLength, refAllele);

        alleles = genotypeAlleles(alleleGroups, samples, parameters.onlyUseInputAlleles, haplotypeLength);

        // require all complete observations to effectively cover the same amount of sequence
        // basically, the "probe" length should be the same or we will incur bias when generating likelihoods
        // should these be put into the partial observations bin?

        int maxAlleleLength = haplotypeLength;
        for (vector<Allele>::iterator a = alleles.begin(); a != alleles.end(); ++a) {
            // get max allele length
            if (a->alternateSequence.size() > maxAlleleLength) maxAlleleLength = a->alternateSequence.size();
        }

        // bound this to 50bp so as to not drop out reference obs when we have long insertions directly encoded in the reads
        maxAlleleLength = min(50, maxAlleleLength);
        //cerr << "max allele length is " << maxAlleleLength << " but haplotype length = " << haplotypeLength << endl;
        // XXX make work for deletions as well
        if (maxAlleleLength > haplotypeLength) {
            //cerr << "max allele length = " << maxAlleleLength << endl;
            removeAllelesWithoutReadSpan(registeredAlleles, maxAlleleLength, haplotypeLength);
            samples.clear();
            // require that reference obs are over an equivalent amount of sequence as the max allele length
            getAlleles(samples, allowedAlleleTypes, haplotypeLength, false, true);
            alleleGroups.clear();
            groupAlleles(samples, alleleGroups);  // groups by alternate sequence
            // establish alleles again, now that we've filtered observations which don't have the required probe length
            alleles = genotypeAlleles(alleleGroups, samples, parameters.onlyUseInputAlleles, haplotypeLength);
        }

        // force the ref allele into the analysis, if it somehow isn't supported
        // this can happen where we don't have sufficient read span, such as in long deletions
        // or where our samples are homozygous for an alternate
        if (!parameters.useRefAllele) {
            vector<Allele> refAlleleVector;
            refAlleleVector.push_back(refAllele);
            alleles = alleleUnion(alleles, refAlleleVector);
        }

        // this is where we have established our genotype alleles
        /*
        for (vector<Allele>::iterator a = alleles.begin(); a != alleles.end(); ++a) {
            cerr << "genotype allele " << &*a << " " << *a << endl;
        }
        */

        // pick up observations that are potentially partial (not unambiguous)
        // the way to do this is to test the full observations as if they are partial, and if they
        // end up partially supporting multiple observations, removing them from the "complete" observations

        if (parameters.usePartialObservations && haplotypeLength > 1) {
            // check this out
            // here we are going to pass a set of full haplotype observations
            // and we'll remove now-partial obs from the full set
            samples.assignPartialSupport(alleles,
                                         haplotypeObservations,
                                         partialObservationGroups,
                                         partialObservationSupport,
                                         currentPosition,
                                         haplotypeLength);

            vector<Allele*> pureHaplotypeObservations;
            for (vector<Allele*>::iterator h = haplotypeObservations.begin(); h != haplotypeObservations.end(); ++h) {
                //if (partialObservationSupport.find(*h) != partialObservationSupport.end())
                //cerr << "partials for " << **h << " are " << partialObservationSupport[*h].size() << endl;
                if (partialObservationSupport.find(*h) != partialObservationSupport.end()
                    && partialObservationSupport[*h].size() > 0) {
                    DEBUG("full obs " << **h << " is actually partial and supports "
                          << partialObservationSupport[*h].size() << " alleles");
                    partialObservationSupport.erase(*h);
                    // and remove from partial observation groups?
                } else {
                    //cerr << "saving " << *h << endl;
                    pureHaplotypeObservations.push_back(*h);
                }
            }

            // now regenerate partial observation groups using updated partial support
            partialObservationGroups.clear();
            for (map<Allele*, set<Allele*> >::iterator p = partialObservationSupport.begin();
                 p != partialObservationSupport.end(); ++p) {
                set<Allele*>& supported = p->second;
                for (set<Allele*>::iterator s = supported.begin(); s != supported.end(); ++s) {
                    partialObservationGroups[(*s)->currentBase].push_back(p->first);
                }
            }

            // and keep only the pure haplotype observations for further use
            haplotypeObservations = pureHaplotypeObservations;

            addToRegisteredAlleles(haplotypeObservations);

            // clean up potential duplicates
            sort(registeredAlleles.begin(), registeredAlleles.end());
            registeredAlleles.erase(unique(registeredAlleles.begin(), registeredAlleles.end()), registeredAlleles.end());

            samples.clearFullObservations();
            getAlleles(samples, allowedAlleleTypes, haplotypeLength, false, true);
            alleleGroups.clear();
            groupAlleles(samples, alleleGroups);

            // stash partials for later
            addToRegisteredAlleles(partialHaplotypeObservations);

            for (vector<Allele*>::iterator p = partialHaplotypeObservations.begin(); p != partialHaplotypeObservations.end(); ++p) {
                (*p)->currentBase = (*p)->alternateSequence;
                (*p)->setQuality();
                (*p)->update(haplotypeLength);
            }

            // now add in partial observations collected from partially-overlapping reads
            if (!partialHaplotypeObservations.empty()) {
                samples.assignPartialSupport(alleles,
                                             partialHaplotypeObservations,
                                             partialObservationGroups,
                                             partialObservationSupport,
                                             currentPosition,
                                             haplotypeLength);
            }
        }

        registeredAlleles.clear();

        // reset registered alleles
        for (map<long unsigned int, deque<RegisteredAlignment> >::iterator ras = registeredAlignments.begin(); ras != registeredAlignments.end(); ++ras) {
            deque<RegisteredAlignment>& rq = ras->second;
            for (deque<RegisteredAlignment>::iterator rai = rq.begin(); rai != rq.end(); ++rai) {
                RegisteredAlignment& ra = *rai;
                for (vector<Allele>::iterator a = ra.alleles.begin(); a != ra.alleles.end(); ++a) {
                    registeredAlleles.push_back(&*a);
                }
            }
        }

        if (!parameters.useRefAllele) {
            vector<Allele> refAlleleVector;
            refAlleleVector.push_back(refAllele);
            alleles = alleleUnion(alleles, refAlleleVector);
        }

        //removeDuplicateAlleles(samples, alleleGroups, allowedAlleleTypes, haplotypeLength);
        //alleles = genotypeAlleles(alleleGroups, samples, parameters.onlyUseInputAlleles, haplotypeLength);

    }

    // hack......... TODO unhack this and set in Sample class
    samples.setSupportedAlleles();

    // processed flag..
    //unsetAllProcessedFlags();

    // redundant?

    // remove alleles which should no longer be considered
    //removePreviousAlleles(registeredAlleles, currentPosition);

    lastHaplotypeLength = haplotypeLength;

}

void AlleleParser::getCompleteObservationsOfHaplotype(Samples& samples, int haplotypeLength, vector<Allele*>& haplotypeObservations) {
    for (map<long unsigned int, deque<RegisteredAlignment> >::iterator ras = registeredAlignments.begin(); ras != registeredAlignments.end(); ++ras) {
        deque<RegisteredAlignment>& rq = ras->second;
        for (deque<RegisteredAlignment>::iterator rai = rq.begin(); rai != rq.end(); ++rai) {
            RegisteredAlignment& ra = *rai;
            Allele* aptr;
            // this guard prevents trashing allele pointers when getting partial observations
            //cerr << ra.start << " <= " << currentPosition << " && " << ra.end << " >= " << currentPosition + haplotypeLength << endl;
            if (ra.start <= currentPosition && ra.end >= currentPosition + haplotypeLength) {
                if (ra.fitHaplotype(currentPosition, haplotypeLength, aptr)) {
                    for (vector<Allele>::iterator a = ra.alleles.begin(); a != ra.alleles.end(); ++a) {
                        //cerr << a->position << " == " << currentPosition << " && " << a->referenceLength << " == " << haplotypeLength << endl;
                        if (a->position == currentPosition && a->referenceLength == haplotypeLength) {
                            haplotypeObservations.push_back(&*a);
                        }
                    }
                } /*else {
                    DEBUG("could not fit observation " << ra.name << " with alleles " << ra.alleles);
                    // the alleles have (possibly) been changed in fithaplotype, so add them to the registered alleles again
                    for (vector<Allele>::iterator a = ra.alleles.begin(); a != ra.alleles.end(); ++a) {
                        registeredAlleles.push_back(&*a);
                    }
                    }*/
            }
        }
    }
    DEBUG("got complete observations of haplotype");
}

void AlleleParser::unsetAllProcessedFlags(void) {
    for (map<long unsigned int, deque<RegisteredAlignment> >::iterator ras = registeredAlignments.begin(); ras != registeredAlignments.end(); ++ras) {
        deque<RegisteredAlignment>& rq = ras->second;
        for (deque<RegisteredAlignment>::iterator rai = rq.begin(); rai != rq.end(); ++rai) {
            RegisteredAlignment& ra = *rai;
            Allele* aptr;
            for (vector<Allele>::iterator a = ra.alleles.begin(); a != ra.alleles.end(); ++a) {
                a->processed = false; // re-trigger use of all alleles
            }
        }
    }
}


// process the next length bp of alignments, so as to get allele observations partially overlapping our calling window
void AlleleParser::getPartialObservationsOfHaplotype(Samples& samples, int haplotypeLength, vector<Allele*>& partials) {
    //cerr << "getting partial observations of haplotype from " << currentPosition << " to " << currentPosition + haplotypeLength << endl;
    vector<Allele*> newAlleles;

    bool gettingPartials = true;
    DEBUG("in AlleleParser::getPartialObservationsOfHaplotype, updating alignment queue");
    updateAlignmentQueue(currentPosition + haplotypeLength, newAlleles, gettingPartials);
    DEBUG("in AlleleParser::getPartialObservationsOfHaplotype, done updating alignment queue");

    vector<Allele*> otherObs;
    vector<Allele*> partialObs;
    // now get the partial obs
    // get the max alignment end position, iterate to there
    long int maxAlignmentEnd = registeredAlignments.rbegin()->first;
    for (long int i = currentPosition+1; i < maxAlignmentEnd; ++i) {
        DEBUG("getting partial observations of haplotype @" << i);
        deque<RegisteredAlignment>& ras = registeredAlignments[i];
        for (deque<RegisteredAlignment>::iterator r = ras.begin(); r != ras.end(); ++r) {
            RegisteredAlignment& ra = *r;
            if ((ra.start > currentPosition && ra.start < currentPosition + haplotypeLength)
		 || (ra.end > currentPosition && ra.end < currentPosition + haplotypeLength)) {
                Allele* aptr;
                bool allowPartials = true;
                ra.fitHaplotype(currentPosition, haplotypeLength, aptr, allowPartials);
                for (vector<Allele>::iterator a = ra.alleles.begin(); a != ra.alleles.end(); ++a) {
                    if (a->position >= currentPosition
                        && a->position < currentPosition+haplotypeLength
                        && !a->isNull()) {
                        //a->processed = false; // re-trigger use of all alleles
                        partials.push_back(&*a);
                    } else {
                        //a->processed = false;
                        otherObs.push_back(&*a);
                    }
                }
            } else {
                for (vector<Allele>::iterator a = ra.alleles.begin(); a != ra.alleles.end(); ++a) {
                    //a->processed = false;
                    otherObs.push_back(&*a);
                }
            }
        }
    }
    //addToRegisteredAlleles(partialObs);
    addToRegisteredAlleles(otherObs);
}

bool AlleleParser::getNextAlleles(Samples& samples, int allowedAlleleTypes) {
    long int nextPosition = currentPosition + lastHaplotypeLength;
    while (currentPosition < nextPosition) {
        if (!toNextPosition()) {
            return false;
        } else {
            // triggers cleanup
            if (justSwitchedTargets) {
                nextPosition = 0;
                justSwitchedTargets = false;
            }
            getAlleles(samples, allowedAlleleTypes);
        }
    }
    lastHaplotypeLength = 1;
    return true;
}

void AlleleParser::getAlleles(Samples& samples, int allowedAlleleTypes,
                              int haplotypeLength, bool getAllAllelesInHaplotype,
                              bool ignoreProcessedFlag) {
    Samples gvcf_held; // make some samples that by bass filtering for gvcf lines
    DEBUG2("getting alleles");
    samples.clear();
    // Commenting this out and replacinf with .clear() to relly empty it, it is more aloc, but no major change
    //for (Samples::iterator s = samples.begin(); s != samples.end(); ++s)
    //    s->second.clear();
    // TODO ^^^ this should be optimized for better scanning performance

    // if we have targets and are outside of the current target, don't return anything

    // add the reference allele to the analysis
    if (parameters.useRefAllele) {
        if (currentReferenceAllele) delete currentReferenceAllele; // clean up after last position
        currentReferenceAllele = referenceAllele(parameters.MQR, parameters.BQR);
        samples[referenceSampleName].clear();
        samples[referenceSampleName][currentReferenceAllele->currentBase].push_back(currentReferenceAllele);
        //alleles.push_back(currentReferenceAllele);
    }

    // get the variant alleles *at* the current position
    // and the reference alleles *overlapping* the current position
    for (vector<Allele*>::const_iterator a = registeredAlleles.begin(); a != registeredAlleles.end(); ++a) {
        Allele& allele = **a;
        //cerr << "getting alleles at position " << currentPosition << " with length " << haplotypeLength << " " << allele << endl;
        if (!ignoreProcessedFlag && allele.processed) continue;
        //cerr << "allele " << allele << endl;
        if (allowedAlleleTypes & allele.type
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
            allele.update(haplotypeLength);
            if(parameters.gVCFout){
                gvcf_held[allele.sampleID][allele.currentBase].push_back(*a); // store things incase
            }
            if (allele.quality >= parameters.BQL0 && allele.currentBase != "N"
                && (allele.isReference() || !allele.alternateSequence.empty())) { // filters haplotype construction chaff
                //cerr << "keeping allele " << allele << endl;
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
    if(samples.size() == 0 && parameters.gVCFout){
        samples = gvcf_held;  // if there are no non reference vals try to recover any allined values for gvcf if doing gvcf output!!
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
                                1,
                                currentPosition + 1,
                                0,
                                0,
                                base,
                                name,
                                name,
                                name,
                                sequencingTech,
                                true,
                                baseQ,
                                baseQstr,
                                mapQ,
                                false,
                                false,
                                false,
                                "1M",
                                NULL,
                                currentPosition,
                                currentPosition+1); // pair information
    allele->genotypeAllele = true;
    allele->baseQualities.push_back(baseQ);
    allele->update();
    return allele;
}

vector<Allele> AlleleParser::genotypeAlleles(
    map<string, vector<Allele*> >& alleleGroups, // alleles grouped by equivalence
    Samples& samples, // alleles grouped by sample
    bool useOnlyInputAlleles,
    int haplotypeLength
    ) {

    vector<pair<Allele, int> > unfilteredAlleles;

    DEBUG("getting genotype alleles");

    for (map<string, vector<Allele*> >::iterator group = alleleGroups.begin(); group != alleleGroups.end(); ++group) {
        // for each allele that we're going to evaluate, we have to have at least one supporting read with
        // map quality >= MQL1 and the specific quality of the allele has to be >= BQL1
        DEBUG("allele group " << group->first);
        vector<Allele*>& alleles = group->second;
        DEBUG(alleles);
        if (!allATGC(group->second.front()->alternateSequence)) {
            DEBUG("allele group contains partially-null observations, skipping");
            continue;
        }
        if (alleles.size() < parameters.minAltTotal) {
            DEBUG("allele group lacks sufficient observations in the whole population (min-alternate-total)");
            continue;
        }
        bool passesFilters = false;
        int qSum = 0;
        int mqSum = 0;
        for (vector<Allele*>::iterator a = alleles.begin(); a != alleles.end(); ++a) {
            DEBUG2("allele " << **a);
            Allele& allele = **a;
            qSum += allele.quality;
            mqSum += allele.mapQuality;
        }
        if (qSum >= parameters.minSupportingAlleleQualitySum && mqSum >= parameters.minSupportingMappingQualitySum) {
            Allele& allele = *(alleles.front());
            int length = allele.length;
            int reflength = allele.referenceLength;
            string altseq = allele.alternateSequence;
            if (allele.type == ALLELE_REFERENCE) {
                length = haplotypeLength;
                reflength = haplotypeLength;
                if (haplotypeLength == 1) {
                    altseq = currentReferenceBase;
                } else {
                    altseq = reference.getSubSequence(currentSequenceName, currentPosition, haplotypeLength);
                }
            }
            unfilteredAlleles.push_back(make_pair(genotypeAllele(allele.type,
                                                                 altseq,
                                                                 length,
                                                                 allele.cigar,
                                                                 reflength,
                                                                 allele.position,
                                                                 allele.repeatRightBoundary), qSum));
        }
    }
    DEBUG("found genotype alleles");

    map<Allele, int> filteredAlleles;

    DEBUG("filtering genotype alleles which are not supported by at least " << parameters.minAltCount
           << " observations comprising at least " << parameters.minAltFraction << " of the observations in a single individual");
    for (vector<pair<Allele, int> >::iterator p = unfilteredAlleles.begin();
         p != unfilteredAlleles.end(); ++p) {

        Allele& genotypeAllele = p->first;
        int qSum = p->second;
        DEBUG("genotype allele: " << genotypeAllele << " qsum " << qSum);

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
                      << " in individual " << s->first << " (" << observationCount << " obs)" <<  " and fraction "
                      << (float) alleleCount / (float) observationCount);
                filteredAlleles[genotypeAllele] = qSum;
                break;
                //out << *genotypeAllele << endl;
            }
        }
    }
    DEBUG("filtered genotype alleles");


    vector<Allele> resultAlleles;
    vector<Allele> resultIndelAndMNPAlleles;

    //string refBase = currentReferenceBaseString();
    // XXX XXX XXX
    string refBase = currentReferenceHaplotype();

    if (parameters.useBestNAlleles == 0) {
        // this means "use everything"
        bool hasRefAllele = false;
        for (map<Allele, int>::iterator p = filteredAlleles.begin();
             p != filteredAlleles.end(); ++p) {
            if (p->first.currentBase == refBase)
                hasRefAllele = true;
            DEBUG("adding allele to result alleles " << p->first.currentBase);
            resultAlleles.push_back(p->first);
        }
        // and add the reference allele if we need it
        if (parameters.forceRefAllele && !hasRefAllele) {
            DEBUG("including reference allele");
            // XXX TODO change to get the haplotype of the reference sequence
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

        DEBUG("getting " << parameters.useBestNAlleles << " best SNP alleles, and all other alleles");
        bool hasRefAllele = false;
        for (vector<pair<Allele, int> >::iterator a = sortedAlleles.begin(); a != sortedAlleles.end(); ++a) {
            Allele& allele = a->first;
            if (allele.currentBase == refBase) {
                hasRefAllele = true;
            }
            /*            if (allele.type & (ALLELE_DELETION | ALLELE_INSERTION | ALLELE_MNP | ALLELE_COMPLEX)) {
                DEBUG("adding allele to result alleles " << allele.currentBase);
                resultIndelAndMNPAlleles.push_back(allele);
            } else {
                DEBUG("adding allele to SNP alleles " << allele.currentBase);
            }
            */
            DEBUG("adding allele to result alleles " << allele.currentBase);
            resultAlleles.push_back(allele);
            DEBUG("allele quality sum " << a->second);
        }
        DEBUG("found " << sortedAlleles.size() << " SNP/ref alleles of which we now have " << resultAlleles.size() << endl
               << "and " << resultIndelAndMNPAlleles.size() << " INDEL and MNP alleles");

        // if we have reached the limit of allowable alleles, and still
        // haven't included the reference allele, include it
        if (parameters.forceRefAllele && !hasRefAllele) {
            DEBUG("including reference allele in analysis");
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

    // this needs to be fixed in a big way
    // the alleles have to be put into the local haplotype structure
    if (inputVariantAlleles.find(currentRefID) != inputVariantAlleles.end()) {
        map<long int, vector<Allele> >::iterator v = inputVariantAlleles[currentRefID].find(currentPosition);
        if (v != inputVariantAlleles[currentRefID].end()) {
            vector<Allele>& inputalleles = v->second;
            for (vector<Allele>::iterator a = inputalleles.begin(); a != inputalleles.end(); ++a) {
                DEBUG("evaluating input allele " << *a);
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
    }
    // remove non-unique alleles after

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

map<string, int> AlleleParser::repeatCounts(long int position, const string& sequence, int maxsize) {
    map<string, int> counts;
    for (int i = 1; i <= maxsize; ++i) {
        // subseq here i bases
        string seq = sequence.substr(position, i);
        // go left.

        int j = position - i;
        int leftsteps = 0;
        while (j >= 0 && seq == sequence.substr(j, i)) {
            j -= i;
            ++leftsteps;
        }

        // go right.
        j = position;

        int rightsteps = 0;
        while (j + i <= sequence.size() && seq == sequence.substr(j, i)) {
            j += i;
            ++rightsteps;
        }
        // if we went left and right a non-zero number of times,
        if (leftsteps + rightsteps > 1) {
            counts[seq] = leftsteps + rightsteps;
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

bool AlleleParser::isRepeatUnit(const string& seq, const string& unit) {

    if (seq.size() % unit.size() != 0) {
        return false;
    } else {
        int maxrepeats = seq.size() / unit.size();
        for (int i = 0; i < maxrepeats; ++i) {
            if (seq.substr(i * unit.size(), unit.size()) != unit) {
                return false;
            }
        }
        return true;
    }

}

bool AlleleParser::hasInputVariantAllelesAtCurrentPosition(void) {
    if (inputVariantAlleles.find(currentRefID) != inputVariantAlleles.end()) {
        map<long int, vector<Allele> >::iterator v = inputVariantAlleles[currentRefID].find(currentPosition);
        if (v != inputVariantAlleles[currentRefID].end()) {
            return true;
        }
    }
    return false;
}

bool operator<(const AllelicPrimitive& a, const AllelicPrimitive& b) {
    return a.ref < b.ref && a.alt < b.alt;
}

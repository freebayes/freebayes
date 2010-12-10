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


// XXX TODO change these void functions to bool

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
        bamMultiReader.SetIndexCacheMode(BamToolsIndex::NoIndexCaching);
    }

    if (parameters.useStdin) {
        if (!bamMultiReader.Open(parameters.bams, false, false, false)) {
            ERROR("Could not read BAM data from stdin");
            exit(1);
        }
    } else {
        if ( !bamMultiReader.Open(parameters.bams, true, false, true) ) {
            if ( !bamMultiReader.Open(parameters.bams, false, false, false) ) {
                ERROR("Could not open input BAM files");
                exit(1);
            } else {
                ERROR("Opened BAM reader without index file, jumping is disabled.");
                if (!targets.empty()) {
                    ERROR("Targets specified but no BAM index file provided.");
                    ERROR("FreeBayes cannot jump through targets in BAM files without BAM index files, exiting.");
                    ERROR("Please generate a BAM index file either .bai (standard) or .bti (bamtools), e.g.:");
                    ERROR("    \% bamtools index -bti -in <bam_file>   # for high-performance .bti index");
                    ERROR("    \% bamtools index -in <bam_file>        # for standard .bai index");
                    ERROR("    \% samtools index <bam_file>            # for standard .bai index");
                    //exit(1);
                }
            }
        }
    }
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

    // retrieve header information

    string bamHeader = bamMultiReader.GetHeaderText();

    vector<string> headerLines = split(bamHeader, '\n');

    for (vector<string>::const_iterator it = headerLines.begin(); it != headerLines.end(); ++it) {

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
        ERROR(string(80, '-') << endl
             //--------------------------------------------------------------------------------
           << "Warning: No sample file given, and no @RG tags found in BAM header." << endl
           << "All alignments from all input files will be assumed to come from the same" << endl
           << "individual.  To group alignments by sample, you must add read groups and sample" << endl 
           << "names to your alignments.  You can do this using ./scripts/sam_add_rg.pl in the" << endl
           << "freebayes source tree, or by specifying read groups and sample names when you" << endl
           << "prepare your sequencing data for alignment." << endl
           << string(80, '-'));
        sampleList.push_back("unknown");
        readGroupToSampleNames["unknown"] = "unknown";
        oneSampleAnalysis = true;
    }

}

void AlleleParser::writeVcfHeader(ostream& out) {

    time_t rawtime;
    struct tm * timeinfo;
    char datestr [80];

    time(&rawtime);
    timeinfo = localtime(&rawtime);

    strftime(datestr, 80, "%Y%m%d %X", timeinfo);

    out << "##format=VCFv4.0" << endl
            << "##fileDate=" << datestr << endl
            << "##source=freebayes" << endl
            << "##reference=" << parameters.fasta << endl
            << "##phasing=none" << endl
            << "##commandline=\"" << parameters.commandline << "\"" << endl
            << "##INFO=NS,1,Integer,\"total number of samples\"" << endl
            << "##INFO=ND,1,Integer,\"total number of non-duplicate samples\"" << endl
            << "##INFO=DP,1,Integer,\"total read depth at this base\"" << endl
            << "##INFO=AC,1,Integer,\"total number of alternate alleles in called genotypes\"" << endl

            // these are req'd
            << "##FORMAT=GT,1,String,\"Genotype\"" << endl // g
            << "##FORMAT=GQ,1,Integer,\"Genotype Quality\"" << endl // phred prob of genotype
            << "##FORMAT=DP,1,Integer,\"Read Depth\"" << endl // NiBAll[ind]
            << "##FORMAT=HQ,2,Integer,\"Haplotype Quality\"" << endl
            << "##FORMAT=QiB,1,Integer,\"Total base quality\"" << endl
            << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
            << join(sampleList, "\t")
            << endl;

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
    assert(parameters.useStdin); // this should only be used in the case that we are reading from stdin
    DEBUG2("loading reference sequence overlapping first alignment");
    currentPosition = alignment.Position;
    currentSequenceStart = alignment.Position;
    currentSequenceName = referenceIDToName[alignment.RefID];
    currentRefID = alignment.RefID;
    if (currentTarget == NULL) {
        currentTarget = new BedTarget(currentSequenceName, currentSequenceStart + 1, 0);
    } else {
        if (currentSequenceName != currentTarget->seq)
            currentTarget->seq = currentSequenceName;
        if (currentSequenceStart + 1 != currentTarget->left)
            currentTarget->left = currentSequenceStart + 1;
    }
    DEBUG2("reference.getSubSequence("<< currentSequenceName << ", " << currentSequenceStart << ", " << alignment.AlignedBases.length() << ")");
    currentSequence = uppercase(reference.getSubSequence(currentSequenceName, currentSequenceStart, alignment.AlignedBases.length()));
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
    currentSequence += uppercase(reference.getSubSequence(reference.sequenceNameStartingWith(currentTarget->seq), 
                                                 (currentTarget->right - 1) + basesAfterCurrentTarget,
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

    // if we have a region specified, use it to generate a target
    if (!parameters.region.empty()) {
        // drawn from bamtools_utilities.cpp, modified to suit 1-based context, no end sequence

        string region = parameters.region;
        string startSeq;
        int startPos;
        int stopPos;

        size_t foundFirstColon = region.find(":");

        // we only have a single string, use the whole sequence as the target
        if (foundFirstColon == string::npos) {
            startSeq = region;
            startPos = 1;
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
                    stopPos = atoi(region.substr(foundRangeDots + 2).c_str()) + 1; // end-inclusive
                } else {
                    stopPos = reference.sequenceLength(startSeq);
                }
            }
        }

        BedTarget bd(startSeq,
                    (startPos == 1) ? 1 : startPos,
                    (stopPos == -1) ? reference.sequenceLength(startSeq) : stopPos);
        DEBUG2("will process reference sequence " << startSeq << ":" << bd.left << ".." << bd.right);
        targets.push_back(bd);
    }

    // if we have a targets file, use it...
    // if target file specified use targets from file
    if (!parameters.targets.empty()) {

        DEBUG("Making BedReader object for target file: " << parameters.targets << " ...");

        BedReader bedReader(parameters.targets);

        if (!bedReader.is_open()) {
            ERROR("Unable to open target file: " << parameters.targets << "... terminating.");
            exit(1);
        }

        targets = bedReader.entries();

        // check validity of targets wrt. reference
        for (vector<BedTarget>::iterator e = targets.begin(); e != targets.end(); ++e) {
            BedTarget& bd = *e;
            if (bd.left < 1 || bd.right >= reference.sequenceLength(bd.seq)) {
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

        if (targets.size() == 0) {
            ERROR("Could not load any targets from " << parameters.targets);
            exit(1);
        }

        bedReader.close();

        DEBUG("done");

    }

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
        BedTarget bd(refName, 1, refData.RefLength);
        DEBUG2("will process reference sequence " << refName << ":" << bd.left << ".." << bd.right);
        targets.push_back(bd);
    }
}

// meant to be used when we are reading from stdin, to check if we are within targets
bool AlleleParser::inTarget(void) {
    if (targets.empty()) {
        return true;  // everything is in target if we don't have targets
    } else {
        // otherwise, scan through the targets to establish if the coordinate lies within any
        for (vector<BedTarget>::iterator e = targets.begin(); e != targets.end(); ++e) {
            BedTarget& bd = *e;
            if (currentTarget->seq == bd.seq && bd.left <= currentPosition && currentPosition <= bd.right) {
                return true;
            }
        }
    }
    return false;
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


    // initialization
    openTraceFile();
    openFailedFile();
    openOutputFile();

    loadFastaReference();
    // check how many targets we have specified
    loadTargets();
    // when we open the bam files we can use the number of targets to decide if
    // we should load the indexes
    openBams();
    loadBamReferenceSequenceNames();
    getSampleNames();

}

AlleleParser::~AlleleParser(void) {
    // close trace file?  seems to get closed properly on object deletion...
    if (currentReferenceAllele != NULL) delete currentReferenceAllele;

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

string AlleleParser::referenceSubstr(long double pos, unsigned int len) {
    return uppercase(reference.getSubSequence(currentSequenceName, floor(pos), len));
}

bool AlleleParser::isCpG(string& altbase) {
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

RegisteredAlignment& AlleleParser::registerAlignment(BamAlignment& alignment, RegisteredAlignment& ra, string sampleName) {

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
                         << currentTarget->seq << ":" << (long unsigned int) currentPosition + 1 << endl;
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
                         << currentTarget->seq << ":" << (long unsigned int) currentPosition + 1 << endl;
                    abort();
                }

                // this regex-generated nonsense comment maintained for laughs:
                //
                // THIS DEBUGIC IS ILDEBUGICAL!!!
                //
                //
                // record mismatch if we have a mismatch here
                if (b != sb) {
                    if (firstMatch < csp) {
                        // TODO ; verify that the read and reference sequences *do* match
                        int length = csp - firstMatch;
                        string matchingSequence = currentSequence.substr(csp - length, length);
                        //cerr << "matchingSequence " << matchingSequence << endl;
                        string readSequence = rDna.substr(rp - length, length);
                        //cerr << "readSequence " << readSequence << endl;
                        string qualstr = rQual.substr(rp - length, length);
                        // record 'reference' allele for last matching region
                        ra.alleles.push_back(Allele(ALLELE_REFERENCE,
                                currentTarget->seq, sp - length, &currentPosition, &currentReferenceBase, length, 
                                matchingSequence, readSequence, sampleName, alignment.Name,
                                !alignment.IsReverseStrand(), alignment.MapQuality, qualstr,
                                alignment.MapQuality));
                        DEBUG2(ra.alleles.back());
                    }
                    // register mismatch
                    if (qual >= parameters.BQL2) {
                        ra.mismatches++;  // increment our mismatch counter if we're over BQL2
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
                    string matchingSequence = currentSequence.substr(csp - length, length);
                    string readSequence = rDna.substr(rp - length, length);
                    string qualstr = rQual.substr(rp - length, length);
                    AlleleType mismatchtype = (length == 1) ? ALLELE_SNP : ALLELE_MNP;
                    long double lqual = averageQuality(qualstr);
                    ra.alleles.push_back(Allele(mismatchtype, currentTarget->seq, sp - length, &currentPosition,
                                &currentReferenceBase, length, matchingSequence, readSequence,
                                sampleName, alignment.Name, !alignment.IsReverseStrand(), lqual, qualstr, alignment.MapQuality));
                    DEBUG2(ra.alleles.back());
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
                string matchingSequence = currentSequence.substr(csp - length, length);
                string readSequence = rDna.substr(rp - length, length);
                string qualstr = rQual.substr(rp - length, length);
                AlleleType mismatchtype = (length == 1) ? ALLELE_SNP : ALLELE_MNP;
                long double lqual = averageQuality(qualstr);
                ra.alleles.push_back(Allele(mismatchtype, currentTarget->seq, sp - length, &currentPosition,
                            &currentReferenceBase, length, matchingSequence, readSequence,
                            sampleName, alignment.Name, !alignment.IsReverseStrand(), lqual, qualstr, alignment.MapQuality));
                DEBUG2(ra.alleles.back());
            // or, if we are not in a mismatch, construct the last reference allele of the match
            } else if (firstMatch < csp) {
                int length = csp - firstMatch;
                string matchingSequence = currentSequence.substr(csp - length, length);
                string readSequence = rDna.substr(rp - length, length);
                string qualstr = rQual.substr(rp - length, length);
                ra.alleles.push_back(Allele(ALLELE_REFERENCE,
                        currentTarget->seq, sp - length, &currentPosition, &currentReferenceBase, length, 
                        matchingSequence, readSequence, sampleName, alignment.Name,
                        !alignment.IsReverseStrand(), alignment.MapQuality, qualstr,
                        alignment.MapQuality));
                DEBUG2(ra.alleles.back());
            }
        } else if (t == 'D') { // deletion

            // because deletions have no quality information,
            // use the surrounding sequence quality as a proxy
            // to provide quality scores of equivalent magnitude to insertions,
            // take N bp, centered on the position of the deletion
            // this logic prevents overflow of the read
            int spanstart;
            if (l > rQual.size()) {
                l = rQual.size();
                spanstart = 0;
            } else {
                // set lower bound to 0
                if (rp < (l / 2)) {
                    spanstart = 0;
                } else {
                    spanstart = rp - (l / 2);
                }
                // set upper bound to the string length
                if (spanstart + l > rQual.size()) {
                    spanstart = rQual.size() - l;
                }
            }

            string qualstr = rQual.substr(spanstart, l);

            long double qual = sumQuality(qualstr);

            if (qual >= parameters.BQL2) {
                ra.mismatches += l;
            }

            if (qual >= parameters.BQL2) {
                for (int i=0; i<l; i++) {
                    indelMask[sp - alignment.Position + i] = true;
                }
            }

            ra.alleles.push_back(Allele(ALLELE_DELETION,
                    currentTarget->seq, sp, &currentPosition, &currentReferenceBase, l,
                    currentSequence.substr(csp, l), "", sampleName, alignment.Name,
                    !alignment.IsReverseStrand(), qual, qualstr,
                    alignment.MapQuality));
            DEBUG2(ra.alleles.back());

            sp += l;  // update sample position
            csp += l;

        } else if (t == 'I') { // insertion

            string qualstr = rQual.substr(rp, l);

            //long double qual = jointQuality(qualstr);
            // calculate average quality of the insertion
            long double qual = sumQuality(qualstr); //(long double) *max_element(quals.begin(), quals.end());

            if (qual >= parameters.BQL2) {
                ra.mismatches += l;
            }
            if (qual >= parameters.BQL2) {
                indelMask[sp - alignment.Position] = true;
                indelMask[sp - alignment.Position + 1] = true; // TODO wrong, but mirrors bambayes behavior
            }

            ra.alleles.push_back(Allele(ALLELE_INSERTION,
                    currentTarget->seq, sp - 0.5, &currentPosition, &currentReferenceBase, l, "", rDna.substr(rp, l),
                    sampleName, alignment.Name, !alignment.IsReverseStrand(), qual,
                    qualstr, alignment.MapQuality));
            DEBUG2(ra.alleles.back());

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
            int start = (allele.position - alignment.Position);
            int end = start + allele.length;
            vector<bool>::iterator im = indelMask.begin();
            // if there is anything masked, store it, otherwise just leave the
            // indelMask on this alignment empty, which means, "no masking" in
            // Allele::masked()
            for (vector<bool>::iterator q = im + start; q != im + end; ++q) {
                if (*q) { // there is a masked element
                    allele.indelMask.resize(allele.length);
                    copy(im + start, im + end, allele.indelMask.begin());
                    break;
                }
            }
        }
    }

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
            if (!currentAlignment.IsPrimaryAlignment())
                continue;

            // otherwise, get the sample name and register the alignment to generate a sequence of alleles
            // we have to register the alignment to acquire some information required by filters
            // such as mismatches

            // initially skip reads with low mapping quality (what happens if MapQuality is not in the file)
            if (currentAlignment.MapQuality >= parameters.MQL0) {
                // extend our cached reference sequence to allow processing of this alignment
                extendReferenceSequence(currentAlignment);
                string sampleName = readGroupToSampleNames[readGroup];
                // decomposes alignment into a set of alleles
                // here we get the deque of alignments ending at this alignment's end position
                deque<RegisteredAlignment>& rq = registeredAlignments[currentAlignment.GetEndPosition(true)];
                // and insert the registered alignment into that deque
                rq.push_front(RegisteredAlignment(currentAlignment));
                RegisteredAlignment& ra = rq.front();
                registerAlignment(currentAlignment, ra, sampleName);
                // backtracking if we have too many mismatches
                if (ra.mismatches > parameters.RMU) {
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
        if (currentPosition > position + (*allele)->length) {
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

void AlleleParser::removeNonOverlappingAlleles(vector<Allele*>& alleles) {
    for (vector<Allele*>::iterator allele = alleles.begin(); allele != alleles.end(); ++allele) {
        if ((*allele)->type == ALLELE_REFERENCE) {
            //if (currentPosition >= (*allele)->position + (*allele)->length) {
            if (!(((*allele)->position + (*allele)->length - currentPosition) > 0.5)) {
                *allele = NULL;
            }
        } else { // snps, insertions, deletions
            if (currentPosition >= (*allele)->position) {
                *allele = NULL;
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

    // load first target if we have targets and have not loaded the first
    if (currentTarget == NULL && !parameters.useStdin && !targets.empty()) {
        if (!loadTarget(&targets.front())) {
            ERROR("Could not load first target");
            return false;
        }
        if (!getFirstAlignment()) {
            ERROR("Could not get first alignment from target");
            return false;
        }
        // XXX hack
        clearRegisteredAlignments();
        currentSequenceStart = currentAlignment.Position;
        currentSequenceName = referenceIDToName[currentAlignment.RefID];
        currentRefID = currentAlignment.RefID;
        currentPosition = (currentPosition < currentAlignment.Position) ? currentAlignment.Position : currentPosition;
        currentSequence = uppercase(reference.getSubSequence(currentSequenceName, currentSequenceStart, currentAlignment.Length));
    // stdin, no targets cases
    } else if (currentTarget == NULL && (parameters.useStdin || targets.empty())) {
        if (!getFirstAlignment()) {
            ERROR("Could not get first alignment from target");
            return false;
        }
        clearRegisteredAlignments();
        loadReferenceSequence(currentAlignment); // this seeds us with new reference sequence
    // we've reached the end of file, or stdin, before we ever started?
    } else if (parameters.useStdin || targets.empty()) {
        ERROR("could not read any alignments");
        return false;
    } else {
        bool ok = false;
        // step through targets until we get to one with alignments
        while (currentTarget != &targets.back()) {
            if (!loadTarget(++currentTarget)) {
                continue;
            }
            if (ok = getFirstAlignment()) {
                break;
            }
        }
        if (!ok) return false; // last target and couldn't get alignment

        clearRegisteredAlignments();
        currentSequenceStart = currentAlignment.Position;
        currentSequenceName = referenceIDToName[currentAlignment.RefID];
        currentRefID = currentAlignment.RefID;
        currentPosition = (currentPosition < currentAlignment.Position) ? currentAlignment.Position : currentPosition;
        currentSequence = uppercase(reference.getSubSequence(currentSequenceName, currentSequenceStart, currentAlignment.Length));
        //clearRegisteredAlignments();
        //loadReferenceSequence(currentAlignment);
    }

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

    int refSeqID = bamMultiReader.GetReferenceID(currentTarget->seq);

    DEBUG2("reference sequence id " << refSeqID);

    DEBUG2("setting new position " << currentTarget->left);
    currentPosition = currentTarget->left - 1;

    if (!bamMultiReader.SetRegion(refSeqID, currentTarget->left - 1, refSeqID, currentTarget->right - 1)) {
        ERROR("Could not SetRegion to " << currentTarget->seq << ":" << currentTarget->left << ".." << currentTarget->right);
        return false;
    }

    // now that we've jumped, reset the hasMoreAlignments counter
    hasMoreAlignments = true;

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
        ERROR("Could not find any mapped reads in target region " << currentTarget->seq << ":" << currentTarget->left << ".." << currentTarget->right);
        return false;
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
bool AlleleParser::toNextPosition(void) {

    if (currentTarget == NULL) {
        DEBUG2("loading first target");
        if (!toNextTarget()) {
            return false;
        }
    } 
    else {
        if (!parameters.allowIndels) {
            ++currentPosition;
        } else {
            currentPosition += 0.5;
        }
    }

    if (!parameters.useStdin && !targets.empty() && currentPosition >= currentTarget->right - 1) { // time to move to a new target
        DEBUG2("next position " << currentPosition + 1 <<  " outside of current target right bound " << currentTarget->right);
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
        if (registeredAlignments.empty() && currentRefID != currentAlignment.RefID) {
            clearRegisteredAlignments();
            loadReferenceSequence(currentAlignment);
            justSwitchedTargets = true;
        } else if (!hasMoreAlignments && registeredAlignments.empty()) {
            DEBUG("no more alignments in input");
            return false;
        }
    }

    currentReferenceBase = currentReferenceBaseChar();

    // handle the case in which we don't have targets but in which we've switched reference sequence

    DEBUG2("processing position " << (long unsigned int) currentPosition + 1 << " in sequence " << currentTarget->seq);
    updateAlignmentQueue();
    DEBUG2("updating registered alleles");
    updateRegisteredAlleles(); // this removes unused left-flanking sequence

    // if we have alignments which ended at the previous base, erase them and their alleles
    map<long unsigned int, deque<RegisteredAlignment> >::iterator f = registeredAlignments.find(currentPosition - 2);
    if (f != registeredAlignments.end()) {
        registeredAlignments.erase(f);
    }
    // so we have to make sure it's still there (this matters in low-coverage)
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

    //DEBUG2("processing position " << currentPosition + 1 << " in sequence " << currentTarget->seq);
    return true;
}

bool AlleleParser::getNextAlleles(Samples& samples, int allowedAlleleTypes) {
    if (toNextPosition()) {
        getAlleles(samples, allowedAlleleTypes);
        return true;
    } else {
        return false;
    }
}

void AlleleParser::getAlleles(Samples& samples, int allowedAlleleTypes) {

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
                removeNonOverlappingAlleles(alleles); // removes alleles which no longer overlap our current position
                updateAllelesCachedData(alleles);  // calls allele.update() on each Allele*
                removeFilteredAlleles(alleles); // removes alleles which are filtered at this position, 
                                                  // and requeues them for processing by unsetting their 'processed' flag
            }
            sample.sortAlleles();
        }
    }

    // if we have targets and are outside of the current target, don't return anything


    // add the reference allele to the analysis
    if (parameters.useRefAllele) {
        if (currentReferenceAllele != NULL) delete currentReferenceAllele; // clean up after last position
        currentReferenceAllele = referenceAllele(parameters.MQR, parameters.BQR);
        samples[currentTarget->seq].clear();
        samples[currentTarget->seq][currentReferenceAllele->currentBase].push_back(currentReferenceAllele);
        //alleles.push_back(currentReferenceAllele);
    }

    // get the variant alleles *at* the current position
    // and the reference alleles *overlapping* the current position
    for (vector<Allele*>::const_iterator a = registeredAlleles.begin(); a != registeredAlleles.end(); ++a) {
        Allele& allele = **a;
        if (!allele.processed
                && allowedAlleleTypes & allele.type
                && (
                    (allele.type == ALLELE_REFERENCE 
                      && currentPosition >= allele.position 
                      //&& currentPosition < allele.position + allele.length // 0-based, means position + length - 1 is the last included base
                      && (allele.position + allele.length - currentPosition) > 0.5)
                  || 
                    (allele.position == currentPosition)
                    ) 
                ) {
            allele.update();
            if (allele.quality >= parameters.BQL0 && !allele.masked() && allele.currentBase != "N") {
                samples[allele.sampleID][allele.currentBase].push_back(*a);
                allele.processed = true;
            }
        }
    }

    vector<string> samplesToErase;
    // now remove empty alleles from our return so as to not confuse processing
    for (Samples::iterator s = samples.begin(); s != samples.end(); ++s) {

        const string& name = s->first;
        Sample& sample = s->second;

        // move updated alleles to the right bin
        sample.sortAlleles();

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
        if (empty) {
            samplesToErase.push_back(name);
        }
    }

    for (vector<string>::iterator name = samplesToErase.begin(); name != samplesToErase.end(); ++name) {
        samples.erase(*name);
    }

    DEBUG2("done getting alleles");

}

Allele* AlleleParser::referenceAllele(int mapQ, int baseQ) {
    string base = string(1, currentReferenceBase);
    //string name = reference.filename;
    string name = currentTarget->seq; // this behavior matches old bambayes
    string baseQstr = "";
    //baseQstr += qualityInt2Char(baseQ);
    Allele* allele = new Allele(ALLELE_REFERENCE, 
            currentTarget->seq,
            currentPosition,
            &currentPosition, 
            &currentReferenceBase,
            1, base, base, name, name,
            true, baseQ,
            baseQstr,
            mapQ);
    allele->genotypeAllele = true;
    allele->baseQualities.push_back(baseQ);
    allele->update();
    return allele;
}

vector<Allele> AlleleParser::genotypeAlleles(
        map<string, vector<Allele*> >& alleleGroups, // alleles grouped by equivalence
        Samples& samples, // alleles grouped by sample
        vector<Allele>& allGenotypeAlleles      // all possible genotype alleles, 
                                                // to add back alleles if we don't have enough to meet our minimum allele count
        ) {

    vector<pair<Allele, int> > unfilteredAlleles;

    DEBUG2("getting genotype alleles");

    for (map<string, vector<Allele*> >::iterator group = alleleGroups.begin(); group != alleleGroups.end(); ++group) {
        // for each allele that we're going to evaluate, we have to have at least one supporting read with
        // map quality >= MQL1 and the specific quality of the allele has to be >= BQL1
        DEBUG2("allele group " << group->first);
        bool passesFilters = false;
        int qSum = 0;
        vector<Allele*>& alleles = group->second;
        for (vector<Allele*>::iterator a = alleles.begin(); a != alleles.end(); ++a) {
            DEBUG2("allele " << **a);
            Allele& allele = **a;
            qSum += allele.quality;
        }
        Allele& allele = *(alleles.front());
        int length = (allele.type == ALLELE_REFERENCE || allele.type == ALLELE_SNP) ? 1 : allele.length;
        unfilteredAlleles.push_back(make_pair(genotypeAllele(allele.type, allele.currentBase, length), qSum));
    }
    DEBUG2("found genotype alleles");

    map<Allele, int> filteredAlleles;

    DEBUG2("filtering genotype alleles which are not supported by at least " << parameters.minAltCount 
            << " observations comprising at least " << parameters.minAltFraction << " of the observations in a single individual");
    for (vector<pair<Allele, int> >::iterator p = unfilteredAlleles.begin();
            p != unfilteredAlleles.end(); ++p) {

        Allele& genotypeAllele = p->first;
        int qSum = p->second;

        for (Samples::iterator s = samples.begin(); s != samples.end(); ++s) {
            Sample& sample = s->second; 
            int alleleCount = 0;
            Sample::iterator c = sample.find(genotypeAllele.currentBase);
            if (c != sample.end())
                alleleCount = c->second.size();
            int observationCount = sample.observationCount();
            if (alleleCount >= parameters.minAltCount 
                    && ((float) alleleCount / (float) observationCount) >= parameters.minAltFraction) {
                DEBUG2(genotypeAllele << " has support of " << alleleCount 
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

    string refBase = string(1, currentReferenceBase);

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
            resultAlleles.insert(resultAlleles.begin(), genotypeAllele(ALLELE_REFERENCE, refBase, 1));
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

        DEBUG2("getting N best alleles");
        bool hasRefAllele = false;
        for (vector<pair<Allele, int> >::iterator a = sortedAlleles.begin(); a != sortedAlleles.end(); ++a) {
            if (a->first.currentBase == refBase)
                hasRefAllele = true;
            DEBUG2("adding allele to result alleles " << a->first.currentBase);
            resultAlleles.push_back(a->first);
            if (resultAlleles.size() >= parameters.useBestNAlleles) {
                DEBUG2("found enough alleles, breaking out");
                break;
            }
        }
        DEBUG2("found " << sortedAlleles.size() << " alleles of which we now have " << resultAlleles.size());

        // if we have reached the limit of allowable alleles, and still
        // haven't included the reference allele, include it
        if (parameters.forceRefAllele && !hasRefAllele) {
            DEBUG2("including reference allele");
            resultAlleles.insert(resultAlleles.begin(), genotypeAllele(ALLELE_REFERENCE, refBase, 1));
        }

        // if we now have too many alleles (most likely one too many), get rid of some
        while (resultAlleles.size() > parameters.useBestNAlleles) {
            resultAlleles.pop_back();
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

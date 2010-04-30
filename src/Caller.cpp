#include "Caller.h"
#include "multichoose.h" // includes generic functions, so it must be included here
                         // otherwise we will get a linker error
                         // see: http://stackoverflow.com/questions/36039/templates-spread-across-multiple-files
                         // http://www.cplusplus.com/doc/tutorial/templates/ "Templates and Multi-file projects"

// local helper debugging macros to improve code readability
#define LOG(msg) \
    if (parameters->record) { logFile << msg << endl; } \
    if (parameters->debug) { cerr << msg << endl; }

// lower-priority messages
#define LOG2(msg) \
    if (parameters->record) { logFile << msg << endl; } \
    if (parameters->debug2) { cerr << msg << endl; }

// must-see error messages
#define ERROR(msg) \
    if (parameters->record) { logFile << msg << endl; } \
    cerr << msg << endl;

using namespace std;


// XXX TODO change these void functions to bool

// open BAM input file
void Caller::openBam(void) {

    LOG("Opening BAM fomat alignment input file: " << parameters->bam << " ...")
  
    bamReader.Open(parameters->bam.c_str(), (parameters->bam + ".bai").c_str());

    LOG(" done");
}

void Caller::openLogFile(void) {

    logFile.open(parameters->log.c_str(), ios::out);
    if (parameters->record) {
        if (parameters->debug) {cerr << "Opening log file: " << parameters->log << " ...";}

        if (!logFile) {
            ERROR(" unable to open file: " << parameters->log);
            exit(1);
        }
        if (parameters->debug) {cerr << " done." << endl;}
    }
    // previously we wrote the command to the logfile here
}

// read sample list file or get sample names from bam file header
void Caller::getSampleNames(void) {
    // If a sample file is given, use it.  But otherwise process the bam file
    // header to get the sample names.
    //
    if (parameters->samples != "") {
        ifstream sampleFile(parameters->samples.c_str(), ios::in);
        if (! sampleFile) {
            ERROR("unable to open file: " << parameters->samples);
            exit(1);
        }
        boost::regex patternSample("^(\\S+)\\s*(.*)$");
        boost::regex re("\\s+");
        boost::regex pr("^(\\S+):(\\S+)$");
        boost::smatch match;
        string line;
        while (getline(sampleFile, line)) {
            // if proper line
            if (boost::regex_search(line, match, patternSample)) {
                // assign content
                string s = match[1];
                LOG2("found sample " << s);
                sampleList.push_back(s);
            }
        }
    } else { // no samples file given, read from BAM file header for sample names
        // retrieve header information
        LOG("no sample list file given, attempting to read sample names from bam file");

        string bamHeader = bamReader.GetHeaderText();

        vector<string> headerLines;
        boost::split(headerLines, bamHeader, boost::is_any_of("\n"));

        for (vector<string>::const_iterator it = headerLines.begin(); it != headerLines.end(); ++it) {

            // get next line from header, skip if empty
            string headerLine = *it;
            if ( headerLine.empty() ) { continue; }

            // lines of the header look like:
            // "@RG     ID:-    SM:NA11832      CN:BCM  PL:454"
            //                     ^^^^^^^\ is our sample name
            if ( headerLine.find("@RG") == 0 ) {
                vector<string> readGroupParts;
                boost::split(readGroupParts, headerLine, boost::is_any_of("\t "));
                vector<string> nameParts;
                boost::split(nameParts, readGroupParts.at(2), boost::is_any_of(":"));
                string name = nameParts.back();
                //mergedHeader.append(1, '\n');
                LOG2("found sample " << name);
                sampleList.push_back(name);
            }
        }
    }
}

void Caller::loadBamReferenceSequenceNames(void) {

    //--------------------------------------------------------------------------
    // read reference sequences from input file
    //--------------------------------------------------------------------------

    // store the names of all the reference sequences in the BAM file
    referenceSequences = bamReader.GetReferenceData();

    LOG("Number of ref seqs: " << bamReader.GetReferenceCount());

}


void Caller::loadFastaReference(void) {

    // XXX we don't technically need to load the whole thing into memory
    // the FastaReference class will let us grab individual sequences and subsequences 
    // very fast from the file
    // thus cutting down on memory requirements...
    //
    // perhaps a good approach is to open the file here
    // and then get long subsequences at runtime
    // 
    // this keeps our memory requirements low, and will allow us to operate unmodified on more systems

    LOG("processing fasta reference " << parameters->fasta);

    //--------------------------------------------------------------------------
    // process input fasta file
    //--------------------------------------------------------------------------
    // This call loads the reference and reads any index file it can find.
    // If it can't find an index file for the reference, it will attempt to
    // generate one alongside it.

    reference = new FastaReference(parameters->fasta);

    fastaReferenceSequenceCount = 0;

    int id = 0;

    //--------------------------------------------------------------------------
    // load ref seq names into hash
    //--------------------------------------------------------------------------
    for(map<string, FastaIndexEntry>::const_iterator it = reference->index->begin(); 
          it != reference->index->end(); ++it) {

        FastaIndexEntry entry = it->second;

        // we split out the first word of the full sequence name for keying our sequences
        // as often the full sequence names are overkill
        vector<string> sequenceNameParts;
        boost::split(sequenceNameParts, entry.name, boost::is_any_of(" "));
        string name = sequenceNameParts.front();

        // get the reference names in this vector
        referenceSequenceNames.push_back(name);  // WARNING: these are broken; no order guarantees
        referenceSequenceNameToID[name] = id++;
        ++fastaReferenceSequenceCount;

    }

    LOG(" done.");

}

void Caller::loadReferenceSequence(int seqID) {
    LOG2("loading reference sequence " << seqID);
    string name = reference->sequenceNameStartingWith(referenceSequenceNames[seqID]);
    currentSequence = reference->getSequence(name);
}

void Caller::loadReferenceSequence(string seqName, int start, int length) {
    LOG2("loading reference subsequence " << seqName << " from " << start << " to " << start + length);
    string name = reference->sequenceNameStartingWith(seqName);
    currentSequence = reference->getSubSequence(name, start, length);
}

void Caller::loadReferenceSequence(BedData* target) {
    loadReferenceSequence(target->seq, target->left, target->right - target->left);
}

// intended to load all the sequence covered by reads which overlap our current target
// this lets us process the reads fully, checking for suspicious reads, etc.
// but does not require us to load the whole sequence
void Caller::loadReferenceSequence(BedData* target, int before, int after) {
    basesBeforeCurrentTarget = before;
    basesAfterCurrentTarget = after;
    loadReferenceSequence(target->seq, target->left - before, target->right - target->left + after + before);
}

void Caller::loadTargets(void) {

  // process target region input file
  
  // if target file specified use targets from file
  int targetCount = 0;
  // if we have a targets file, use it...
  if (parameters->targets != "") {
    
    LOG("Making BedReader object for target file: " << parameters->targets << " ...");
    
    BedReader bedReader(parameters->targets);
    
    if (! bedReader.isOpen()) {
      ERROR("Unable to open target file: " << parameters->targets << "... terminating.");
      exit(1);
    }
    
    BedData bd;
    while (bedReader.getNextEntry(bd)) {
        if (parameters->debug2) {
            cerr << bd.seq << "\t" << bd.left << "\t" << bd.right << "\t" << bd.desc << endl;
        }
        // TODO add back check that the right bound isn't out of bounds
        string seqName = reference->sequenceNameStartingWith(bd.seq);
        if (bd.left < 1 || bd.right < bd.left || bd.right >= reference->sequenceLength(seqName)) {
          ERROR("Target region coordinate outside of reference sequence bounds... terminating.");
          exit(1);
        }
        targets.push_back(bd);
        targetCount++;
    }
    
    bedReader.close();

    LOG("done");

  // otherwise analyze all reference sequences from BAM file
  } else {
    RefVector::iterator refIter = referenceSequences.begin();
    RefVector::iterator refEnd  = referenceSequences.end();
    for( ; refIter != refEnd; ++refIter) {
      RefData refData = *refIter;
      string refName = refData.RefName;
      BedData bd;
      bd.seq = refName;
      bd.left = 1;
      bd.right = refData.RefLength;
      targets.push_back(bd);
      targetCount++;
    }
  }

  LOG("Number of target regions: " << targetCount);

}

void Caller::initializeOutputFiles(void) {

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // open output file(s)
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  // open report output file
  //----------------------------------------------------------------------------

  // report
  LOG("opening report output file for writing: " << parameters->rpt << "...");

  // open output streams
  bool outputRPT, outputVCF; // for legibility

  if (parameters->rpt != "") {
      outputRPT = true;
      rptFile.open(parameters->rpt.c_str());
      if (!rptFile) {
        ERROR(" unable to open file: " << parameters->rpt);
        exit(1);
      }
  } else { outputRPT = false; }

  if (parameters->vcf != "") {
      outputVCF = true;
      vcfFile.open(parameters->vcf.c_str());
      if (!vcfFile) {
        ERROR(" unable to open file: " << parameters->vcf);
        exit(1);
      }
  } else { outputVCF = false; }
  LOG(" done.");

  //----------------------------------------------------------------------------
  // write header information
  //----------------------------------------------------------------------------
  if (outputRPT) {
      rptFile << "# Complete list of parameter values:" << endl;
      rptFile << "#   --bam = " << parameters->bam << endl;
      rptFile << "#   --fasta = " << parameters->fasta << endl;
      rptFile << "#   --targets = " << parameters->targets << endl;
      rptFile << "#   --samples = " << parameters->samples << endl;
      rptFile << "#   --rpt = " << parameters->rpt << endl;
      rptFile << "#   --log = " << parameters->log << endl;
      rptFile << "#   --useRefAllele = " <<  ( parameters->useRefAllele ? "true" : "false" ) << endl;
      rptFile << "#   --forceRefAllele = " <<  ( parameters->forceRefAllele ? "true" : "false" ) << endl;
      rptFile << "#   --MQR = " << parameters->MQR << endl;
      rptFile << "#   --BQR = " << parameters->BQR << endl;
      rptFile << "#   --ploidy = " << parameters->ploidy << endl;
      rptFile << "#   --sampleNaming = " << parameters->sampleNaming << endl;
      rptFile << "#   --sampleDel = " << parameters->sampleDel << endl;
      rptFile << "#   --BQL0 = " << parameters->BQL0 << endl;
      rptFile << "#   --MQL0 = " << parameters->MQL0 << endl;
      rptFile << "#   --BQL1 = " << parameters->BQL1 << endl;
      rptFile << "#   --MQL1 = " << parameters->MQL1 << endl;
      rptFile << "#   --BQL2 = " << parameters->BQL2 << endl;
      rptFile << "#   --RMU = " << parameters->RMU << endl;
      rptFile << "#   --IDW = " << parameters->IDW << endl;
      rptFile << "#   --TH = " << parameters->TH << endl;
      rptFile << "#   --PVL = " << parameters->PVL << endl;
      rptFile << "#   --algorithm = " << parameters->algorithm << endl;
      rptFile << "#   --RDF = " << parameters->RDF << endl;
      rptFile << "#   --WB = " << parameters->WB << endl;
      rptFile << "#   --TB = " << parameters->TB << endl;
      rptFile << "#   --includeMonoB = " <<  ( parameters->includeMonoB ? "true" : "false" ) << endl;
      rptFile << "#   --TR = " << parameters->TR << endl;
      rptFile << "#   --I = " << parameters->I << endl;
      rptFile << "#   --debug = " <<  ( parameters->debug ? "true" : "false" ) << endl;
      rptFile << "#   --debug2 = " <<  ( parameters->debug2 ? "true" : "false" ) << endl;
      rptFile << "#" << endl;
  }

  
  if (outputVCF) {
      time_t rawtime;
      struct tm * timeinfo;
      char datestr [80];

      time(&rawtime);
      timeinfo = localtime(&rawtime);

      strftime(datestr, 80, "%Y%m%d %X", timeinfo);

      vcfFile << "##format=VCFv3.3" << endl
              << "##fileDate=" << datestr << endl
              << "##source=gigabayes" << endl
              << "##reference=1000GenomesPilot-NCBI36" << endl
              << "##phasing=none" << endl
              << "##notes=\"All FORMAT fields matching *i* (e.g. NiBAll, NiA) refer to individuals.\"" << endl
             
              << "##INFO=NS,1,Integer,\"total number of samples\"" << endl
              << "##INFO=ND,1,Integer,\"total number of non-duplicate samples\"" << endl
              << "##INFO=DP,1,Integer,\"total read depth at this base\"" << endl
              << "##INFO=AC,1,Integer,\"total number of alternate alleles in called genotypes\"" << endl
              //<< "##INFO=AN,1,Integer,\"total number of alleles in called genotypes\"" << endl

              // these are req'd
              << "##FORMAT=GT,1,String,\"Genotype\"" << endl // g
              << "##FORMAT=GQ,1,Integer,\"Genotype Quality\"" << endl // phred prob of genotype
              << "##FORMAT=DP,1,Integer,\"Read Depth\"" << endl // NiBAll[ind]
              << "##FORMAT=HQ,2,Integer,\"Haplotype Quality\"" << endl
              << "##FORMAT=QiB,1,Integer,\"Total base quality\"" << endl
              << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" 
              << boost::algorithm::join(sampleList, "\t")
              << endl;

  }
}

// initialization function
// sets up environment so we can start registering alleles
Caller::Caller(int argc, char** argv)
{
    parameters = new Parameters(argc, argv);

    // mathematical constants
    LOGFACTOR = log((long double)10.0) / ((long double)-10.0); 
    LN3 = log((long double)3.0);

    // initialization
    // NOTE: these void functions have side effects, and currently have to be called in this order
    // this separation is done to improve legibility and debugging
    // perhaps it will just increase confusion
    openLogFile();
    openBam();
    getSampleNames();
    loadFastaReference();
    loadBamReferenceSequenceNames();
    loadTargets();
    initializeOutputFiles();

    currentRefID = 0; // will get set properly via toNextRefID
    //toNextRefID(); // initializes currentRefID
    toFirstTargetPosition(); // initializes currentTarget, currentAlignment
}

Caller::~Caller(void) {
    delete parameters;
    delete reference;
}

// position of alignment relative to current sequence
int Caller::currentSequencePosition(const BamAlignment& alignment) {
    return (alignment.Position - currentTarget->left) + basesBeforeCurrentTarget;
}

// registeredalignment friend
ostream& operator<<(ostream& out, RegisteredAlignment& ra) {
    out << ra.alignment.Name << " " << ra.alignment.Position << endl
        << ra.alignment.QueryBases << endl
        << ra.alignment.Qualities << endl;
    for (vector<Allele>::const_iterator a = ra.alleles.begin(); a != ra.alleles.end(); ++a) {
        Allele allele = *a;
        out << allele;
    }
    return out;
}

RegisteredAlignment Caller::registerAlignment(BamAlignment& alignment) {

    RegisteredAlignment ra = RegisteredAlignment(alignment); // result

    string rDna = alignment.QueryBases;
    string rQual = alignment.Qualities;
    int rp = 0;  // read position, 0-based relative to read
    int csp = currentSequencePosition(alignment); // current sequence position, 0-based relative to currentSequence
    int sp = alignment.Position + 1;  // sequence position
              //   ^^^ conversion between 0 and 1 based index

    
    // extract sample name and information
    string readName = alignment.Name;
    string sampleName;
    if (! alignment.GetReadGroup(sampleName)) {
        /*cerr << "WARNING: Couldn't find read group id (@RG tag) for BAM Alignment " << alignment.Name
            << " ... attempting to read from read name" << endl;
            */
        SampleInfo sampleInfo = extractSampleInfo(readName, parameters->sampleNaming, parameters->sampleDel);
        sampleName = sampleInfo.sampleId;
    }

    LOG2("registering alignment " << rp << " " << csp << " " << sp << endl <<
         "alignment readName " << readName << endl <<
         "alignment sampleID " << sampleName << endl << 
         "alignment position " << alignment.Position << endl <<
         "alignment length " << alignment.Length << endl <<
         "alignment AlignedBases.size() " << alignment.AlignedBases.size() << endl <<
         "alignment end position " << alignment.Position + alignment.AlignedBases.size());

    LOG2(endl << rDna << endl << alignment.AlignedBases << endl << currentSequence.substr(csp, alignment.AlignedBases.size()));


    vector<CigarOp>::const_iterator cigarIter = alignment.CigarData.begin();
    vector<CigarOp>::const_iterator cigarEnd  = alignment.CigarData.end();
    for ( ; cigarIter != cigarEnd; ++cigarIter ) {
        unsigned int l = (*cigarIter).Length;
        char t = (*cigarIter).Type;
        LOG2("cigar item: " << t << l);

        if (t == 'S') { // soft clip
            rp += l;
        } else if (t == 'M') { // match or mismatch
            int lastMismatch = csp; // track the last mismatch, for recording 'reference' alleles
            for (int i=0; i<l; i++) {

                // extract aligned base
                string b;
                TRY { b = rDna.at(rp); } CATCH;

                // convert base quality value into short int
                short qual = qualityChar2ShortInt(rQual[rp]);

                // get reference allele
                string sb;
                TRY { sb = currentSequence.at(csp); } CATCH;

                // record match if we have a mismatch here
                // TODO instrument this to validate behavior
                if (b != sb) {
                    if (lastMismatch < csp) {
                        // TODO ; verify that the read and reference sequences *do* match
                        int length = csp - lastMismatch;
                        string matchingSequence = currentSequence.substr(lastMismatch, length);
                        string qualstr = rQual.substr(rp - length, length);
                        ra.alleles.push_back(Allele(ALLELE_REFERENCE,
                                    currentTarget->seq, sp - length, length, 
                                    matchingSequence, "", sampleName, alignment.Name,
                                    !alignment.IsReverseStrand(), -1, qualstr,
                                    alignment.MapQuality));
                    }
                    lastMismatch = csp;
                }

                // register mismatch
                if (b != sb && qual >= parameters->BQL2) {
                    // record 'reference' allele for last matching region
                    ra.mismatches++;
                    ra.alleles.push_back(Allele(ALLELE_SNP, currentTarget->seq, sp, 1, sb, b,
                            sampleName, alignment.Name, !alignment.IsReverseStrand(), qual, "", alignment.MapQuality));
                }

                // update positions
                ++sp;
                ++csp;
                ++rp;
            }
            if (lastMismatch < csp) {  // TODO use some param flag to trigger this operation
                // TODO ; verify that the read and reference sequences *do* match
                int length = csp - lastMismatch;
                string matchingSequence = currentSequence.substr(lastMismatch, length);
                string qualstr = rQual.substr(rp - length, length);
                ra.alleles.push_back(Allele(ALLELE_REFERENCE,
                            currentTarget->seq, sp - length, length, 
                            matchingSequence, "", sampleName, alignment.Name,
                            !alignment.IsReverseStrand(), -1, qualstr,
                            alignment.MapQuality));
            }
            // XXX what about 'N' s?
        } else if (t == 'D') { // deletion

            // extract base quality of left and right flanking non-deleted bases
            string qualstr = rQual.substr(rp, 2);

            // calculate joint quality of the two flanking qualities
            short qual = jointQuality(qualstr); // XXX was max, but joint makes more sense, right ?
            if (qual >= parameters->BQL2) {
                ra.alleles.push_back(Allele(ALLELE_DELETION,
                            currentTarget->seq, sp, l,
                            currentSequence.substr(csp, l), "", sampleName, alignment.Name,
                            !alignment.IsReverseStrand(), qual, qualstr,
                            alignment.MapQuality));
            }

            sp += l;  // update sample position
            csp += l;

        } else if (t == 'I') { // insertion
            
            string qualstr = rQual.substr(rp, l);

            // calculate joint quality, which is the probability that there are no errors in the observed bases
            short qual = jointQuality(qualstr);
            // register insertion + base quality with reference sequence
            if (qual >= parameters->BQL2) { // XXX this cutoff may not make sense for long indels... the joint quality is much lower than the        'average' quality
                ra.alleles.push_back(Allele(ALLELE_INSERTION,
                            currentTarget->seq, sp, l, "", rDna.substr(rp, l),
                            sampleName, alignment.Name, !alignment.IsReverseStrand(), qual,
                            qualstr, alignment.MapQuality));
            }

            rp += l;

        } // not handled, skipped region 'N's
    } // end cigar iter loop
    return ra;
}


void Caller::updateAlignmentQueue(void) {

    LOG2("updating alignment queue");

    // push to the front until we get to an alignment that doesn't overlap our
    // current position or we reach the end of available alignments
    // filter input reads; only allow mapped reads with a certain quality
    bool moreAlignments = true; // flag to catch BAM EOF
    int i = 0;
    while (moreAlignments && currentAlignment.Position <= currentPosition) {
    //while (moreAlignments && currentTarget->left < currentAlignment.Position) {
        // only process if mapped
        /*
        cerr << i++ << endl;
        cerr << "alignment position " << currentAlignment.Position << endl;
        cerr << "alignment end " << currentAlignment.Position + currentAlignment.Length << endl;
        cerr << "target limits " << currentTarget->left << " " << currentTarget->right << endl;
        cerr << ((currentAlignment.Position >= currentTarget->right) ? "out of bounds" : "in bounds") << endl;
        */
        if (currentAlignment.IsMapped()) {
            RegisteredAlignment ra = registerAlignment(currentAlignment);
            // filters
            // 'read mask' --- this just means "don't consider snps right next to potential indels
            //                 ... but it should be implemented
            // low mapping quality --- hmmm ... could calculate it with the jointQuality function?
            // duplicates --- tracked via each BamAlignment
            //for (vector<Allele>::const_iterator it = ra.alleles.begin(); it != ra.alleles.end(); ++it) {
            //    Allele a = *it; cout << a << endl << endl;
           // }
            if (!(ra.mismatches > parameters->RMU)) {
                registeredAlignmentQueue.push_front(ra);
            }
            // TODO collect statistics here...
        }
        moreAlignments &= bamReader.GetNextAlignment(currentAlignment);
    }

    LOG2("... finished pushing new alignments");

    // pop from the back until we get to an alignment that overlaps our current position
    if (registeredAlignmentQueue.size() > 0) {
        BamAlignment* alignment = &registeredAlignmentQueue.back().alignment;
        // is indexing (0 or 1 based) oK here?
        while (currentPosition > alignment->Position + alignment->Length && registeredAlignmentQueue.size() > 0) {
            LOG2("popping alignment");
            registeredAlignmentQueue.pop_back();
            if (registeredAlignmentQueue.size() > 0) {
                alignment = &registeredAlignmentQueue.back().alignment;
            } else {
                break;
            }
        }
    }

    LOG2("... finished popping old alignments");
}

// initialization function, should only be called via constructor
bool Caller::toFirstTargetPosition(void) {
    return loadTarget(&targets.front());
}

// steps our position/beddata/reference pointers through all positions in all
// targets, returns false when we are finished
//
// pushes and pulls alignments out of our queue of overlapping alignments via
// updateAlignmentQueue() as we progress

/// XXX TODO  ... we MUST deal with updating the target reference sequence here
//  ... also need to have various checks here; 
//  1) invalid reference id
//  2) no alignments in region
//  3) failed jump to target start
//  ...
//  TODO might want to generalize this into a jump function and a step function

bool Caller::toNextTarget(void) {

    LOG2("seeking to next valid target...");

    // if we are at the end of the list of targets in this reference sequence
    if (currentTarget == &targets.back()) {
        return false;
    } else {
        return loadTarget(++currentTarget);
    }

}

bool Caller::loadTarget(BedData* target) {

    currentTarget = target;

    LOG("processing target " << currentTarget->desc << " " <<
            currentTarget->seq << " " << currentTarget->left << " " <<
            currentTarget->right);
    LOG2("clearing alignment queue");
    registeredAlignmentQueue.clear(); // clear our alignment deque on jumps

    LOG2("loading target reference subsequence");
    BamAlignment alignment;
    int refSeqID = referenceSequenceNameToID[currentTarget->seq];
    LOG2("reference sequence id " << refSeqID);

    bool r = bamReader.Jump(refSeqID, currentTarget->left);
    r &= bamReader.GetNextAlignment(alignment);
    int left_gap = currentTarget->left - alignment.Position;

    r &= bamReader.Jump(refSeqID, currentTarget->right - 1);

    int maxPos = 0;
    do {
        r &= bamReader.GetNextAlignment(alignment);
        int newPos = alignment.Position + alignment.AlignedBases.size();
        maxPos = (newPos > maxPos) ? newPos : maxPos;
    } while (alignment.Position < currentTarget->right);

    int right_gap = maxPos - currentTarget->right;
    //cerr << "left_gap " << left_gap << " right gap " << right_gap << endl;
    loadReferenceSequence(currentTarget,
            (left_gap > 0) ? left_gap : 0,
            (right_gap > 0) ? right_gap : 0);

    LOG2("setting new position " << currentTarget->left);
    currentPosition = currentTarget->left;

    LOG2("jumping to first alignment in new target");
    r &= bamReader.Jump(refSeqID, currentTarget->left);
    r &= bamReader.GetNextAlignment(currentAlignment);

    return r;

}

// stepping
//
// if the next position is outside of target region
// seek to next target which is in-bounds for its sequence
// if none exist, return false
bool Caller::toNextTargetPosition(void) {

    ++currentPosition;
    // ... 0-base , 1-base ??? XXX this is borked
    if (currentPosition > currentTarget->right) { // time to move to a new target
        LOG2("next position " << currentPosition <<  " outside of current target right bound " << currentTarget->right);
        if (!toNextTarget()) {
            LOG("no more valid targets, finishing");
            return false;
        }
    }
    LOG2("processing position " << currentPosition << " in sequence " << currentTarget->seq);
    updateAlignmentQueue();
    return true;
}

bool Caller::getNextAlleles(vector<Allele>& alleles) {
    if (toNextTargetPosition()) {
        getAlleles(alleles);
        return true;
    } else {
        return false;
    }
}

// updates the passed vector with the current alleles at the caller's target position
void Caller::getAlleles(vector<Allele>& alleles) {

    // we used to just clear, but this seems inefficient (?)
    // alleles.clear();

    // remove alleles which aren't reference alleles
    // remove reference alleles which are no longer overlapping the current position
    for (vector<Allele>::iterator ai = alleles.begin(); ai != alleles.end(); ++ai) {
        Allele a = *ai;
        if (!(a.type == ALLELE_REFERENCE && currentPosition >= a.position && currentPosition < a.position + a.length)) {
            ai = alleles.erase(ai);
            if (ai == alleles.end())
                break;
        }
    }

    // get the variant alleles *at* the current position
    // and the reference alleles *overlapping* the current position
    
    for (deque<RegisteredAlignment>::const_iterator it = registeredAlignmentQueue.begin(); it != registeredAlignmentQueue.end(); ++it) {
        RegisteredAlignment ra = *it;
        for (vector<Allele>::const_iterator ai = ra.alleles.begin(); ai != ra.alleles.end(); ++ai) {
            Allele a = *ai;
            // for now only record the allele if it is at exactly the current position
            if (a.position == currentPosition)
                alleles.push_back(a);
        }
    }
}

// TODO allele sorting by sample
// which will enable tho following to work


///////////////////////////////////////////////////////////////////////////////
// math
///////////////////////////////////////////////////////////////////////////////

// p( observed alleles | genotypes )
// for all genotypes
//
// overview:
//
// the alleles vector is assumed to be from a single sample
// the genotype vector is the genotype at this position
// we are estimating p( observed allele | genotype )
//
// note that the observed alleles are not true alleles in that they are, by
// definition, only representative of a single base, whereas the true
// alleles have can be hetero-  or homozygous
//
// sum, for all possible true alleles in our observedAlleles
// the product p( observed allele | true allele ) * p( true allele | genotype )
// 
// (1) "observed allele probabilities"
//      This relates the estimated observation error rate and the
//      probability that we observe what is actually present in the sample.
// p( observed allele | true allele ) = 1 - 10 ^(-Q/10)  if observed allele == true allele
//                          ... and   = 10 ^(-Q/10) if observed allele != true allele
//             where Q is the quality value associated with the observed allele
//             nb: this looks like the benoulli distribution to me -eg
//
// (2) "true allele probabilities"
//     This describes the probability that our 'true' alleles accurately
//     sample the underlying genotype.
//  in the case of diploid individuals, this is approximated as a binomial probability
// p( true allele | genotype ) = m! / f!(m - f) * p ^ f
//  where m is the count of observations and f is successes ? TODO
//  however, for n-ploid individuals, this must be approximated as a multinomial probability:
//  TODO
//

// computational steps:
//
// (A) for each observed allele
//   step through the allele set vector
//   if we match a set of alleles, push back on that set
//   else, push back on a new set and push it onto the allele set vector
//
// (B) generate all possible true allele combinations
//  choose all possible allele combinations of n-ploidy.
//  Reduces to generating all multiset combinations, or k multichoos n from our
//  groups of observed alleles.
//  http://stackoverflow.com/questions/127704/algorithm-to-return-all-combinations-of-k-elements-from-n
//  http://stackoverflow.com/questions/561/using-combinations-of-sets-as-test-data#794
//
// (C) for all possible true allele combinations
//   sum the product of (1) and (2)
//
vector<double> Caller::probObservedAllelesGivenGenotype(vector<Allele> &observedAlleles, vector<vector<Allele> > &genotypes) {

    int ploidy = genotypes.front().size(); // ploidy is determined by the number of alleles in the genotypes
    vector<double> results;

    // (A)
    vector<vector<Allele> > alleleGroups = groupAlleles(observedAlleles, allelesEquivalent);
    
    // (B)
    // k multichoose n, or ploidy multichoose alleleGroups
    //vector<vector<vector<Allele> > > alleleMultiCombinations = multichoose(ploidy, alleleGroups);

    // (C) 
    for (vector<vector<Allele > >::iterator g = genotypes.begin(); g != genotypes.end(); ++g) {
        results.push_back(probAlleleComboGivenGenotype(alleleGroups, *g));
    }

    return results;

}

// the product p( observed allele | true allele ) * p( true allele | genotype )
// (1) "observed allele probabilities"
//      This relates the estimated observation error rate and the
//      probability that we observe what is actually present in the sample.
// p( observed allele | true allele ) = 1 - 10 ^(-Q/10)  if observed allele == true allele
//                          ... and   = 10 ^(-Q/10) if observed allele != true allele
//             where Q is the quality value associated with the observed allele
//
// (2) "true allele probabilities"
// multinomial probability
//   p( true allele | genotype ) = [ n! / ( n1! * n2! * ... nk! ) ] * ( p1n1 * p2n2 * . . . * pknk )
//   P = ( ploidy! / product( alleleCount! for alleleType in alleleCombo ) ) * product( (1/ploidy)^alleleCount for alleleType in alleleCombo )

double Caller::probAlleleComboGivenGenotype(vector<vector<Allele> > &alleleCombo, vector<Allele> &genotype) {

    int ploidy = genotype.size();

    LOG2(stringForAlleles(genotype));

    // counts indexed by genotype position
    //vector<int> sumQ (ploidy, 0);
    //vector<int> counts (ploidy, 0);
    int outCount = 0, inCount = 0;
    double inProb = 1;
    double outProb = 1;
    int alleleComboObsCount = 0;

    // (1)
    // p( observed allele | true allele ) = 1 - 10 ^(-Q/10)  if observed allele in true allele
    //                          ... and   = 10 ^(-Q/10) if observed allele != true allele
    //             where Q is the quality value associated with the observed allele

    double factAlleleCountProduct = 1;
    double probObsAllelesProduct = 1;
    for (vector<vector<Allele> >::iterator alleleObservations = alleleCombo.begin();
            alleleObservations != alleleCombo.end(); ++alleleObservations) {
        
        // NB: Not the greatest solution.  This should be improved.
        // Avoids processing the same observations more than once.
        // We can safely do this because multichoose provides the allele
        // combination multisets in lexographic order.
        if ((alleleObservations + 1) != alleleCombo.end() 
                && (alleleObservations + 1)->front().equivalent(alleleObservations->front())) {
            continue;
        }
    // (2)
    // P = ( ploidy! / product( alleleCount! for alleleType in alleleCombo ) ) * product( p_^alleleCount for alleleType in alleleCombo )
        alleleComboObsCount += alleleObservations->size();
        factAlleleCountProduct *= factorial(alleleObservations->size());
        probObsAllelesProduct *= pow(probChooseAlleleFromAlleles(alleleObservations->front(), genotype), 
                alleleObservations->size());

        for (vector<Allele>::iterator observation = alleleObservations->begin();
                observation != alleleObservations->end(); ++observation) {
            bool in = false;
            int i = 0;
            for (vector<Allele>::iterator g = genotype.begin(); g != genotype.end(); g++) {
                // process genotypes only once
                if ((g + 1) != genotype.end() 
                        && (g + 1)->equivalent(*g)) {
                    continue;
                }
                if (g->equivalent(*observation)) {
                    inProb *= 1 - phred2float(observation->Quality(currentPosition));
                    inCount++;
                    in = true;
                }
                ++i;
            }
            if (!in) {
                outProb *= phred2float(observation->Quality(currentPosition));
                outCount++;
            }
        }
    }

    // (1)
    LOG2("outCount = " << outCount);
    LOG2("outProb = " << outProb);
    LOG2("inCount = " << inCount);
    LOG2("inProb = " << inProb);
    double probAlleleObsGivenGenotype = outProb * inProb;
    LOG2("probAlleleObsGivenGenotype = " << probAlleleObsGivenGenotype);
    
    // (2)
    LOG2("factAlleleCountProduct = " << factAlleleCountProduct);
    LOG2("alleleComboObsCount = " << alleleComboObsCount);
    LOG2("probObsAllelesProduct = " << probObsAllelesProduct);
    double trueAlleleProbability = ( (double) factorial(alleleComboObsCount) / factAlleleCountProduct ) * probObsAllelesProduct;
    LOG2("trueAlleleProbability = " << trueAlleleProbability);

    // (1) * (2)
    return probAlleleObsGivenGenotype * trueAlleleProbability;
}

double probChooseAlleleFromAlleles(Allele &allele, vector<Allele> &alleles) {

    int matches = 0;

    for (vector<Allele>::iterator ai = alleles.begin(); ai != alleles.end(); ai++) {
        if (ai->equivalent(allele))
            ++matches;
    }

    return (double) matches / (double) alleles.size();

}

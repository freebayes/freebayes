#include "Caller.h"

// local helper macro to improve code readability
#define DEBUG_LOG(msg) \
    if (parameters->record) { logFile << msg; } \
    if (parameters->debug) { cerr << msg; }


using namespace std;


// XXX TODO change these void functions to bool

// open BAM input file
void Caller::openBam(void) {

    DEBUG_LOG("Opening BAM fomat alignment input file: " << parameters->bam << " ...")
  
    bamReader.Open(parameters->bam.c_str());

    DEBUG_LOG(" done" << endl);
}

void Caller::openLogFile(void) {

    logFile.open(parameters->log.c_str(), ios::out);
    if (parameters->record) {
        if (parameters->debug) {cerr << "Opening log file: " << parameters->log << " ...";}

        if (!logFile) {
            cerr << " unable to open file: " << parameters->log << endl;
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
            cerr << "unable to open file: " << parameters->samples << endl;
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
                if (parameters->debug) cerr << "found sample " << s << endl;
                sampleList.push_back(s);
            }
        }
    } else { // no samples file given, read from BAM file header for sample names
        // retrieve header information
        if (parameters->debug) cerr << "no sample list file given, attempting to read sample names from bam file" << endl;

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
                if (parameters->debug) cerr << "found sample " << name << endl;
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

    DEBUG_LOG("Number of ref seqs: " << bamReader.GetReferenceCount() << endl);

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

    DEBUG_LOG("processing fasta reference " << parameters->fasta);

    //--------------------------------------------------------------------------
    // process input fasta file
    //--------------------------------------------------------------------------
    // This call loads the reference and reads any index file it can find.
    // If it can't find an index file for the reference, it will attempt to
    // generate one alongside it.

    reference = new FastaReference(parameters->fasta);

    referenceSequenceCount = 0;

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
        referenceSequenceNames.push_back(entry.name);
        referenceSequenceNameToID[entry.name] = id++;
        ++referenceSequenceCount;

    }

    DEBUG_LOG(" done." << endl);

}

void Caller::loadReferenceSequence(int seqID) {
    currentSequence = reference->getSequence(referenceSequenceNames[seqID]);
}

void Caller::loadTargetRegions(void) {

  // process target region input file
  
  // if target file specified use targets from file
  int targetCount = 0;
  if (parameters->targets != "") {
    
    //------------------------------------------------------------------------
    // open input BED file if required
    //------------------------------------------------------------------------
    
    // report
    DEBUG_LOG("Making BedReader object for target file: " << parameters->targets << " ...");
    
    // make
    BedReader bedReader(parameters->targets);
    
    if (! bedReader.isOpen()) {
      cerr << "Unable to open target file: " << parameters->targets << "... terminating." << endl;
      exit(1);
    }
    
    //------------------------------------------------------------------------
    // iterate through entries
    //------------------------------------------------------------------------
    BedData bd;
    while (bedReader.getNextEntry(bd)) {
        if (parameters->debug2) {
            cerr << bd.seq << "\t" << bd.left << "\t" << bd.right << "\t" << bd.desc << endl;
        }
        if (bd.left < 1 || bd.right > referenceSequences[referenceSequenceNameToID[bd.seq]].RefLength) {
          cerr << "Target region coordinate outside of reference sequence bounds... terminating." << endl;
          exit(1);
        }
        if (targetsByRefseq.count(bd.seq) > 0) {
            targetsByRefseq[bd.seq].push_back(bd);
        } else {
            vector<BedData> bv;
            bv.push_back(bd);
            targetsByRefseq[bd.seq] = bv;
        }
        targetCount++;
    }
    
    //------------------------------------------------------------------------
    // close
    //------------------------------------------------------------------------
    bedReader.close();

    DEBUG_LOG("done" << endl);

  }

  // otherwise analyze all reference sequences from BAM file
  else {
    RefVector::iterator refIter = referenceSequences.begin();
    RefVector::iterator refEnd  = referenceSequences.end();
    for( ; refIter != refEnd; ++refIter) {
      RefData refData = *refIter;
      string refName = refData.RefName;
      BedData bd;
      bd.seq = refName;
      bd.left = 1;
      bd.right = refData.RefLength;
      targetsByRefseq[bd.seq].push_back(bd);
      targetCount++;
    }
  }

  DEBUG_LOG("Number of target regions: " << targetCount << endl);

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
  DEBUG_LOG("opening report output file for writing: " << parameters->rpt << "...");

  // open output streams
  bool outputRPT, outputVCF; // for legibility

  if (parameters->rpt != "") {
      outputRPT = true;
      rptFile.open(parameters->rpt.c_str());
      if (!rptFile) {
        DEBUG_LOG(" unable to open file: " << parameters->rpt << endl);
        exit(1);
      }
  } else { outputRPT = false; }

  if (parameters->vcf != "") {
      outputVCF = true;
      vcfFile.open(parameters->vcf.c_str());
      if (!vcfFile) {
        DEBUG_LOG(" unable to open file: " << parameters->vcf << endl);
        exit(1);
      }
  } else { outputVCF = false; }
  DEBUG_LOG(" done." << endl);

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
      rptFile << "#   --parameters->BQL2 = " << parameters->BQL2 << endl;
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

    // initialization
    // NOTE: these void functions have side effects, and currently have to be called in this order
    // this separation is done to improve legibility and debugging
    // perhaps it will just increase confusion
    openLogFile();
    openBam();
    getSampleNames();
    loadFastaReference();
    loadBamReferenceSequenceNames();
    loadTargetRegions();
    initializeOutputFiles();

    currentRefID = -1; // will get set properly via toNextRefID
    toNextRefID(); // initializes currentRefID
    toFirstTargetPosition(); // initializes currentTarget, currentAlignment
}


RegisteredAlignment& Caller::registerAlignment(BamAlignment& alignment) {

    RegisteredAlignment ra = RegisteredAlignment(alignment); // result

    string rDna = alignment.QueryBases;
    string rQual = alignment.Qualities;
    int rp = 1;  // read position
    int sp = alignment.Position + 1;  // sequence position
              //   ^^^ conversion between 0 and 1 based index

    string sampleName;
    if (!alignment.GetReadGroup(sampleName)) {
        cerr << "ERROR: Couldn't find read group id for BAM Alignment " << alignment.Name << endl;
    }

    vector<CigarOp>::const_iterator cigarIter = alignment.CigarData.begin();
    vector<CigarOp>::const_iterator cigarEnd  = alignment.CigarData.end();
    for ( ; cigarIter != cigarEnd; ++cigarIter ) {
        unsigned int l = (*cigarIter).Length;
        char t = (*cigarIter).Type;
      
        if (t == 'S') { // soft clip
            rp += l;
        } else if (t == 'M') { // match or mismatch
            for (int i=1; i<=l; i++) {
      
                // extract aligned base
                string b;
                TRY { b = rDna.substr(rp-1, 1); } CATCH;

                // convert base quality value into short int
                short qual = qualityChar2ShortInt(rQual[rp-1]);

                // get reference allele
                string sb;
                TRY { sb = currentSequence.substr(sp-1, 1); } CATCH;

                // register mismatch
                if (b != sb && qual >= parameters->BQL2)
                    ++ra.mismatches;
      
                // update positions
                ++sp;
                ++rp;
            }
            // XXX what about 'N' s?
        } else if (t == 'D') { // deletion

            // extract base quality of left and right flanking non-deleted bases
            short qL = qualityChar2ShortInt(rQual[rp-1]);
            short qR = qualityChar2ShortInt(rQual[rp]);

            // calculate maximum of the two qualities values
            short qual = min(qL, qR); // XXX was max, but min makes more sense, right ?
            if (qual >= parameters->BQL2) {
                //Allele::Allele(AlleleType, ReferenceID, Position, int, std::string, SampleID, Strand, short int)
                ra.alleles.push_back(Allele(ALLELE_DELETION, currentRefID, sp, l, "", sampleName,
                      (!alignment.IsReverseStrand()) ? STRAND_FORWARD : STRAND_REVERSE, qual));
            }

            sp += l;  // update sample position

        } else if (t == 'I') { // insertion

            vector<short> quals;
            for (int i=1; i<=l; i++) {

                // extract base quality of inserted base
                quals.push_back(qualityChar2ShortInt(rQual[rp-1]));

                rp += 1; // update read position
            }        

            // calculate joint quality, which is the probability that there are no errors in the observed bases
            short qual = jointQuality(quals);
            // register insertion + base quality with reference sequence
            if (qual >= parameters->BQL2) // XXX this cutoff may not make sense for long indels... the joint quality is much lower than the 'average' quality
                ra.alleles.push_back(Allele(ALLELE_INSERTION, currentRefID, sp, l, rDna.substr(rp, l), sampleName,
                                        (!alignment.IsReverseStrand()) ? STRAND_FORWARD : STRAND_REVERSE, qual));

        } // not handled, skipped region 'N's
    } // end cigar iter loop
}

void Caller::updateAlignmentQueue(void) {

    // BamAlignments are 0-based ... does this jive?
    BamAlignment* alignment = &registeredAlignmentQueue.back().alignment;
    // pop from the back until we get to an alignment that overlaps our current position
    while (!(currentTarget->left < alignment->Position <= currentTarget->right ||
            currentTarget->left < alignment->Position + alignment->Length <= currentTarget->right)) {
        registeredAlignmentQueue.pop_back();
        alignment = &registeredAlignmentQueue.back().alignment;
    }
    // push to the front until we get to an alignment that doesn't overlap our
    // current position or we reach the end of available alignments
    // filter input reads; only allow mapped reads with a certain quality
    bool moreAlignments = true; // flag to catch BAM EOF
    while (moreAlignments && currentTarget->left < currentAlignment.Position <= currentTarget->right ||
            currentTarget->left < currentAlignment.Position + currentAlignment.Length <= currentTarget->right) {
        // only process if mapped
        if (currentAlignment.IsMapped()) {
            RegisteredAlignment& ra = registerAlignment(currentAlignment);
            // filters
            // 'read mask' --- this just means "don't consider snps right next to potential indels
            //                 ... but it should be implemented
            // low mapping quality --- hmmm ... could calculate it with the jointQuality function?
            // duplicates --- tracked via each BamAlignment
            if (!(ra.mismatches > parameters->RMU))
                registeredAlignmentQueue.push_front(ra);
            // TODO collect statistics here...
        }
        moreAlignments &= bamReader.GetNextAlignment(currentAlignment);
    }
}

vector<BedData> Caller::targetsInCurrentRefSeq(void) {
    return targetsByRefseq[referenceSequenceNames[currentRefID]];
}

bool Caller::toNextRefID(void) {
    while (targetsInCurrentRefSeq().size() == 0 && currentRefID < referenceSequenceNames.size()) {
        ++currentRefID;
    }
    if (currentRefID >= referenceSequenceNames.size())
        return false;
    else
        return true;
}

// initialization function, should only be called via constructor
bool Caller::toFirstTargetPosition(void) {
    currentTarget = &targetsInCurrentRefSeq().front();
    currentPosition = currentTarget->left;
    bamReader.Jump(currentRefID, currentPosition);
    if (!bamReader.GetNextAlignment(currentAlignment)) {
        // probably indicates an error in the bam file, as this should be our first pass
        cerr << "Bam file has no alignments??" << endl;
        exit(1);
    }
    loadReferenceSequence(currentRefID);
    updateAlignmentQueue();
    DEBUG_LOG("  Processing target: " << currentTarget->seq << ":"
            << currentTarget->left << "-" << currentTarget->right <<
            endl);
    return true;
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
bool Caller::toNextTargetPosition(void) {

    if (currentPosition + 1 > currentTarget->right) {  // time to move targets
        // if we are at the end of the targets in the current refseq
        if (currentTarget == &targetsInCurrentRefSeq().back()) {  
            // if there are more reference sequences to process
            if (currentRefID + 1 < referenceSequenceCount) {
                do {
                    ++currentRefID;
                } while (currentRefID + 1 < referenceSequenceCount && targetsInCurrentRefSeq().size() == 0);
                if (currentRefID + 1 == referenceSequenceCount && targetsInCurrentRefSeq().size() == 0) {
                    return false; // ... done
                } // otherwise we have remaining targets
                currentTarget = &targetsInCurrentRefSeq().front();
                currentPosition = currentTarget->left;
                if (!bamReader.Jump(currentRefID, currentPosition)) {
                    cerr << "ERROR: cannot jump to refseq:" << currentRefID << ", position:" << currentPosition << " in bam file!" << endl;
                    exit(1);
                }
                if (!bamReader.GetNextAlignment(currentAlignment))
                    return false; // done
                loadReferenceSequence(currentRefID);
                updateAlignmentQueue();
                DEBUG_LOG("  Processing target: " << currentTarget->seq << ":"
                        << currentTarget->left << "-" << currentTarget->right <<
                        endl);
                return true;
            // if none, we are ...
            } else {
                return false; // ... done
            }
        }
    } else { // if we are still in the current target, just step the position
        ++currentPosition;
        updateAlignmentQueue();
    }
}

bool Caller::getNextAlleles(vector<Allele>& alleles) {
    bool more = toNextTargetPosition();
    if (more) {
        getAlleles(alleles);
        return true;
    } else {
        return false;
    }
}

// updates the passed vector with the current alleles at the caller's target position
void Caller::getAlleles(vector<Allele>& alleles) {

    // clear the allele vector
    alleles.clear();

    // get the alleles overlapping the current position
    // 
    // NB: if we ignore the differences between deletions and insertions, we
    // will report them at a number of positions, as several positions in the
    // reference can correspond to a deletion.  and we will report insertions
    // only at one position.  for now this is fine; it can be fixed when the
    // i/o systems are verified to be working.
    
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

#include "Caller.h"

// local helper macro to improve code readability
#define DEBUG_LOG(msg) \
    if (parameters.record) { logFile << msg; } \
    if (parameters.debug) { cerr << msg; }


using namespace std;


// XXX TODO change these void functions to bool

// open BAM input file
void BayesCaller::openBam(void) {

    DEBUG_LOG("Opening BAM fomat alignment input file: " << bam << " ...")
  
    bamReader.Open(parameters.bam.c_str());

    DEBUG_LOG(" done" << endl);
}

void BayesCaller::openLogFile(void) {

    logFile.open(log.c_str(), ios::out);
    if (parameters.record) {
        if (parameters.debug) {cerr << "Opening log file: " << log << " ...";}

        if (!logFile) {
            cerr << " unable to open file: " << log << endl;
            exit(1);
        }
        if (parameters.debug) {cerr << " done." << endl;}
    }

    if (parameters.record) {
        logFile << "Command line:";
        for (int i=0; i<argc; i++) {
          logFile << " " << argv[i];
        }
        logFile << endl << *this;
    }

    if (parameters.debug) {
        cerr << "Command line:";
        for (int i=0; i<argc; i++) {
          cerr << " " << argv[i];
        }
        cerr << endl << *this;
    }
}

// read sample list file or get sample names from bam file header
void BayesCaller::getSampleNames(void) {
    // If a sample file is given, use it.  But otherwise process the bam file
    // header to get the sample names.
    //
    if (parameters.samples != "") {
        ifstream sampleFile(samples.c_str(), ios::in);
        if (! sampleFile) {
            cerr << "unable to open file: " << samples << endl;
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
                if (parameters.debug) cerr << "found sample " << s << endl;
                sampleList.push_back(s);
            }
        }
    } else { // no samples file given, read from BAM file header for sample names
        // retrieve header information
        if (parameters.debug) cerr << "no sample list file given, attempting to read sample names from bam file" << endl;

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
                if (parameters.debug) cerr << "found sample " << name << endl;
                sampleList.push_back(name);
            }
        }
    }
}

void BayesCaller::loadBamReferenceSequenceNames(void) {

    //--------------------------------------------------------------------------
    // read reference sequences from input file
    //--------------------------------------------------------------------------

    // store the names of all the reference sequences in the BAM file
    referenceSequences = bamReader.GetReferenceData();

    DEBUG_LOG("Number of ref seqs: " << bamReader.GetReferenceCount() << endl);

}


void BayesCaller::loadFastaReference(void) {

    // XXX we don't technically need to load the whole thing into memory
    // the FastaReference class will let us grab individual sequences and subsequences 
    // very fast from the file
    // thus cutting down on memory requirements...
    //
    // perhaps a good approach is to open the file here
    // and then get long subsequences at runtime
    // 
    // this keeps our memory requirements low, and will allow us to operate unmodified on more systems

    DEBUG_LOG("processing fasta reference " << parameters.fasta);

    //--------------------------------------------------------------------------
    // process input fasta file
    //--------------------------------------------------------------------------
    // This call loads the reference and reads any index file it can find.
    // If it can't find an index file for the reference, it will attempt to
    // generate one alongside it.

    reference = FastaReference(parameters.fasta);

    referenceSequenceCount = 0;

    int id = 0;

    //--------------------------------------------------------------------------
    // load ref seq names into hash
    //--------------------------------------------------------------------------
    for(map<string, FastaIndexEntry>::const_iterator it = fastaReference->index->begin(); 
          it != fastaReference->index->end(); ++it) {

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

void BayesCaller::loadReferenceSequence(int seqID) {
    currentSequence = reference.getSequence(referenceSequenceNameToID[seqID]);
}

void BayesCaller::loadTargetRegions(void) {

  // process target region input file
  
  // if target file specified use targets from file
  int targetCount = 0;
  if (parameters.targets != "") {
    
    //------------------------------------------------------------------------
    // open input BED file if required
    //------------------------------------------------------------------------
    
    // report
    DEBUG_LOG("Making BedReader object for target file: " << parameters.targets << " ...");
    
    // make
    bedReader = BedReader(parameters.targets);
    
    if (! bedReader.isOpen()) {
      cerr << "Unable to open target file: " << parameters.targets << "... terminating." << endl;
      exit(1);
    }
    
    //------------------------------------------------------------------------
    // iterate through entries
    //------------------------------------------------------------------------
    BedData bd;
    while (bedReader.getNextEntry(bd)) {
      if (parameters.debug2) {
        cerr << bd.seq << "\t" << bd.left << "\t" << bd.right << "\t" << bd.desc << endl;
      }
      if (referenceTargetMap.count(bd.seq) > 0) {
        if (bd.left < 1 || bd.right > referenceSequences[referenceSequenceNameToID[bd.seq]].RefLength) {
          cerr << "Target region coordinate outside of reference sequence bounds... terminating." << endl;
          exit(1);
        }
        targetsByRefseq[bd.seq].push_back(bd);
        targetCount++;
      }
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
      referenceTargetMap[bd.seq].push_back(bd);
      targetCount++;
    }
  }

  DEBUG_LOG("Number of target regions: " << targetCount << endl);

}

void BayesCaller::initializeOutputFiles(void) {

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // open output file(s)
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  // open report output file
  //----------------------------------------------------------------------------

  // report
  DEBUG_LOG("opening report output file for writing: " << rpt << "...");

  // open output streams
  bool outputRPT, outputVCF; // for legibility

  if (parameters.rpt != "") {
      outputRPT = true;
      rptFile.open(rpt.c_str());
      if (!rptFile) {
        DEBUG_LOG(" unable to open file: " << rpt << endl);
        exit(1);
      }
  } else { outputRPT = false; }

  if (parameters.vcf != "") {
      outputVCF = true;
      vcfFile.open(vcf.c_str());
      if (!vcfFile) {
        DEBUG_LOG(" unable to open file: " << rpt << endl);
        exit(1);
      }
  } else { outputVCF = false; }
  DEBUG_LOG(" done." << endl);

  //----------------------------------------------------------------------------
  // write header information
  //----------------------------------------------------------------------------
  if (outputRPT) {
      rptFile << "# Command line that generated this output:";
      for (int i=0; i<argc; i++) {
        rptFile << " " << argv[i];
      }
      rptFile << endl;
      rptFile << "#" << endl;
      rptFile << "# Complete list of parameter values:" << endl;
      rptFile << "#   --bam = " << bam << endl;
      rptFile << "#   --fasta = " << fasta << endl;
      rptFile << "#   --targets = " << targets << endl;
      rptFile << "#   --samples = " << samples << endl;
      rptFile << "#   --rpt = " << rpt << endl;
      rptFile << "#   --log = " << log << endl;
      rptFile << "#   --useRefAllele = " <<  bool2String[useRefAllele] << endl;
      rptFile << "#   --forceRefAllele = " <<  bool2String[forceRefAllele] << endl;
      rptFile << "#   --MQR = " << MQR << endl;
      rptFile << "#   --BQR = " << BQR << endl;
      rptFile << "#   --ploidy = " << ploidy << endl;
      rptFile << "#   --sampleNaming = " << sampleNaming << endl;
      rptFile << "#   --sampleDel = " << sampleDel << endl;
      rptFile << "#   --BQL0 = " << BQL0 << endl;
      rptFile << "#   --MQL0 = " << MQL0 << endl;
      rptFile << "#   --BQL1 = " << BQL1 << endl;
      rptFile << "#   --MQL1 = " << MQL1 << endl;
      rptFile << "#   --parameters.BQL2 = " << parameters.BQL2 << endl;
      rptFile << "#   --RMU = " << RMU << endl;
      rptFile << "#   --IDW = " << IDW << endl;
      rptFile << "#   --TH = " << TH << endl;
      rptFile << "#   --PVL = " << PVL << endl;
      rptFile << "#   --algorithm = " << algorithm << endl;
      rptFile << "#   --RDF = " << RDF << endl;
      rptFile << "#   --WB = " << WB << endl;
      rptFile << "#   --TB = " << TB << endl;
      rptFile << "#   --includeMonoB = " <<  bool2String[includeMonoB] << endl;
      rptFile << "#   --TR = " << TR << endl;
      rptFile << "#   --I = " << I << endl;
      rptFile << "#   --debug = " <<  bool2String[debug] << endl;
      rptFile << "#   --debug2 = " <<  bool2String[debug2] << endl;
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
BayesCaller::BayesCaller(Parameters& params) {

    // init program parameters
    parameters = params;

    // set up position tracking
    currentTarget = NULL;  // hmmm

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

}


RegisteredAlignment& BayesCaller::registerAlignment(BamAlignment& alignment) {

    RegisteredAlignment ra = RegisteredAlignment(alignment); // result

    string rDna = alignment.QueryBases;
    string rQual = alignment.Qualities;
    int rp = 1;  // read position
    int sp = alignment.Position + 1;  // sequence position
              //   ^^^ conversion between 0 and 1 based index

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
                if (b != sb && qual >= parameters.BQL2)
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
            if (qual >= parameters.BQL2)
                ra.alleles.push_back(Allele(ALLELE_DELETION, currentRefID, sp, l, '',
                      (!alignment.IsReverseStrand()) ? STRAND_FORWARD : STRAND_REVERSE, qual));

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
            if (qual >= parameters.BQL2) // XXX this cutoff may not make sense for long indels... the joint quality is much lower than the 'average' quality
                ra.alleles.push_back(Allele(ALLELE_INSERTION, currentRefID, sp, l, rDna.substr(rp, l),
                                        (!alignment.IsReverseStrand()) ? STRAND_FORWARD : STRAND_REVERSE));

        } // not handled, skipped region 'N's
    } // end cigar iter loop
}

void BayesCaller::updateAlignmentQueue(void) {

    // BamAlignments are 0-based ... does this jive?
    BamAlignment* firstAlignment = registeredAlignmentQueue.front().alignment;
    // pop from the back until we get to an alignment that overlaps our current position
    while (!(currentTarget->left < front->Position <= currentTarget->right ||
            currentTarget->left < front->Position + front->Length <= currentTarget->right)) {
        registeredAlignmentQueue.pop();
        firstAlignment = registeredAlignmentQueue.front().alignment;
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
            // low mapping quality --- hmmm ... could calculate it with the jointQuality function?
            // duplicates --- handled via BamAlignment
            if (!(ra.mismatches > parameters.RMU))
                registeredAlignmentQueue.push(ra);
            // TODO collect statistics here...
        }
        moreAlignments &= bamReader.GetNextAlignment(currentAlignment);
    }
}

vector<BedData>* BayesCaller::targetsInCurrentRefSeq(void) {
    return targetsByRefseq[referenceSequenceNames[currentRefID]];
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
bool BayesCaller::toNextTargetPosition(void) {

    // first pass
    if (currentTarget == NULL) {
        currentTarget = targetsInCurrentRefSeq().begin();
        currentRefID = referenceSequenceNameToID[currentTarget->seq];
        currentPosition = currentTarget->left;
        bamReader.Jump(currentRefID, currentPosition);
        if (!bamReader.GetNextAlignment(currentAlignment)) {
            // probably indicates an error in the bam file, as this should be our first pass
            cerr << "Bam file has no alignments ???" << endl;
            exit(1);
        }
        loadReferenceSequence(currentRefID);
        updateAlignmentQueue();
        DEBUG_LOG("  Processing target: " << currentTarget.seq << ":"
                << currentTarget.left << "-" << currentTarget.right <<
                endl);
        return true;
    }

    if (currentPosition + 1 > currentTarget.right) {  // time to move targets
        // if we are at the end of the targets in the current refseq
        if (currentTarget == targetsInCurrentRefSeq().end()) {  
            // if there are more reference sequences to process
            if (currentRefID + 1 < referenceSequenceCount) {
                do {
                    ++currentRefID;
                } while (currentRefID + 1 < referenceSequenceCount && targetsInCurrentRefSeq().size() == 0);
                if (currentRefID + 1 == referenceSequenceCount && targetsInCurrentRefSeq().size() == 0) {
                    return false; // ... done
                } // otherwise we have remaining targets
                currentTarget = targetsInCurrentRefSeq().begin();
                currentPosition = currentTarget->left;
                if (!bamReader.Jump(currentRefID, currentPosition)) {
                    cerr << "ERROR: cannot jump to refseq:" << currentRefID << ", position:" << currentPosition << " in bam file!" << endl;
                    exit(1);
                }
                if (!bamReader.GetNextAlignment(currentAlignment))
                    return false; // done
                loadReferenceSequence(currentRefID);
                updateAlignmentQueue();
                DEBUG_LOG("  Processing target: " << currentTarget.seq << ":"
                        << currentTarget.left << "-" << currentTarget.right <<
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

bool BayesCaller::getNextAlleles(vector<Allele>& alleles) {
    return (toNextTargetPosition() && getAlleles(alleles));
}

/*
 *
 * plan:
 *
 * presently basecalls are stored in a big map->map->map thing
 * for the whole target region under consideration:
 *
 *      pos      sample        basecalls
 * map <int, map<string, vector<Basecall>>
 *
 * but this isn't really necessary because we only consider the putative
 * alleles at a single given position at a time.
 *
 * so alternatively:
 *
 * vector <int, vector<int, vector<Basecall>>
 *
 * indexed access...
 *
 * the outside int is probably unnecessary; as we're just reporting putative alleles
 * for each position we then get:
 *
 *         sample     alleles
 * vector <int, vector<Allele>>
 *
 */

// updates the passed vector with the current alleles at the caller's target position
bool BayesCaller::getAlleles(vector<Allele>& alleles) {

    // clear the allele vector
    alleles.clear();

    // get the alleles overlapping the current position
    // what this means is unclear
    // here are the problems
    // with deletions, we have bases in the reference that correspond to the deletion
    // but with insertions, we only have one reference base that corresponds
    // ... i suppose this is OK
    //

    // loop over allele queue items
    // find alleles that 

}




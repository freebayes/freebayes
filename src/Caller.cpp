#include "Caller.h"
#include "multichoose.h" // includes generic functions, so it must be included here
                         // otherwise we will get a linker error
                         // see: http://stackoverflow.com/questions/36039/templates-spread-across-multiple-files
                         // http://www.cplusplus.com/doc/tutorial/templates/ "Templates and Multi-file projects"
#include "multipermute.h"

// local helper debugging macros to improve code readability
#define LOG(msg) \
    if (parameters.record) { logFile << msg << endl; } \
    if (parameters.debug) { cerr << msg << endl; }

// lower-priority messages
#define LOG2(msg) \
    if (parameters.record) { logFile << msg << endl; } \
    if (parameters.debug2) { cerr << msg << endl; }

// must-see error messages
#define ERROR(msg) \
    if (parameters.record) { logFile << msg << endl; } \
    cerr << msg << endl;

using namespace std;


// XXX TODO change these void functions to bool

// open BAM input file
void Caller::openBams(void) {

    if (parameters.bams.size() == 1) {
        LOG("Opening BAM fomat alignment input file: " << parameters.bams.front() << " ...");
    } else if (parameters.bams.size() > 1) {
        LOG("Opening " << parameters.bams.size() << " BAM fomat alignment input files");
        for (vector<string>::const_iterator b = parameters.bams.begin(); 
                b != parameters.bams.end(); ++b) {
            LOG2(*b);
        }
    }
    bamMultiReader.Open(parameters.bams);
    LOG(" done");

}

void Caller::openLogFile(void) {

    logFile.open(parameters.log.c_str(), ios::out);
    if (parameters.record) {
        if (parameters.debug) {cerr << "Opening log file: " << parameters.log << " ...";}

        if (!logFile) {
            ERROR(" unable to open file: " << parameters.log);
            exit(1);
        }
        if (parameters.debug) {cerr << " done." << endl;}
    }
    // previously we wrote the command to the logfile here
}

// read sample list file or get sample names from bam file header
void Caller::getSampleNames(void) {
    // If a sample file is given, use it.  But otherwise process the bam file
    // header to get the sample names.
    //

  //--------------------------------------------------------------------------
  //--------------------------------------------------------------------------
  // read sample list file or get sample names from bam file header
  //--------------------------------------------------------------------------
  //--------------------------------------------------------------------------
  //
  // If a sample file is given, use it.  But otherwise process the bam file
  // header to get the sample names.
  //
    if (parameters.samples != "") {
        ifstream sampleFile(parameters.samples.c_str(), ios::in);
        if (! sampleFile) {
            cerr << "unable to open file: " << parameters.samples << endl;
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
                LOG("found sample " << s);
                sampleList.push_back(s);
            }
        }
    } 

    // retrieve header information

    string bamHeader = bamMultiReader.GetHeaderText();

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
            string name = "";
            string readGroupID = "";
            for (vector<string>::const_iterator r = readGroupParts.begin(); r != readGroupParts.end(); ++r) {
                vector<string> nameParts;
                boost::split(nameParts, *r, boost::is_any_of(":"));
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
            LOG("found read group id " << readGroupID << " containing sample " << name);
            sampleListFromBam.push_back(name);
            readGroupToSampleNames[readGroupID] = name;
        }
    }
    //cout << sampleListFromBam.size() << endl;
     // no samples file given, read from BAM file header for sample names
    if (sampleList.size() == 0) {
        LOG("no sample list file given, reading sample names from bam file");
        for (vector<string>::const_iterator s = sampleListFromBam.begin(); s != sampleListFromBam.end(); ++s) {
            LOG("found sample " << *s);
            if (!stringInVector(*s, sampleList)) {
                sampleList.push_back(*s);
            }
        }
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
}

void Caller::loadBamReferenceSequenceNames(void) {

    //--------------------------------------------------------------------------
    // read reference sequences from input file
    //--------------------------------------------------------------------------

    // store the names of all the reference sequences in the BAM file
    referenceSequences = bamMultiReader.GetReferenceData();

    LOG("Number of ref seqs: " << bamMultiReader.GetReferenceCount());

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

    LOG("processing fasta reference " << parameters.fasta);

    //--------------------------------------------------------------------------
    // process input fasta file
    //--------------------------------------------------------------------------
    // This call loads the reference and reads any index file it can find.
    // If it can't find an index file for the reference, it will attempt to
    // generate one alongside it.

    reference = new FastaReference(parameters.fasta);

    fastaReferenceSequenceCount = 0;

    int id = 0;

    //--------------------------------------------------------------------------
    // load ref seq names into hash
    //--------------------------------------------------------------------------
    for(vector<FastaIndexEntry>::const_iterator entry = reference->index->begin(); 
          entry != reference->index->end(); ++entry) {

        // we split out the first word of the full sequence name for keying our sequences
        // as often the full sequence names are overkill
        vector<string> sequenceNameParts;
        boost::split(sequenceNameParts, entry->name, boost::is_any_of(" "));
        string name = sequenceNameParts.front();
        LOG2("sequence name " << name << " id = " << id);

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
    LOG2("loading reference subsequence " << target->seq << " from " << target->left << " - " << before << " to " << target->right << " + " << after << " + before");
    loadReferenceSequence(target->seq, target->left - before, target->right - target->left + after + before);
}

void Caller::loadTargets(void) {

  // process target region input file
  
  // if target file specified use targets from file
  int targetCount = 0;
  // if we have a targets file, use it...
  if (parameters.targets != "") {
    
    LOG("Making BedReader object for target file: " << parameters.targets << " ...");
    
    BedReader bedReader(parameters.targets);
    
    if (! bedReader.isOpen()) {
      ERROR("Unable to open target file: " << parameters.targets << "... terminating.");
      exit(1);
    }
    
    BedData bd;
    while (bedReader.getNextEntry(bd)) {
        if (parameters.debug2) {
            cerr << bd.seq << "\t" << bd.left << "\t" << bd.right << "\t" << bd.desc << endl;
        }
        // TODO add back check that the right bound isn't out of bounds
        string seqName = reference->sequenceNameStartingWith(bd.seq);
        if (bd.left < 1 || bd.right < bd.left || bd.right >= reference->sequenceLength(seqName)) {
          ERROR("Target region coordinates (" << bd.seq << " " << bd.left << " " << bd.right << ") outside of reference sequence bounds (" << bd.seq << " " << reference->sequenceLength(seqName) << ") terminating.");
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

// initialization function
// sets up environment so we can start registering alleles
Caller::Caller(int argc, char** argv) : parameters(Parameters(argc, argv))
{

    // initialization
    // NOTE: these void functions have side effects, and currently have to be called in this order
    // this separation is done to improve legibility and debugging
    // perhaps it will just increase confusion
    openLogFile();
    openBams();
    getSampleNames();
    loadFastaReference();
    loadBamReferenceSequenceNames();
    loadTargets();

    currentRefID = 0; // will get set properly via toNextRefID
    //toNextRefID(); // initializes currentRefID
    //toFirstTargetPosition(); // initializes currentTarget, currentAlignment
    currentTarget = NULL; // to be initialized on first call to getNextAlleles
}

Caller::~Caller(void) {
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
        << ra.alignment.Qualities << endl
        << ra.alleles << endl;
    return out;
}

RegisteredAlignment Caller::registerAlignment(BamAlignment& alignment, string sampleName) {

    RegisteredAlignment ra = RegisteredAlignment(alignment); // result

    string rDna = alignment.QueryBases;
    string rQual = alignment.Qualities;
    int rp = 0;  // read position, 0-based relative to read
    int csp = currentSequencePosition(alignment); // current sequence position, 0-based relative to currentSequence
    int sp = alignment.Position;  // sequence position
    //   ^^^ conversion between 0 and 1 based index

    string readName = alignment.Name;

    LOG2("registering alignment " << rp << " " << csp << " " << sp << endl <<
            "alignment readName " << readName << endl <<
            "alignment sampleID " << sampleName << endl << 
            "alignment position " << alignment.Position << endl <<
            "alignment length " << alignment.Length << endl <<
            "alignment AlignedBases.size() " << alignment.AlignedBases.size() << endl <<
            "alignment end position " << alignment.Position + alignment.AlignedBases.size());

    LOG2(rDna << endl << alignment.AlignedBases << endl << currentSequence.substr(csp, alignment.AlignedBases.size()));

    /*
     * This should be simpler, but isn't because the cigar only records matches
     * for sequences that have embedded mismatches.
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

    vector<CigarOp>::const_iterator cigarIter = alignment.CigarData.begin();
    vector<CigarOp>::const_iterator cigarEnd  = alignment.CigarData.end();
    for ( ; cigarIter != cigarEnd; ++cigarIter ) {
        unsigned int l = (*cigarIter).Length;
        char t = (*cigarIter).Type;
        LOG2("cigar item: " << t << l);

        if (t == 'S') { // soft clip
            rp += l;
        } else if (t == 'M') { // match or mismatch
            int firstMatch = csp; // track the first match after a mismatch, for recording 'reference' alleles
                                        // we start one back because the first position in a match should match...
            for (int i=0; i<l; i++) {

                // extract aligned base
                string b;
                TRY { b = rDna.at(rp); } CATCH;

                // convert base quality value into short int
                short qual = qualityChar2ShortInt(rQual.at(rp));

                // get reference allele
                string sb;
                TRY { sb = currentSequence.at(csp); } CATCH;

                // XXX XXX XXX THIS LOGIC IS ILLOGICAL!!!
                //
                // record mismatch if we have a mismatch here
                if (b != sb) {
                    if (firstMatch < csp) {
                        // TODO ; verify that the read and reference sequences *do* match
                        int length = csp - firstMatch - 1;
                        string matchingSequence = currentSequence.substr(csp - length - 1, length);
                        string readSequence = rDna.substr(rp - length - 1, length);
                        string qualstr = rQual.substr(rp - length - 1, length);
                        Allele* allele = new Allele(ALLELE_REFERENCE,
                                currentTarget->seq, sp - length - 1, length, 
                                matchingSequence, readSequence, sampleName, alignment.Name,
                                !alignment.IsReverseStrand(), -1, qualstr,
                                alignment.MapQuality);
                        ra.alleles.push_back(allele);
                    }
                    firstMatch = csp + 1;
                }

                // register mismatch
                if (b != sb) {
                    if (qual >= parameters.BQL2) {
                        // record 'reference' allele for last matching region
                        ra.mismatches++;
                        Allele* allele = new Allele(ALLELE_SNP, currentTarget->seq, sp, 1, sb, b,
                                sampleName, alignment.Name, !alignment.IsReverseStrand(), qual, "", alignment.MapQuality);
                        ra.alleles.push_back(allele);
                    }
                    firstMatch = csp + 1;
                }

                // update positions
                ++sp;
                ++csp;
                ++rp;
            }
            if (firstMatch < csp) {  // TODO use some param flag to trigger this operation
                // TODO ; verify that the read and reference sequences *do* match
                int length = csp - firstMatch - 1;
                string matchingSequence = currentSequence.substr(csp - length - 1, length);
                string readSequence = rDna.substr(rp - length - 1, length);
                string qualstr = rQual.substr(rp - length - 1, length);
                Allele* allele = new Allele(ALLELE_REFERENCE,
                        currentTarget->seq, sp - length - 1, length, 
                        matchingSequence, readSequence, sampleName, alignment.Name,
                        !alignment.IsReverseStrand(), -1, qualstr,
                        alignment.MapQuality);
                ra.alleles.push_back(allele);
            }
            // XXX what about 'N' s?
        } else if (t == 'D') { // deletion

            // extract base quality of left and right flanking non-deleted bases
            string qualstr = rQual.substr(rp, 2);

            // calculate joint quality of the two flanking qualities
            short qual = jointQuality(qualstr); // XXX was max, but joint makes more sense, right ?
            if (qual >= parameters.BQL2) {
                Allele* allele = new Allele(ALLELE_DELETION,
                        currentTarget->seq, sp, l,
                        currentSequence.substr(csp, l), "", sampleName, alignment.Name,
                        !alignment.IsReverseStrand(), qual, qualstr,
                        alignment.MapQuality);
                ra.alleles.push_back(allele);
            }

            sp += l;  // update sample position
            csp += l;

        } else if (t == 'I') { // insertion

            string qualstr = rQual.substr(rp, l);

            // calculate joint quality, which is the probability that there are no errors in the observed bases
            short qual = jointQuality(qualstr);
            // register insertion + base quality with reference sequence
            // XXX this cutoff may not make sense for long indels... the joint
            // quality is much lower than the 'average' quality
            if (qual >= parameters.BQL2) {
                Allele* allele = new Allele(ALLELE_INSERTION,
                        currentTarget->seq, sp, l, "", rDna.substr(rp, l),
                        sampleName, alignment.Name, !alignment.IsReverseStrand(), qual,
                        qualstr, alignment.MapQuality);
                ra.alleles.push_back(allele);
            }

            rp += l;

        } // not handled, skipped region 'N's
    } // end cigar iter loop
    //cerr << ra << endl;
    return ra;
}


void Caller::updateAlignmentQueue(void) {

    LOG2("updating alignment queue");

    // push to the front until we get to an alignment that doesn't overlap our
    // current position or we reach the end of available alignments
    // filter input reads; only allow mapped reads with a certain quality
    int i = 0;
    if (bamMultiReader.GetNextAlignment(currentAlignment) && currentAlignment.Position <= currentPosition) {
        do {
            // get read group, and map back to a sample name
            string readGroup;
            if (!currentAlignment.GetTag("RG", readGroup)) {
                ERROR("Couldn't find read group id (@RG tag) for BAM Alignment " <<
                        currentAlignment.Name << " at position " << currentPosition
                        << " in sequence " << currentSequence << " EXITING!");
                exit(1);
            }
            // skip this alignment if we are not analyzing the sample it is drawn from
            if (readGroupToSampleNames.find(readGroup) == readGroupToSampleNames.end()) {
                continue;
            }
            // otherwise, get the sample name and register the alignment to generate a sequence of alleles
            string sampleName = readGroupToSampleNames[readGroup];
            if (currentAlignment.IsMapped()) {
                // we have to register the alignment to acquire some information required by filters
                // such as mismatches
                // filter low mapping quality (what happens if MapQuality is not in the file)
                if (currentAlignment.MapQuality > parameters.MQL0) {
                    RegisteredAlignment ra = registerAlignment(currentAlignment, sampleName);
                    // TODO filters to implement:
                    // duplicates --- tracked via each BamAlignment
                    if (ra.mismatches <= parameters.RMU) {
                        registeredAlignmentQueue.push_front(ra);
                        for (vector<Allele*>::const_iterator allele = ra.alleles.begin(); allele != ra.alleles.end(); ++allele) {
                            registeredAlleles.push_back(*allele);
                        }
                    }
                }
                // TODO collect statistics here...
            }
        } while (bamMultiReader.GetNextAlignment(currentAlignment) && currentAlignment.Position <= currentPosition);
    }

    LOG2("... finished pushing new alignments");

    // pop from the back until we get to an alignment that overlaps our current position
    if (registeredAlignmentQueue.size() > 0) {
        BamAlignment* alignment = &registeredAlignmentQueue.back().alignment;
        // is indexing (0 or 1 based) oK here?
        while (currentPosition >= alignment->Position + alignment->Length && registeredAlignmentQueue.size() > 0) {
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

void Caller::updateRegisteredAlleles(void) {

    // remove reference alleles which are no longer overlapping the current position
    // http://stackoverflow.com/questions/347441/erasing-elements-from-a-vector
    /*
       for (list<Allele*>::iterator allele = registeredAlleles.begin(); allele != registeredAlleles.end(); ) {
       if (currentPosition > (*allele)->position + (*allele)->length) {
       delete *allele;
       allele = registeredAlleles.erase(allele);
       } else {
       ++allele;
       }
       }
       */
    vector<Allele*>& alleles = registeredAlleles;
    for (vector<Allele*>::iterator allele = alleles.begin(); allele != alleles.end(); ++allele) {
        if (currentPosition >= (*allele)->position + (*allele)->length) {
            delete *allele;
            *allele = NULL;
        }
    }
    alleles.erase(remove(alleles.begin(), alleles.end(), (Allele*)NULL), alleles.end());

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

    bool r = bamMultiReader.Jump(refSeqID, currentTarget->left);
    r &= bamMultiReader.GetNextAlignment(alignment);
    int left_gap = currentTarget->left - alignment.Position;

    r &= bamMultiReader.Jump(refSeqID, currentTarget->right - 1);

    int maxPos = 0;
    do {
        r &= bamMultiReader.GetNextAlignment(alignment);
        int newPos = alignment.Position + alignment.AlignedBases.size();
        maxPos = (newPos > maxPos) ? newPos : maxPos;
    } while (alignment.Position <= currentTarget->right);
    // NB ^^^ if all bed files were 0-based start 1-based end, then the above could be
    //        alignment.Position < currentTarget->right
    //        However, with the 'Marth Lab' BED format and other potentially variant formats this may not be the case.
    //        Reading in slightly more data here is not a serious problem, so this quirk allows it.

    int right_gap = maxPos - currentTarget->right;
    //cerr << "left_gap " << left_gap << " right gap " << right_gap << endl;
    loadReferenceSequence(currentTarget,
            (left_gap > 0) ? left_gap : 0,
            (right_gap > 0) ? right_gap : 0);

    LOG2("setting new position " << currentTarget->left);
    currentPosition = currentTarget->left; // bed targets are always 0-based at the left

    LOG2("jumping to first alignment in new target");
    r &= bamMultiReader.Jump(refSeqID, currentTarget->left);
    r &= bamMultiReader.GetNextAlignment(currentAlignment);

    LOG2("clearing registered alignments and alleles");
    registeredAlignmentQueue.clear();
    registeredAlleles.clear();

    return r;

}

// stepping
//
// if the next position is outside of target region
// seek to next target which is in-bounds for its sequence
// if none exist, return false
bool Caller::toNextTargetPosition(void) {

    if (currentTarget == NULL) {
        toFirstTargetPosition();
    } else {
        ++currentPosition;
    }
    if (currentPosition >= currentTarget->right) { // time to move to a new target
        LOG2("next position " << currentPosition <<  " outside of current target right bound " << currentTarget->right);
        if (!toNextTarget()) {
            LOG("no more valid targets, finishing");
            return false;
        }
    }
    LOG2("processing position " << currentPosition << " in sequence " << currentTarget->seq);
    updateAlignmentQueue();
    LOG2("updating registered alleles");
    updateRegisteredAlleles();
    return true;
}

bool Caller::getNextAlleles(list<Allele*>& alleles) {
    if (toNextTargetPosition()) {
        getAlleles(alleles);
        return true;
    } else {
        return false;
    }
}

// updates the passed vector with the current alleles at the caller's target position
void Caller::getAlleles(list<Allele*>& alleles) {

    // this is inefficient but no method besides shared_ptr exists to guarantee
    // that we don't end up with corrupted alleles in this structure-- we won't
    // be able to remove them if we recycle them elsewhere
    alleles.clear();

    // get the variant alleles *at* the current position
    // and the reference alleles *overlapping* the current position
    for (vector<Allele*>::const_iterator a = registeredAlleles.begin(); a != registeredAlleles.end(); ++a) {
        Allele* allele = *a;
        if (((allele->type == ALLELE_REFERENCE 
                 && currentPosition >= allele->position 
                 && currentPosition < allele->position + allele->length) // 0-based, means position + length - 1 is the last included base
                || (allele->position == currentPosition)) 
                && allele->Quality(currentPosition) >= parameters.BQL0
                )
            alleles.push_back(allele);
    }

    alleles.sort();

    // TODO allele sorting by sample on registration
    // for another potential perf boost
    // as we always sort them by sample later
}



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

vector<pair<Genotype, long double> > Caller::probObservedAllelesGivenPossibleGenotypes(vector<Allele*> &observedAlleles, int ploidy) {

    vector<vector<Allele*> > alleleGroups = groupAlleles(observedAlleles, allelesEquivalent);
    vector<tuple<long double, long double, vector<Allele*> > > alleleGroupsAndQualities;
    calculateAlleleGroupProbabilities(alleleGroups, alleleGroupsAndQualities);
    
    vector<Allele> genotypeAlleles =genotypeAllelesFromAlleleGroups(alleleGroups);
    vector<Genotype> genotypes = multichoose(ploidy, genotypeAlleles);

    vector<pair<Genotype, long double> > results;

    for (vector<Genotype>::iterator g = genotypes.begin(); g != genotypes.end(); ++g) {
        long double prob = 0;
        results.push_back(make_pair(*g, probObservedAllelesGivenGenotype(alleleGroupsAndQualities, *g)));
    }

    return results;

}

/*
vector<pair<Genotype, long double> > 
Caller::probObservedAllelesGivenGenotypes(
        vector<Allele*> &observedAlleles,
        vector<Genotype> &genotypes) {

    int ploidy = genotypes.front().size(); // ploidy is determined by the number of alleles in the genotypes
    vector<pair<Genotype, long double> > results;

    vector<vector<Allele*> > alleleGroups = groupAlleles(observedAlleles, allelesEquivalent);
    vector<tuple<long double, long double, vector<Allele*> > > alleleGroupsAndQualities;
    calculateAlleleGroupProbabilities(alleleGroups, alleleGroupsAndQualities); // caches for faster calc

    for (vector<Genotype>::iterator g = genotypes.begin(); g != genotypes.end(); ++g) {
        long double prob = 0;
        results.push_back(make_pair(*g, probObservedAllelesGivenGenotype(alleleGroupsAndQualities, *g)));
    }

    return results;

}
*/

// caches in and out probabilities per allele grouping in alleleGroups
// results are pushed into alleleGroupsAndQualities
// position indicates position at which qualities are to be calculated
void Caller::calculateAlleleGroupProbabilities(vector<vector<Allele*> >& alleleGroups, 
        vector<tuple<long double, long double, vector<Allele*> > >& alleleGroupsAndQualities) {

    for (vector<vector<Allele*> >::iterator group = alleleGroups.begin(); group != alleleGroups.end(); ++group) {
        long double inProb = 0;
        long double outProb = 0;
        for (vector<Allele*>::iterator obs = group->begin(); obs != group->end(); obs++) {
            int q = (*obs)->Quality(currentPosition);
            long double p = phred2float(q);
            if (p == 1) {
                ERROR("Encountered phred Q = 0 when processing allele " << *obs << " at " << currentPosition);
            }
            inProb += log(1 - p);
            outProb += phred2ln(q);
        }
        alleleGroupsAndQualities.push_back(make_tuple(inProb, outProb, *group));
    }

}

vector<pair<Genotype, long double> > 
Caller::probObservedAllelesGivenGenotypes(
        vector<Allele*> &observedAlleles,
        vector<Genotype> &genotypes) {

    vector<pair<Genotype, long double> > results;

    vector<tuple<long double, long double, Allele* > > observedAllelesAndProbs;
    calculateAlleleBinaryProbabilities(observedAlleles, observedAllelesAndProbs); // in and out probability in log space for each allele at this position

    for (vector<Genotype>::iterator g = genotypes.begin(); g != genotypes.end(); ++g) {
        // XXX this following statement is where you should target pruning optimizations:
        vector<vector<Allele*> > trueCombos = trueAlleleCombinations(observedAlleles.size(), *g);
        vector<long double> probs;
        for (vector<vector<Allele*> >::iterator combo = trueCombos.begin(); combo != trueCombos.end(); ++combo) {
            probs.push_back(probObservedAllelesGivenGenotypeAndTrueAlleles(observedAllelesAndProbs, *combo, *g));
        }
        results.push_back(make_pair(*g, logsumexp(probs)));
    }

    return results;

}

//returns prob in log space
long double
Caller::probObservedAllelesGivenGenotypeAndTrueAlleles(
        vector<tuple<long double, long double, Allele* > >& observedAllelesAndProbs,
        vector<Allele*>& trueAlleles,
        vector<Allele>& genotype) {

    int ploidy = genotype.size();
    long double ploidyFactorial = factorialln(ploidy);

    // summation variables for counting our sampling distribution
    long double alleleGenotypeFrequencyProd = 0;
    long double genotypeAlleleFrequencyFactorialProd = 0;
    long double alleleCountFactorialSum;

    vector<vector<Allele*> > trueAlleleGroups = groupAlleles(trueAlleles, allelesEquivalent);
    vector<vector<Allele*> > genotypeAlleleGroups = groupAlleles(genotype, allelesEquivalent);

    for (vector<vector<Allele*> >::const_iterator alleles = genotypeAlleleGroups.begin(); alleles != genotypeAlleleGroups.end(); ++alleles) {
        genotypeAlleleFrequencyFactorialProd += factorialln(alleles->size());
    }

    for (vector<vector<Allele*> >::const_iterator alleles = trueAlleleGroups.begin(); alleles != trueAlleleGroups.end(); ++alleles) {
        alleleGenotypeFrequencyProd += log(probChooseAlleleFromAlleles(*(alleles->front()), genotype));
        alleleCountFactorialSum += factorial(alleles->size());
    }

    //long double samplingProb = (ploidyFactorial + ploidyFactorial) - (genotypeFreqProd + obsFreqSum) + alleleProbProd;
    long double samplingProb = 
        (ploidyFactorial + ploidyFactorial) // squared in log space
            - (genotypeAlleleFrequencyFactorialProd + log(alleleCountFactorialSum)) 
        + alleleGenotypeFrequencyProd;

    vector<long double> obsProbs;
    vector<Allele*>::const_iterator t = trueAlleles.begin();
    for (vector<tuple<long double, long double, Allele* > >::const_iterator p = observedAllelesAndProbs.begin();
         p != observedAllelesAndProbs.end() && t != trueAlleles.end();
         ++p, ++t) {
        if (p->get<2>()->equivalent(**t)) {
            obsProbs.push_back(p->get<0>());
        } else {
            obsProbs.push_back(p->get<1>());
        }
    }
    
    return logsumexp(obsProbs) + samplingProb;

}

// the true allele combinations / permutations that can be drawn out of a given genotype
// this will be a very, very large set for even modest n and small g
vector<vector<Allele*> > trueAlleleCombinations(int n, Genotype& g) {

    vector<vector<Allele*> > results;

    vector<vector<Allele*> > combinations = multichoose_ptr(n, g);
    for(vector<vector<Allele*> >::iterator av = combinations.begin(); av != combinations.end(); ++av) {
        vector<vector<Allele*> > perms = multipermute(*av);
        for(vector<vector<Allele*> >::const_iterator p = perms.begin(); p != perms.end(); ++p) {
            results.push_back(*p);
        }
    }

    return results;

}

// log in and out probabilities per allele, caching step
void Caller::calculateAlleleBinaryProbabilities(
        vector<Allele*>& alleles, 
        vector<tuple<long double, long double, Allele* > >& allelesAndProbs) {

    for (vector<Allele*>::const_iterator a = alleles.begin(); a != alleles.end(); ++a) {
        Allele *allele = *a;
        long double q = allele->Quality(currentPosition);
        allelesAndProbs.push_back(make_tuple(log(1 - phred2float(q)), phred2ln(q), allele));
    }

}


// TODO move to log space
void normalizeGenotypeProbabilities(vector<pair<Genotype, long double> >& genotypeProbabilities) {
    long double sum = 0;
    for (vector<pair<Genotype, long double> >::const_iterator p = genotypeProbabilities.begin(); p != genotypeProbabilities.end(); ++p) {
        sum += p->second;
    }
    if (sum != 0) { // -nan guard
        for (vector<pair<Genotype, long double> >::iterator p = genotypeProbabilities.begin(); p != genotypeProbabilities.end(); ++p) {
            p->second /= sum;
        }
    }
}

void normalizeGenotypeProbabilitiesln(vector<pair<Genotype, long double> >& genotypeProbabilities) {
    long double sum = 0;
    for (vector<pair<Genotype, long double> >::const_iterator p = genotypeProbabilities.begin(); p != genotypeProbabilities.end(); ++p) {
        sum += exp(p->second);
    }
    sum = log(sum);
    //if (sum != 0) { // -nan guard
    for (vector<pair<Genotype, long double> >::iterator p = genotypeProbabilities.begin(); p != genotypeProbabilities.end(); ++p) {
        p->second -= sum;
    }
    //}
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

long double Caller::probObservedAllelesGivenGenotype(vector<vector<Allele*> > &alleleCombo, vector<Allele> &genotype) {

    int ploidy = genotype.size();

    LOG2(genotype);

    // counts indexed by genotype position
    //vector<int> sumQ (ploidy, 0);
    //vector<int> counts (ploidy, 0);
    int outCount = 0, inCount = 0;
    long double inProb = 1;
    long double outProb = 1;
    long double inProbln = 0;
    long double outProbln = 0;
    int obsCount = 0;

    // (1)
    // p( observed allele | true allele ) = 1 - 10 ^(-Q/10)  if observed allele in true allele
    //                          ... and   = 10 ^(-Q/10) if observed allele != true allele
    //             where Q is the quality value associated with the observed allele

    long double factAlleleCountProduct = 1;
    long double probObsAllelesProduct = 1;
    for (vector<vector<Allele*> >::iterator alleleObservations = alleleCombo.begin();
            alleleObservations != alleleCombo.end(); ++alleleObservations) {
        LOG2(*alleleObservations);
        
        // NB: Not the greatest solution.  This should be improved.
        // Avoids processing the same observations more than once.
        // We can safely do this because multichoose provides the allele
        // combination multisets in lexographic order.
        if ((alleleObservations + 1) != alleleCombo.end() 
                && (alleleObservations + 1)->front()->equivalent(*alleleObservations->front())) {
            continue;
        }
    // (2)
    // P = ( ploidy! / product( alleleCount! for alleleType in alleleCombo ) ) * product( p_^alleleCount for alleleType in alleleCombo )
        long double multinomialProb = pow(probChooseAlleleFromAlleles(*alleleObservations->front(), genotype), 
                alleleObservations->size());
        if (multinomialProb > 0) {
            probObsAllelesProduct *= multinomialProb;
            factAlleleCountProduct *= factorial(alleleObservations->size());
            obsCount += alleleObservations->size();
        }

        bool in = false;
        int i = 0;
        for (vector<Allele>::iterator g = genotype.begin(); g != genotype.end(); g++) {
            if (g->equivalent(*alleleObservations->front())) {
                in = true; break;
            }
        }

        if (in) {
            LOG2(*alleleObservations << " in " << genotype);
            for (vector<Allele*>::iterator obs = alleleObservations->begin(); 
                    obs != alleleObservations->end(); obs++) {
                inProb *= 1 - phred2float((*obs)->Quality(currentPosition));
                //inProbln += (*obs)->Quality(currentPosition);
                ++inCount;
            }
        } else {
            LOG2(*alleleObservations << " not in " << genotype);
            for (vector<Allele*>::iterator obs = alleleObservations->begin(); 
                    obs != alleleObservations->end(); obs++) {
                //outProb *= phred2float((*obs)->Quality(currentPosition));
                outProbln += (*obs)->Quality(currentPosition);
                ++outCount;
            }
        }
    }

    // (1)
    long double probAlleleObsGivenGenotype = ( inCount ? inProb : 1 ) * ( outCount ? phred2float(outProbln) : 1 );
    LOG2(genotype);
    LOG2("outCount:" << outCount << ";"
        << "outProb:" << phred2float(outProbln) << ";"
        << "outProbln:" << outProbln << ";"
        << "inCount:" << inCount << ";"
        << "inProb:" << inProb << ";"
        << "inProbln:" << inProbln << ";"
        << "probAlleleObsGivenGenotype:" << probAlleleObsGivenGenotype);
    
    // (2)
    long double alleleSamplingProbability = ( (long double) factorial(obsCount) / factAlleleCountProduct ) * probObsAllelesProduct;
    LOG2("factAlleleCountProduct:" << factAlleleCountProduct << ";"
        << "obsCount:" << obsCount << ";"
        << "probObsAllelesProduct:" << probObsAllelesProduct << ";"
        << "alleleSamplingProbability:" << alleleSamplingProbability);

    // (1) * (2)
    return probAlleleObsGivenGenotype * alleleSamplingProbability;
}

// uses cached probabilities to boost performance when performing repetitive calculations
//
// tuple structure:
//     < probability all observed alleles are true,
//       probability all observed alleles are false,
//       vector of observed alleles >
//
long double Caller::probObservedAllelesGivenGenotype(
        vector<tuple<long double, long double, vector<Allele*> > > &alleleGroupProbs, 
        vector<Allele> &genotype) {

    int ploidy = genotype.size();

    LOG2(genotype);

    // counts indexed by genotype position
    //vector<int> sumQ (ploidy, 0);
    //vector<int> counts (ploidy, 0);
    int outCount = 0, inCount = 0;
    long double inProbln = 0;
    long double outProbln = 0;
    int obsCount = 0;

    // (1)
    // p( observed allele | true allele ) = 1 - 10 ^(-Q/10)  if observed allele in true allele
    //                          ... and   = 10 ^(-Q/10) if observed allele != true allele
    //             where Q is the quality value associated with the observed allele

    long double factAlleleCountProductln = 0;
    long double probObsAllelesProductln = 0;
    for (vector<tuple<long double, long double, vector<Allele*> > >::iterator obs = alleleGroupProbs.begin();
            obs != alleleGroupProbs.end(); ++obs) {

        vector<Allele*>& alleles = obs->get<2>();
        LOG2(alleles);
        
        // NB: Not the greatest solution.  This should be improved.
        // Avoids processing the same observations more than once.
        // We can safely do this because multichoose provides the allele
        // combination multisets in lexographic order.
        if ((obs + 1) != alleleGroupProbs.end() 
                && (obs + 1)->get<2>().front()->equivalent(*alleles.front())) {
            continue;
        }
    // (2)
    // P = ( ploidy! / product( alleleCount! for alleleType in alleleCombo ) ) * product( p_^alleleCount for alleleType in alleleCombo )
        long double multinomialProb = pow(probChooseAlleleFromAlleles(*alleles.front(), genotype), 
                alleles.size());
        if (multinomialProb > 0) {
            probObsAllelesProductln += log(multinomialProb);
            factAlleleCountProductln += log(factorial(alleles.size()));
            obsCount += alleles.size();
        }

        // check if this set of observations supports our genotype
        bool in = false;
        int i = 0;
        for (vector<Allele>::iterator g = genotype.begin(); g != genotype.end(); ++g) {
            if (g->equivalent(*alleles.front())) {
                in = true; break;
            }
        }

        // if it supports, add to our in* prob and count
        if (in) {
            LOG2(alleles << " in " << genotype);
            inProbln += obs->get<0>();
            inCount += alleles.size();
        } else {
            LOG2(alleles << " not in " << genotype);
            outProbln += obs->get<1>();
            outCount += alleles.size();
        }
    }

    // (1)
    long double probAlleleObsGivenGenotypeln = ( inCount ? inProbln : 0 ) + ( outCount ? outProbln : 0 );
    LOG2(genotype);
    LOG2("outCount:" << outCount << ";"
        << "outProbln:" << outProbln << ";"
        << "inCount:" << inCount << ";"
        << "inProbln:" << inProbln << ";"
        << "probAlleleObsGivenGenotypeln:" << probAlleleObsGivenGenotypeln);
    
    // (2)
    long double alleleSamplingProbability = ( (long double) log(factorial(obsCount)) - factAlleleCountProductln ) + probObsAllelesProductln;
    LOG2("factAlleleCountProduct:" << factAlleleCountProductln << ";"
        << "obsCount:" << obsCount << ";"
        << "probObsAllelesProduct:" << probObsAllelesProductln << ";"
        << "alleleSamplingProbability:" << alleleSamplingProbability);

    // (1) * (2)
    return probAlleleObsGivenGenotypeln + alleleSamplingProbability;
}


int observationsInAlleleCombo(vector<vector<Allele> > &combo) {
    int count = 0;
    for (vector<vector<Allele> >::iterator obs = combo.begin(); obs != combo.end(); obs++) {
        if ((obs + 1) != combo.end()
                && (obs + 1)->front().equivalent(obs->front())) {
            continue;
        }
        count += obs->size();
    }
    return count;
}

long double probChooseAlleleFromAlleles(Allele &allele, vector<Allele> &alleles) {

    int matches = 0;

    for (vector<Allele>::iterator ai = alleles.begin(); ai != alleles.end(); ai++) {
        if (ai->equivalent(allele))
            ++matches;
    }

    return (long double) matches / (long double) alleles.size();

}

// gets an approximate normalizer for our bayesian posterior
long double approximateBayesianNormalizationFactor(vector<vector<Allele> > &genotypes,
        vector<vector<long double> > &probGenotypesGivenSampleObs,
        vector<vector<Allele> > &sampleGroups) {

    long double permutationsCount =  factorial(sampleGroups.size()) / factorial(sampleGroups.size() - genotypes.size());
    long double probSum = 0;
    int probsCount = 0;

    for (vector<vector<long double> >::const_iterator s =
            probGenotypesGivenSampleObs.begin(); 
            s != probGenotypesGivenSampleObs.end(); s++) {
        probsCount += s->size();
        for (vector<long double>::const_iterator r = s->begin(); r != s->end(); r++)
            probSum += *r;
    }

    long double averageProb = probSum / probsCount;

    // *gross* approximation

    return averageProb * permutationsCount;

}

// very, very, (impossibly,) computationally intensive
long double bayesianNormalizationFactor(vector<vector<Allele> > &genotypes,
        vector<vector<long double> > &probGenotypesGivenSampleObs,
        vector<vector<Allele> > &sampleGroups) {

    long double probSum = 0;
    int probsCount = 0;
    int permutationsCount = 0;
    
    vector<int> indexes;
    for (int i=0; i<genotypes.size(); ++i)
        indexes.push_back(i);

    vector<vector<int> > genotype_indexes = multichoose(sampleGroups.size(), indexes);

    //int i =0;
    for (vector<vector<int> >::iterator it = genotype_indexes.begin(); it != genotype_indexes.end(); ++it) {
        //cout << i++ << " of " << genotype_indexes.size()  << endl;
        do {
            long double currProb = 1;
            // get the probability of this genotype vector
            int j = 0;
            for (vector<int>::iterator index = it->begin(); index != it->end(); ++index) {
                currProb *= probGenotypesGivenSampleObs.at(j++).at(*index);
            }
            // add it to our probs sum
            probSum += currProb;
            permutationsCount++;
        } while (next_permutation(it->begin(), it->end()));
    }

    return probSum;

}


#include "LeftAlign.h"

//bool debug;

// Attempts to left-realign all the indels represented by the alignment cigar.
//
// This is done by shifting all indels as far left as they can go without
// mismatch, then merging neighboring indels of the same class.  leftAlign
// updates the alignment cigar with changes, and returns true if realignment
// changed the alignment cigar.
//
// To left-align, we move multi-base indels left by their own length as long as
// the preceding bases match the inserted or deleted sequence.  After this
// step, we handle multi-base homopolymer indels by shifting them one base to
// the left until they mismatch the reference.
//
// To merge neighboring indels, we iterate through the set of left-stabilized
// indels.  For each indel we add a new cigar element to the new cigar.  If a
// deletion follows a deletion, or an insertion occurs at the same place as
// another insertion, we merge the events by extending the previous cigar
// element.
//
// In practice, we must call this function until the alignment is stabilized.
//
bool leftAlign(BamAlignment& alignment, string& referenceSequence, bool debug) {

    int arsOffset = 0; // pointer to insertion point in aligned reference sequence
    string alignedReferenceSequence = referenceSequence;
    int aabOffset = 0;
    string alignmentAlignedBases = alignment.QueryBases;

    // store information about the indels
    vector<IndelAllele> indels;

    int rp = 0;  // read position, 0-based relative to read
    int sp = 0;  // sequence position

    string softBegin;
    string softEnd;

    stringstream cigar_before, cigar_after;
    for (vector<CigarOp>::const_iterator c = alignment.CigarData.begin();
        c != alignment.CigarData.end(); ++c) {
        unsigned int l = c->Length;
        char t = c->Type;
        cigar_before << l << t;
        if (t == 'M') { // match or mismatch
            sp += l;
            rp += l;
        } else if (t == 'D') { // deletion
            indels.push_back(IndelAllele(false, l, sp, rp, referenceSequence.substr(sp, l)));
            alignmentAlignedBases.insert(rp + aabOffset, string(l, '-'));
            aabOffset += l;
            sp += l;  // update sample position
        } else if (t == 'I') { // insertion
            indels.push_back(IndelAllele(true, l, sp, rp, alignment.QueryBases.substr(rp, l)));
            alignedReferenceSequence.insert(sp + softBegin.size() + arsOffset, string(l, '-'));
            arsOffset += l;
            rp += l;
        } else if (t == 'S') { // soft clip, clipped sequence present in the read not matching the reference
            // remove these bases from the refseq and read seq, but don't modify the alignment sequence
            if (rp == 0) {
                alignedReferenceSequence = string(l, '*') + alignedReferenceSequence;
                softBegin = alignmentAlignedBases.substr(0, l);
            } else {
                alignedReferenceSequence = alignedReferenceSequence + string(l, '*');
                softEnd = alignmentAlignedBases.substr(alignmentAlignedBases.size() - l, l);
            }
            rp += l;
        } else if (t == 'H') { // hard clip on the read, clipped sequence is not present in the read
        } else if (t == 'N') { // skipped region in the reference not present in read, aka splice
            sp += l;
        }
    }


    int alignedLength = sp;

    DEBUG("| " << cigar_before.str() << endl
       << "| " << alignedReferenceSequence << endl
       << "| " << alignmentAlignedBases << endl);

    // if no indels, return the alignment
    if (indels.empty()) { return false; }

    // for each indel, from left to right
    //     while the indel sequence repeated to the left and we're not matched up with the left-previous indel
    //         move the indel left

    vector<IndelAllele>::iterator previous = indels.begin();
    for (vector<IndelAllele>::iterator id = indels.begin(); id != indels.end(); ++id) {

        // left shift
        //
        // from 1 base to the length of the indel, attempt to shift left
        // if the move would cause no change in alignment optimality (no
        // introduction of mismatches, and by definition no change in gap
        // length), move to the new position.
        // in practice this moves the indel left when we reach the size of
        // the repeat unit.
        //
        int steppos, readsteppos;
        IndelAllele& indel = *id;
        int i = 1;
        if (!indel.homopolymer() && indel.length % 2 != 0) {
            i = indel.length;
        }
        while (i <= indel.length) {

            int steppos = indel.position - i;
            int readsteppos = (indel.insertion ? indel.readPosition - i : indel.readPosition) - 1;

#ifdef VERBOSE_DEBUG
            if (debug) {
                if (steppos >= 0 && readsteppos >= 0) {
                    cerr << referenceSequence.substr(steppos, indel.length) << endl;
                    if (indel.insertion) {
                        cerr << alignment.QueryBases.substr(readsteppos, indel.length) << endl;
                    } else {
                        cerr << indel.sequence << endl;
                    }
                }
            }
#endif
            while (steppos >= 0 && readsteppos >= 0
                   && indel.sequence == referenceSequence.substr(steppos, indel.length)
                   && (!indel.insertion || indel.sequence == alignment.QueryBases.substr(readsteppos, indel.length))
                   && (id == indels.begin()
                       || (previous->insertion && steppos >= previous->position)
                       || (!previous->insertion && steppos >= previous->position + previous->length))) {
                DEBUG("indel " << indel << " shifting " << i << "bp left" << endl);
                indel.position -= i;
                indel.readPosition -= i;
                steppos = indel.position - i;
                readsteppos = (indel.insertion ? indel.readPosition - i : indel.readPosition) - 1;
            }
            ++i;
        }

        // left shift deletions with exchangeable flanking sequence
        //
        // for example:
        //
        // GTTACGTT           GTTACGTT
        // GT-----T  becomes  G-----TT
        //
        steppos = indel.position - 1;
        readsteppos = (indel.insertion ? indel.readPosition - 1 : indel.readPosition) - 1;
        while (steppos >= 0 && readsteppos >= 0
               && indel.sequence.at(indel.sequence.size() - 1) == referenceSequence.at(indel.position - 1)
               && indel.sequence == referenceSequence.substr(steppos, indel.length)
               && (!indel.insertion || indel.sequence == alignment.QueryBases.substr(readsteppos, indel.length))
               && (id == indels.begin()
                   || (previous->insertion && indel.position - 1 >= previous->position)
                   || (!previous->insertion && indel.position - 1 >= previous->position + previous->length))) {
            DEBUG("indel " << indel << " exchanging bases " << 1 << "bp left" << endl);
            indel.sequence = indel.sequence.at(indel.sequence.size() - 1) + indel.sequence.substr(0, indel.sequence.size() - 1);
            indel.position -= 1;
            indel.readPosition -= 1;
            steppos = indel.position - 1;
            readsteppos = (indel.insertion ? indel.readPosition - 1 : indel.readPosition) - 1;
        }
        // tracks previous indel, so we don't run into it with the next shift
        previous = id;
    }

    // bring together floating indels
    // from left to right
    // check if we could merge with the next indel
    // if so, adjust so that we will merge in the next step
    if (indels.size() > 1) {
        previous = indels.begin();
        for (vector<IndelAllele>::iterator id = (indels.begin() + 1); id != indels.end(); ++id) {
            IndelAllele& indel = *id;
            // could we shift right and merge with the previous indel?
            // if so, do it
            int prev_end = previous->insertion ? previous->position : previous->position + previous->length;
            if (previous->insertion == indel.insertion
                    && ((previous->insertion && previous->position < indel.position)
                        ||
                        (!previous->insertion && previous->position + previous->length < indel.position))) {
                if (previous->homopolymer()) {
                    string seq = referenceSequence.substr(prev_end, indel.position - prev_end);
                    if (previous->sequence.at(0) == seq.at(0) && homopolymer(seq)) {
                        DEBUG("moving " << *previous << " right to " 
                                << (indel.insertion ? indel.position : indel.position - previous->length) << endl);
                        previous->position = indel.insertion ? indel.position : indel.position - previous->length;
                    }
                } 
                else {
                    int pos = previous->position;
                    while (pos < (int) referenceSequence.length() &&
                            ((previous->insertion && pos + previous->length <= indel.position)
                            ||
                            (!previous->insertion && pos + previous->length < indel.position))
                            && previous->sequence 
                                == referenceSequence.substr(pos + previous->length, previous->length)) {
                        pos += previous->length;
                    }
                    if (pos < previous->position &&
                        ((previous->insertion && pos + previous->length == indel.position)
                        ||
                        (!previous->insertion && pos == indel.position - previous->length))
                       ) {
                        DEBUG("right-merging tandem repeat: moving " << *previous << " right to " << pos << endl);
                        previous->position = pos;
                    }
                }
            }
        }
    }

    // for each indel
    //     if ( we're matched up to the previous insertion (or deletion) 
    //          and it's also an insertion or deletion )
    //         merge the indels

    vector<CigarOp> newCigar;

    if (!softBegin.empty()) {
        newCigar.push_back(CigarOp('S', softBegin.size()));
    }

    vector<IndelAllele>::iterator id = indels.begin();
    IndelAllele last = *id++;
    if (last.position > 0) {
        newCigar.push_back(CigarOp('M', last.position));
        newCigar.push_back(CigarOp((last.insertion ? 'I' : 'D'), last.length));
    } else {
        newCigar.push_back(CigarOp((last.insertion ? 'I' : 'D'), last.length));
    }
    int lastend = last.insertion ? last.position : (last.position + last.length);
    DEBUG(last << ",");

    for (; id != indels.end(); ++id) {
        IndelAllele& indel = *id;
        DEBUG(indel << ",");
        if (indel.position < lastend) {
            cerr << "impossibility?: indel realigned left of another indel" << endl << alignment.Name
                << " " << alignment.Position << endl << alignment.QueryBases << endl;
            exit(1);
        } else if (indel.position == lastend && indel.insertion == last.insertion) {
            CigarOp& op = newCigar.back();
            op.Length += indel.length;
        } else if (indel.position > lastend) {
            newCigar.push_back(CigarOp('M', indel.position - lastend));
            newCigar.push_back(CigarOp((indel.insertion ? 'I' : 'D'), indel.length));
        }
        last = *id;
        lastend = last.insertion ? last.position : (last.position + last.length);
    }
    
    if ((last.position + last.length) < alignedLength) {
        newCigar.push_back(CigarOp('M', alignedLength - lastend));
    }

    if (!softEnd.empty()) {
        newCigar.push_back(CigarOp('S', softEnd.size()));
    }

    DEBUG(endl);

#ifdef VERBOSE_DEBUG
    if (debug) {
        for (vector<CigarOp>::const_iterator c = alignment.CigarData.begin();
            c != alignment.CigarData.end(); ++c) {
            unsigned int l = c->Length;
            char t = c->Type;
            cerr << l << t;
        }
        cerr << endl;
    }
#endif

    alignment.CigarData = newCigar;

    for (vector<CigarOp>::const_iterator c = alignment.CigarData.begin();
        c != alignment.CigarData.end(); ++c) {
        unsigned int l = c->Length;
        char t = c->Type;
        cigar_after << l << t;
    }
    DEBUG(cigar_after.str() << endl);

    // check if we're realigned
    if (cigar_after.str() == cigar_before.str()) {
        return false;
    } else {
        return true;
    }

}

// Iteratively left-aligns the indels in the alignment until we have a stable
// realignment.  Returns true on realignment success or non-realignment.
// Returns false if we exceed the maximum number of realignment iterations.
//
bool stablyLeftAlign(BamAlignment& alignment, string referenceSequence, int maxiterations, bool debug) {

    if (!leftAlign(alignment, referenceSequence, debug)) {

        DEBUG("did not realign" << endl);
        return true;

    } else {

        while (leftAlign(alignment, referenceSequence, debug) && --maxiterations > 0) {
            DEBUG("realigning ..." << endl);
        }

        if (maxiterations <= 0) {
            return false;
        } else {
            return true;
        }

    }
}

#ifndef __LEFTALIGN_H
#define __LEFTALIGN_H

#include <iostream>
#include <getopt.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <signal.h>
#include <stdlib.h>
#include <cmath>
#include <algorithm>
#include <map>
#include <vector>

#include "Fasta.h"

#ifdef HAVE_BAMTOOLS
#include "api/BamAlignment.h"
#include "api/BamReader.h"
#include "api/BamMultiReader.h"
#include "api/BamWriter.h"
using namespace BamTools;
#define GETNEXT(reader, alignment) reader.GetNextAlignment(alignment)
#define GETERRORSTRING(reader) reader.GetErrorString()
#define ALIGNMENTLENGTH Length
#define ISMAPPED IsMapped()
#define REFID RefID
#define POSITION Position
#define REFNAME RefName
#define REFLEN RefLength
#define REFVEC RefVector
#define REFDATA RefData
#define BAMALIGN BamAlignment
#define QUALITIES Qualities
#define QUERYBASES QueryBases
#define ALIGNEDBASES AlignedBases.size()
#define QNAME Name
#define GETREFDATA GetReferenceData()
#define GETREFNUM GetReferenceCount()
#define GETREFID(name) GetReferenceID(name)
#define ENDPOSITION GetEndPosition()
#define SEQLEN QueryBases.size()
#define MAPPINGQUALITY MapQuality
#define CIGAR std::vector<CigarOp>
#define GETCIGAR CigarData
#define ISDUPLICATE IsDuplicate()
#define ISREVERSESTRAND IsReverseStrand()
#define ISPAIRED IsPaired()
#define ISMATEMAPPED IsMateMapped()
#define ISPROPERPAIR IsProperPair()
#define CIGLEN Length
#define CIGTYPE Type
#define BAMREADER BamMultiReader
#define BAMSINGLEREADER BamReader
#define FILLREADGROUP(rg, align) (align).GetTag("RG", (rg))
#define ADDCIGAR push_back
#define CIGOP CigarOp
#define GETHEADERTEXT GetHeaderText()
#define STDIN "stdin"
#define WRITEALIGNMENT(writer, alignment) writer.SaveAlignment(alignment)
#else

#define GETNEXT(reader, alignment) reader.GetNextRecord(alignment)
#define MAPPINGQUALITY MapQuality()
#define ALIGNMENTLENGTH Length()
#define ISMAPPED MappedFlag()
#define ISPAIRED PairedFlag()
#define ISMATEMAPPED MateMappedFlag()
#define ISPROPERPAIR ProperPair()
#define ISREVERSESTRAND ReverseFlag()
#define SEQLEN Length()
#define BAMALIGN SeqLib::BamRecord
#define REFID ChrID()
#define POSITION Position()
#define REFVEC std::vector<SeqLib::HeaderSequence>
#define REFDATA SeqLib::HeaderSequence
#define REFNAME Name
#define REFLEN Length
#define QUALITIES Qualities()
#define QUERYBASES Sequence()
#define ALIGNEDBASES NumAlignedBases()
#define QNAME Qname()
#define GETREFDATA Header().GetHeaderSequenceVector()
#define GETREFNUM Header().NumSequences()
#define ENDPOSITION PositionEnd()
#define CIGAR SeqLib::Cigar
#define BAMREADER SeqLib::BamReader
#define BAMSINGLEREADER SeqLib::BamReader
#define GETCIGAR GetCigar()
#define GETREFID(name) Header().Name2ID(name)
#define ISDUPLICATE DuplicateFlag()
#define CIGLEN Length()
#define CIGTYPE Type()
#define ADDCIGAR add
#define CIGOP SeqLib::CigarField
#define FILLREADGROUP(rg, align) (rg) = (align).GetZTag("RG")
#define GETHEADERTEXT HeaderConcat()
#include "SeqLib/BamReader.h"
#include "SeqLib/BamWriter.h"
#define STDIN "-"
#define WRITEALIGNMENT(writer, alignment) writer.WriteRecord(alignment)
#endif




#include "IndelAllele.h"

#ifdef VERBOSE_DEBUG
#define LEFTALIGN_DEBUG(msg) \
    if (debug) { cerr << msg; }
#else
#define LEFTALIGN_DEBUG(msg)
#endif

using namespace std;

bool leftAlign(BAMALIGN& alignment, string& referenceSequence, bool debug = false);
bool stablyLeftAlign(BAMALIGN& alignment, string referenceSequence, int maxiterations = 20, bool debug = false);
int countMismatches(BAMALIGN& alignment, string referenceSequence);

#endif

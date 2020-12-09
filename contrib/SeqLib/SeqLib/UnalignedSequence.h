#ifndef SEQLIB_UNALIGNED_SEQ_H
#define SEQLIB_UNALIGNED_SEQ_H

extern "C" {
  #include "bwa/bwa.h"
  #include "bwa/bwt.h"
  #include "bwa/bntseq.h"
  #include "bwa/kseq.h"
  #include <stdlib.h>
  #include "bwa/utils.h"
  #include "bwa/bwamem.h"
  int is_bwt(ubyte_t *T, int n);
  KSEQ_DECLARE(gzFile)
}


#include <cstring>
#include <vector>
#include <iostream>

namespace SeqLib {

  /** Structure to hold unaligned sequence (name and bases)
   */
  struct UnalignedSequence {
  
    /** Construct an empty sequence */
    UnalignedSequence() {}
  
    /** Construct an unaligned sequence with name and sequence
     * @param n Name of the sequence 
     * @param s Sequence, stored as ACTG or N characters
     */
    UnalignedSequence(const std::string& n, const std::string& s) : Name(n), Com(std::string()), Seq(s), Qual(std::string()), Strand('*') {}

    /** Construct an unaligned sequence with name, sequence and quality score
     * @param n Name of the sequence 
     * @param s Sequence, stored as ACTG or N characters
     * @param q Quality string
     */
    UnalignedSequence(const std::string& n, const std::string& s, const std::string& q) : Name(n), Com(std::string()),  Seq(s), Qual(q), Strand('*') {}

    /** Construct an unaligned sequence with name, sequence, quality score and strand
     * @param n Name of the sequence 
     * @param s Sequence, stored as ACTG or N characters
     * @param q Quality string
     * @param t Strand of the sequence, one of '*', '+', '-'
     */
    UnalignedSequence(const std::string& n, const std::string& s, const std::string& q, char t) : Name(n), Com(std::string()), Seq(s), Qual(q), Strand(t) {}

     std::string Name; ///< Name of the contig
     std::string Com;  ///< Comment of the contig
     std::string Seq;  ///< Sequence of the contig (upper-case ACTGN)
     std::string Qual; ///< Quality scores
     char Strand;      ///< Strand of the sequence. Default is '*'
     
    /** Output an unaligned sequence to ostream
     * @param os ostream
     * @param us UnalignedSequence
     */
    friend std::ostream& operator<<(std::ostream& os, const SeqLib::UnalignedSequence& us){
        os << "@" << us.Name << " " << us.Com << "\n";
        os << us.Seq << "\n+\n";
        os << us.Qual << "\n";
        return os;
    }
  };

  typedef std::vector<UnalignedSequence> UnalignedSequenceVector; ///< A collection of unaligned sequences

}

#endif

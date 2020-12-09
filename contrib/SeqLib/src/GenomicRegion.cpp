#include "SeqLib/GenomicRegion.h"

#include <cassert>
#include <stdexcept>
#include <climits>
#ifdef HAVE_C11
#include <regex>
#endif

// 4 billion
#define END_MAX 4000000000

namespace SeqLib {

// return the width of the genomic region
int GenomicRegion::Width() const {
  return pos2 - pos1 + 1;
}

// returns 0 for no overlaps, 1 for partial and 2 for complete
int GenomicRegion::GetOverlap(const GenomicRegion& gr) const {

  if (gr.chr != chr)
    return 0;
  
  // argument pos1 is in
  bool gr1_in = gr.pos1 >= pos1 && gr.pos1 <= pos2;
  // argument pos2 is in
  bool gr2_in = gr.pos2 >= pos1 && gr.pos2 <= pos2;
  // object pos1 is in
  bool pos1_in = pos1 >= gr.pos1 && pos1 <= gr.pos2;
  // object pos2 is in
  bool pos2_in = pos2 >= gr.pos1 && pos2 <= gr.pos2;

  // object is in the argument
  if (pos1_in && pos2_in) 
    return 3;

  // argument is in the oboject
  if ( gr1_in && gr2_in)
    return 2;

  // partial overlap
  if (gr1_in || gr2_in || pos1_in || pos2_in)
    return 1;

  return 0;

}


  std::string GenomicRegion::ChrName(const BamHeader& h) const {
    
    std::string cc;
    if (!h.isEmpty()) {
      if (chr >= h.NumSequences())
	throw std::invalid_argument( "GenomicRegion::ChrName - not enough targets in BamHeader to cover ref id");
      else
	cc = h.IDtoName(chr); // std::string(h->target_name[chr]);
    } else {
      cc = chrToString(chr);
    }
    return cc;
  }

  
  std::string GenomicRegion::PointString(const BamHeader& h) const {
    std::stringstream out;
    out << ChrName(h) << ":" << SeqLib::AddCommas<int>(pos1) << "(" << strand << ")";
    return out.str();
  }

void GenomicRegion::Pad(int32_t pad) {

  if (-pad*2 > Width())
    throw std::out_of_range(
         "GenomicRegion::pad - negative pad values can't obliterate GenomicRegion with val " + 
	 tostring(chr) + ":" + tostring(pos1) + "-" + tostring(pos2) + 
	 " and pad " + tostring(pad));

  pos1 -= pad;
  pos2 += pad;

  //if (pad > pos1)
  //  pos1 = 1;
  //else
  //  pos1 = pos1-pad;

  //const int32_t maxpos = 250000000;
  //pos2 = std::min(pos2+pad, maxpos); // 2500000000 is dummy for now. should be chr end

}

bool GenomicRegion::operator<(const GenomicRegion& b) const {
  return (chr < b.chr) || (chr == b.chr && pos1 < b.pos1) || (chr==b.chr && pos1 == b.pos1 && pos2 < b.pos2);
}

bool GenomicRegion::operator>(const GenomicRegion& b) const {
  return !(*this == b) && !(*this < b);
}

bool GenomicRegion::operator==(const GenomicRegion &b) const {
  return (chr == b.chr && pos1 == b.pos1 && b.pos2 == pos2);
}

bool GenomicRegion::operator!=(const GenomicRegion &b) const {
  return !(*this == b); 
}

bool GenomicRegion::operator<=(const GenomicRegion &b) const {
  return (*this < b || *this == b);
}

bool GenomicRegion::operator>=(const GenomicRegion &b) const {
  return (*this > b || *this == b);
}

  /*  std::string GenomicRegion::ToString() const {
    return chrToString(chr) + ":" + SeqLib::AddCommas<int>(pos1) + "-" + AddCommas<int>(pos2) + "(" +
      strand + ")"; 
      }*/

  std::string GenomicRegion::ToString(const BamHeader& h) const {
    return ChrName(h) + ":" + SeqLib::AddCommas<int>(pos1) + "-" + AddCommas<int>(pos2) + "(" +
      strand + ")"; 
  }


std::ostream& operator<<(std::ostream& out, const GenomicRegion& gr) {
  out << gr.chrToString(gr.chr) << ":" << SeqLib::AddCommas<int>(gr.pos1) << "-" << AddCommas<int>(gr.pos2) << "(" << 
    gr.strand << ")"; 
  return out;
}

  GenomicRegion::GenomicRegion(const std::string& reg, const BamHeader& hdr) {
  
  if (hdr.isEmpty())
    throw std::invalid_argument("GenomicRegion constructor - supplied empty BamHeader");

  // scrub String
  //std::string reg2 = SeqLib::scrubString(reg, "chr");

  // use htslib region parsing code
  int tid, beg, end;
  const char * q = hts_parse_reg(reg.c_str(), &beg, &end);
  if (q) {
    char *tmp = (char*)alloca(q - reg.c_str() + 1); // stack alloc
    strncpy(tmp, reg.c_str(), q - reg.c_str());
    tmp[q - reg.c_str()] = 0;
    tid = hdr.Name2ID(std::string(tmp)); //bam_name2id(h.get(), tmp);
    if (tid < 0) {
      std::string inv = "GenomicRegion constructor: Failed to set region for " + reg;
      throw std::invalid_argument(inv);
    }

    if (end == INT_MAX) { // single chrome
      tid = hdr.Name2ID(reg);
      beg = 0;
      end = hdr.GetSequenceLength(reg);
    }
  } else {
    std::string inv = "GenomicRegion constructor: Failed to set region for " + reg;
    throw std::invalid_argument(inv);
  }
  
  chr = tid;
  pos1 = beg+1;
  pos2 = end;
  strand = '*';

}

// constructor to take a pair of coordinates to define the genomic interval
GenomicRegion::GenomicRegion(int32_t t_chr, int32_t t_pos1, int32_t t_pos2, char t_strand) {

  if (t_pos2 < t_pos1 )
    throw std::invalid_argument( "GenomicRegion constructor: end pos must be >= start pos" );

  if ( !(t_strand == '+' || t_strand == '-' || t_strand == '*') )
    throw std::invalid_argument( "GenomicRegion constructor: strand must be one of +, -, *" );

  chr = t_chr;
  pos1 = t_pos1;
  pos2 = t_pos2;
  strand = t_strand;

}
  //private
std::string GenomicRegion::chrToString(int32_t ref) const {

  std::string ref_id;
  if (ref < 0)
    ref_id = tostring(ref);

  if (ref == 22)
    ref_id = "X";
  else if (ref == 23)
    ref_id = "Y";
  else if (ref == 24)
    ref_id = "M";
  else if (ref >= 0)
    ref_id = tostring(ref+1);
  assert(ref_id != "23");
  return ref_id;
}

// checks whether a GenomicRegion is empty
bool GenomicRegion::IsEmpty() const {
  return chr == -1 && pos1 == 0 && pos2 == 0;
}


int32_t GenomicRegion::DistanceBetweenStarts(const GenomicRegion &gr) const {

  if (gr.chr != chr)
    return -1;
  else
    return std::abs(pos1 - gr.pos1);//((pos1 > gr.pos1) ? (pos1 - gr.pos1) : (gr.pos1 - pos1));

}

int32_t GenomicRegion::DistanceBetweenEnds(const GenomicRegion &gr) const {

  if (gr.chr != chr)
    return -1;
  else
    return std::abs(pos2 - gr.pos2);

}


  /*void GenomicRegion::Random() {
  
  uint32_t big = rand() % SeqLib::genome_size_XY;
  //SeqLib::genRandomValue(big, SeqLib::genome_size_XY, seed);
  
  for (size_t k = 0; k < 25; ++k)
    if (big < SeqLib::CHR_CLEN[k]) {
      assert(k > 0);
      chr = --k;
      assert(big > SeqLib::CHR_CLEN[chr]);
      pos1 = big - SeqLib::CHR_CLEN[chr];
      pos2 = pos1;
      return;
    }
  
  std::cerr << "Value of " << big << " outside of expected range."  << std::endl;
  assert(false);
  
  }*/

  GenomicRegion::GenomicRegion(const std::string& tchr, const std::string& tpos1, const std::string& tpos2, const SeqLib::BamHeader& hdr)
  {
    strand = '*';
    // convert the pos strings
    // throws invalid_argument if conversion can't be performed
    // or throws an out_of_range if it is too big for result
#ifdef HAVE_C11
    pos1 = std::stoi(tpos1);
    pos2 = std::stoi(tpos2);
#else
    pos1 = std::atoi(tpos1.c_str());
    pos2 = std::atoi(tpos2.c_str());
#endif
    
    // if no header, assume that it is "standard"
    if (hdr.isEmpty()) {
      if (tchr == "X" || tchr == "chrX")
	chr = 22;
      else if (tchr == "Y" || tchr == "chrY")
	chr = 23;
      else 
#ifdef HAVE_C11
	chr = std::stoi(SeqLib::scrubString(tchr, "chr")) - 1;
#else
	chr = std::atoi(SeqLib::scrubString(tchr, "chr").c_str());
#endif
      return;
    } else {
      chr = hdr.Name2ID(tchr); //bam_name2id(hdr.get(), tchr.c_str());

      // if tchr was not found in sequence dictionary, it's possible that we're 
      // specifying our contigs in b37 style whereas dict is hg** style;
      // let's attempt to automatically convert [0-9XY]+ -> chr[0-9XY]+
#ifdef HAVE_C11
      static std::regex b37_regex("[0-9XY]+");
      if(chr == -1 && std::regex_match(tchr, b37_regex)) {
	 chr = hdr.Name2ID("chr" + tchr);
      }
#endif
    }
  }
}


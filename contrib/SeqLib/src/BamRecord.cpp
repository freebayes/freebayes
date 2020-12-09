#include "SeqLib/BamRecord.h"

#include <cassert>
#include <bitset>
#include <cctype>
#include <stdexcept>

#include "SeqLib/ssw_cpp.h"

#define TAG_DELIMITER "^"
#define CTAG_DELIMITER '^'

namespace SeqLib {

  const int CigarCharToInt[128] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, //0-9
                                     -1,-1,-1,-1,-1,-1,-1,-1,-1,-1, //10-19
                                     -1,-1,-1,-1,-1,-1,-1,-1,-1,-1, //20
                                     -1,-1,-1,-1,-1,-1,-1,-1,-1,-1, //30
                                     -1,-1,-1,-1,-1,-1,-1,-1,-1,-1, //40
                                     -1,-1,-1,-1,-1,-1,-1,-1,-1,-1, //50
                                     -1,BAM_CEQUAL,-1,-1,-1,-1,BAM_CBACK,-1,BAM_CDEL,-1, //60-69
                                     -1,-1,BAM_CHARD_CLIP,BAM_CINS,-1,-1,-1,BAM_CMATCH,BAM_CREF_SKIP,-1,
                                     BAM_CPAD,-1,-1,BAM_CSOFT_CLIP,-1,-1,-1,-1,BAM_CDIFF,-1,
                                     -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
                                     -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
                                     -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
                                     -1,-1,-1,-1,-1,-1,-1,-1};



  struct free_delete {
    void operator()(void* x) { bam_destroy1((bam1_t*)x); }
  };
  
  void BamRecord::init() {
    bam1_t* f = bam_init1();
    b = SeqPointer<bam1_t>(f, free_delete());
  }

  void BamRecord::assign(bam1_t* a) { 
    b = SeqPointer<bam1_t>(a, free_delete()); 
  }

  int32_t BamRecord::PositionWithSClips() const {
    if(!b) return -1; // to be consistent with BamRecord::Position()

    uint32_t* cig = bam_get_cigar(b);
    return ((*cig) & 0xF) == BAM_CSOFT_CLIP ? b->core.pos - ((*cig) >> 4) : b->core.pos;
  }

  int32_t BamRecord::PositionEnd() const { 
    return b ? (b->core.l_qseq > 0 ? bam_endpos(b.get()) : b->core.pos + GetCigar().NumQueryConsumed()) : -1;
  }

  int32_t BamRecord::PositionEndWithSClips() const {
    if(!b) return -1; // to be consistent with BamRecord::PositionEnd()

    uint32_t* cig_last = bam_get_cigar(b) + b->core.n_cigar - 1;
    if(b->core.l_qseq > 0) {
      return ((*cig_last) & 0xF) == BAM_CSOFT_CLIP ? bam_endpos(b.get()) + ((*cig_last) >> 4) :
                                                     bam_endpos(b.get());
    } else {
      return b->core.pos + GetCigar().NumQueryConsumed();
    }
  }

  int32_t BamRecord::PositionEndMate() const { 
    return b ? (b->core.mpos + (b->core.l_qseq > 0 ? b->core.l_qseq : GetCigar().NumQueryConsumed())) : -1;
  }

  GenomicRegion BamRecord::AsGenomicRegion() const {
    char s = '*';
    if (MappedFlag())
      s = ReverseFlag() ? '-' : '+';
    return GenomicRegion(b->core.tid, b->core.pos, PositionEnd(), s);
  }

  GenomicRegion BamRecord::AsGenomicRegionMate() const {
    char s = '*';
    if (MateMappedFlag())
      s = MateReverseFlag() ? '-' : '+';
    return GenomicRegion(b->core.mtid, b->core.mpos, PositionEndMate(), s);
  }

  std::string BamRecord::Sequence() const {
    uint8_t * p = bam_get_seq(b);
    std::string out(b->core.l_qseq, 'N');
    for (int32_t i = 0; i < b->core.l_qseq; ++i) 
      out[i] = BASES[bam_seqi(p,i)];
    return out;
    
  }

  void BamRecord::SetCigar(const Cigar& c) {

    // case where they are equal, just swap them out
    if (c.size() == b->core.n_cigar) {
      b->core.n_cigar = c.size();
      uint32_t * cigr = bam_get_cigar(b);
      for (size_t i = 0; i < b->core.n_cigar; ++i)
	cigr[i] = c[i].raw();
      return;
    }

    // make the new cigar structure
    uint32_t* new_cig = (uint32_t*)malloc(4 * c.size());
    for (size_t i = 0; i < c.size(); ++i)
      new_cig[i] = c[i].raw();
    
    int new_size = b->l_data - (b->core.n_cigar<<2) + (c.size()<<2);
    int old_seqaux_spot = (b->core.n_cigar<<2) + b->core.l_qname;
    int old_seqaux_len = bam_get_l_aux(b) + ((b->core.l_qseq + 1)>>1) + b->core.l_qseq;
    
    // set the new cigar size
    b->core.n_cigar = c.size();

    // copy out all the old data
    uint8_t* oldd = (uint8_t*)malloc(b->l_data);
    memcpy(oldd, b->data, b->l_data);
    
    // clear out the old data and alloc the new amount
    free(b->data);
    b->data = (uint8_t*)calloc(new_size, sizeof(uint8_t)); 
    
    // add back the qname
    memcpy(b->data, oldd, b->core.l_qname); 
    
    // add in the new cigar
    memcpy(b->data + b->core.l_qname, new_cig, c.size()<<2);

    // add back the rest of the data
    memcpy(b->data + b->core.l_qname + (b->core.n_cigar<<2), oldd + old_seqaux_spot, old_seqaux_len);
    
    // update the sizes
    // >>1 shift is because only 4 bits needed per ATCGN base
    b->l_data = new_size; 
    b->core.n_cigar = c.size();
    
    free(oldd);
    free(new_cig);
  }

  BamRecord::BamRecord(const std::string& name, const std::string& seq, const std::string& ref, const GenomicRegion * gr) {

    StripedSmithWaterman::Aligner aligner;
    // Declares a default filter
    StripedSmithWaterman::Filter filter;
    // Declares an alignment that stores the result
    StripedSmithWaterman::Alignment alignment;
    // Aligns the seq to the ref
    aligner.Align(seq.c_str(), ref.c_str(), ref.size(), filter, &alignment);

    init();
    b->core.tid = gr->chr;
    b->core.pos = gr->pos1 + alignment.ref_begin + 1; // add to make it 1-indexed, not 0-indexed
    b->core.qual = 60; //alignment.sw_score;
    b->core.flag = 0;
    b->core.n_cigar = alignment.cigar.size();
    
    // set dumy mate
    b->core.mtid = -1;
    b->core.mpos = -1;
    b->core.isize = 0;

    // allocate all the data
    b->core.l_qname = name.length() + 1;
    b->core.l_qseq = seq.length(); //(seq.length()>>1) + seq.length() % 2; // 4-bit encoding
    b->l_data = b->core.l_qname + (b->core.n_cigar<<2) + ((b->core.l_qseq+1)>>1) + (b->core.l_qseq);
    b.get()->data = (uint8_t*)malloc(b.get()->l_data);

    // allocate the qname
    memcpy(b->data, name.c_str(), name.length() + 1);

    // allocate the cigar. 32 bits per elem (4 type, 28 length)
    uint32_t * cigr = bam_get_cigar(b);
    for (size_t i = 0; i < alignment.cigar.size(); ++i) {
      cigr[i] = alignment.cigar[i]; //Length << BAM_CIGAR_SHIFT | BAM_CMATCH;
    }

    // allocate the sequence
    uint8_t* m_bases = b->data + b->core.l_qname + (b->core.n_cigar<<2);
    
    // TODO move this out of bigger loop
    int slen = seq.length();
    for (int i = 0; i < slen; ++i) {
      // bad idea but works for now
      uint8_t base = 15;
      if (seq.at(i) == 'A')
	base = 1;
      else if (seq.at(i) == 'C')
	base = 2;
      else if (seq.at(i) == 'G')
	base = 4;
      else if (seq.at(i) == 'T')
	base = 8;
      
      m_bases[i >> 1] &= ~(0xF << ((~i & 1) << 2));   ///< zero out previous 4-bit base encoding
      m_bases[i >> 1] |= base << ((~i & 1) << 2);  ///< insert new 4-bit base encoding
      
    }

    // add in the actual alignment score
    AddIntTag("AS", alignment.sw_score);
      
  }

  void BamRecord::SmartAddTag(const std::string& tag, const std::string& val)
  {
    // get the old tag
    assert(tag.length());
    assert(val.length());
    std::string tmp;
    GetZTag(tag, tmp);
    if (!tmp.length()) 
      {
	AddZTag(tag, val);
	return;
      }
    
    // check that we don't have the tag delimiter in the stirng
    if (val.find(TAG_DELIMITER) != std::string::npos)
      std::cerr << "BamRecord::SmartAddTag -- Tag delimiter " << TAG_DELIMITER << " is in the value to be added. Compile with diff tag delimiter or change val" << std::endl;

    // append the tag
    tmp += TAG_DELIMITER + val;
    
    // remove the old tag
    RemoveTag(tag.c_str());
    
    // add the new one
    assert(tmp.length());
    AddZTag(tag, tmp);
  }

  void BamRecord::ClearSeqQualAndTags() {

    int new_size = b->core.l_qname + ((b)->core.n_cigar<<2);// + 1; ///* 0xff seq */ + 1 /* 0xff qual */;
    b->data = (uint8_t*)realloc(b->data, new_size);
    b->l_data = new_size;
    b->core.l_qseq = 0;
  }

  void BamRecord::SetSequence(const std::string& seq) {

    int new_size = b->l_data - ((b->core.l_qseq+1)>>1) - b->core.l_qseq + ((seq.length()+1)>>1) + seq.length();    
    int old_aux_spot = (b->core.n_cigar<<2) + b->core.l_qname + ((b->core.l_qseq + 1)>>1) + b->core.l_qseq;
    int old_aux_len = bam_get_l_aux(b); //(b->core.n_cigar<<2) + b->core.l_qname + ((b->core.l_qseq + 1)>>1) + b->core.l_qseq;

    // copy out all the old data
    uint8_t* oldd = (uint8_t*)malloc(b->l_data);
    memcpy(oldd, b->data, b->l_data);
    
    // clear out the old data and alloc the new amount
    free(b->data);
    b->data = (uint8_t*)calloc(new_size, sizeof(uint8_t)); 
    
    // add back the qname and cigar
    memcpy(b->data, oldd, b->core.l_qname + (b->core.n_cigar<<2));

    // update the sizes
    // >>1 shift is because only 4 bits needed per ATCGN base
    b->l_data = new_size; //b->l_data - ((b->core.l_qseq + 1)>>1) - b->core.l_qseq + ((seq.length()+1)>>1) + seq.length();
    b->core.l_qseq = seq.length();
    
    // allocate the sequence
    uint8_t* m_bases = b->data + b->core.l_qname + (b->core.n_cigar<<2);
    int slen = seq.length();

    for (int i = 0; i < slen; ++i) {
	
      // bad idea but works for now
      uint8_t base = 15;
      if (seq.at(i) == 'A')
	base = 1;
      else if (seq.at(i) == 'C')
	base = 2;
      else if (seq.at(i) == 'G')
	base = 4;
      else if (seq.at(i) == 'T')
	base = 8;
      
      m_bases[i >> 1] &= ~(0xF << ((~i & 1) << 2));   ///< zero out previous 4-bit base encoding
      m_bases[i >> 1] |= base << ((~i & 1) << 2);  ///< insert new 4-bit base encoding
    }

    // add in a NULL qual
    uint8_t* s = bam_get_qual(b);
    s[0] = 0xff;

    // add the aux data
    uint8_t* t = bam_get_aux(b);
    memcpy(t, oldd + old_aux_spot, old_aux_len);

    // reset the max size
    b->m_data = b->l_data;

    free(oldd); //just added
    
  }
  
  void BamRecord::SetQname(const std::string& n)
  {
    // copy out the non-qname data
    size_t nonq_len = b->l_data - b->core.l_qname;
    uint8_t* nonq = (uint8_t*)malloc(nonq_len);
    memcpy(nonq, b->data + b->core.l_qname, nonq_len);

    // clear the old data and alloc the new amount 
    free(b->data);
    b->data = (uint8_t*)calloc(nonq_len + n.length() + 1, 1);
    
    // add in the new qname
    memcpy(b->data, (uint8_t*)n.c_str(), n.length() + 1); // +1 for \0

    // update the sizes
    b->l_data = b->l_data - b->core.l_qname + n.length() + 1;
    b->core.l_qname = n.length() + 1;    
    
    // copy over the old data
    memcpy(b->data + b->core.l_qname, nonq, nonq_len);
    free(nonq);

    // reset the max size
    b->m_data = b->l_data;
  }

  void BamRecord::SetQualities(const std::string& n, int offset) {

    if (!n.empty() && n.length() != b->core.l_qseq)
      throw std::invalid_argument("New quality score should be same as seq length");
    
    // length of qual is always same as seq. If empty qual, just set first bit of qual to 0
    if (n.empty()) {
      uint8_t* r = bam_get_qual(b); 
      r[0] = 0;
      return;
    }

    char * q = strdup(n.data());
    for (size_t i = 0; i < n.length(); ++i)
      q[i] -= offset;
    memcpy(bam_get_qual(b), q, n.length()); // dont copy /0 terminator
    free(q);

  }

  double BamRecord::MeanPhred() const {

    if (b->core.l_qseq <= 0)
      return -1;

    double s = 0;
    uint8_t* p = bam_get_qual(b);
    for (int32_t i = 0; i < b->core.l_qseq; ++i)
      s += p[i];
    return s / b->core.l_qseq;
  }

  std::string BamRecord::QualitySequence() const {
    std::string seq;
    GetZTag("GV", seq);
    if (!seq.length()) 
      seq = Sequence();
    return seq;
  }

  std::ostream& operator<<(std::ostream& out, const BamRecord &r)
  {
    if (!r.b) {
      out << "empty read";
      return out;
    }
    out << bam_get_qname(r.b) << "\t" << r.b->core.flag
	<< "\t" << (r.b->core.tid+1) << "\t" << r.b->core.pos 
	<< "\t" << r.b->core.qual << "\t" << r.CigarString() 
	<< "\t" << (r.b->core.mtid+1) << "\t" << r.b->core.mpos << "\t" 
        << r.FullInsertSize() //r.b->core.isize 
	<< "\t" << r.Sequence() << "\t*" << std::endl;
    return out;
      
    
  }

  int32_t BamRecord::CountBWASecondaryAlignments() const 
  {
    int xp_count = 0;
    
    // xa tag
    std::string xar_s;
    GetZTag("XA", xar_s);
    if (xar_s.length()) {
      xp_count += std::count(xar_s.begin(), xar_s.end(), ';');
    }

    return xp_count;
    
  }

  int32_t BamRecord::CountBWAChimericAlignments() const 
  {
    int xp_count = 0;
    
    // sa tag (post bwa mem v0.7.5)
    std::string xar_s;
    GetZTag("SA", xar_s);
    if (xar_s.length()) 
      xp_count += std::count(xar_s.begin(), xar_s.end(), ';');

    // xp tag (pre bwa mem v0.7.5)
    std::string xpr_s;
    GetZTag("XP", xpr_s);
    if (xpr_s.length()) 
      xp_count += std::count(xpr_s.begin(), xpr_s.end(), ';');

    return xp_count;
    
  }

  int32_t BamRecord::CountNBases() const {
    uint8_t* p = bam_get_seq(b); 
    int32_t n = 0;
    for (int ww = 0; ww < b->core.l_qseq; ww++)
      if (bam_seqi(p,ww) == 15) 
	++n; 
    return n;
  }

  void BamRecord::QualityTrimmedSequence(int32_t qualTrim, int32_t& startpoint, int32_t& endpoint) const {

    endpoint = -1; //seq.length();
    startpoint = 0;
    int i = 0; 
    
    uint8_t * qual = bam_get_qual(b.get());
    
    // if there is no quality score, return whole thing
    if (qual[0] == 0xff) {
      startpoint = 0;
      return;
      
      //return Sequence();
    }
    
    // get the start point (loop forward)
    while(i < b->core.l_qseq) {
      int ps = qual[i];
      if (ps >= qualTrim) {
          startpoint = i;
          break;
	}
	++i;
    }

    // get the end point (loop backwards)
    i = b->core.l_qseq - 1; //seq.length() - 1;
    while(i >= 0) {

      int ps = qual[i];
      
      if (ps >= qualTrim) { //ps >= qualTrim) {
	endpoint = i + 1; // endpoint is one past edge
	break;
      }
      --i;
    }
  }

  void BamRecord::AddZTag(std::string tag, std::string val) {
    if (tag.empty() || val.empty())
      return;
    bam_aux_append(b.get(), tag.data(), 'Z', val.length()+1, (uint8_t*)val.c_str());
  }

  bool BamRecord::GetTag(const std::string& tag, std::string& s) const {

    if (GetZTag(tag, s))
      return true;

    int32_t t;
    if (GetIntTag(tag, t)) {
      std::stringstream ss; 
      ss << t;
      s = ss.str();
      return true;
    } 

    float f;
    if (GetFloatTag(tag, f)) {
      std::stringstream ss; 
      ss << f;
      s = ss.str();
      return true;
    }     

    return false;

  }

  bool BamRecord::GetZTag(const std::string& tag, std::string& s) const {
    uint8_t* p = bam_aux_get(b.get(),tag.c_str());
    if (!p)
      return false;

    int type = *p; 
    if (type != 'Z')
      return false;
    
    char* pp = bam_aux2Z(p);
    if (!pp) 
      return false;
    s = std::string(pp);
    return true;
  }

  
  // get a string tag that might be separted by "x"
  std::vector<std::string> BamRecord::GetSmartStringTag(const std::string& tag) const {
    
    std::vector<std::string> out;
    std::string tmp;
    GetZTag(tag, tmp);

    if (tmp.empty())
      return std::vector<std::string>();
    
    if (tmp.find(TAG_DELIMITER) != std::string::npos) {
      std::istringstream iss(tmp);
      std::string line;
      while (std::getline(iss, line, CTAG_DELIMITER)) {
	out.push_back(line);
      }
    } else {
      out.push_back(tmp);
    }
    
    assert(out.size());
    return out;
    
  }
  
  
  std::vector<int> BamRecord::GetSmartIntTag(const std::string& tag) const {
    
    std::vector<int> out;
    std::string tmp;
    
    GetZTag(tag, tmp);
    if (tmp.empty())
      return std::vector<int>();
    
    if (tmp.find(TAG_DELIMITER) != std::string::npos) {
      std::istringstream iss(tmp);
      std::string line;
      while (std::getline(iss, line, CTAG_DELIMITER))
	out.push_back(atoi(line.c_str())); 
    } else {
      out.push_back(atoi(tmp.c_str())); 
    }
    
    assert(out.size());
    return out;
    
  }

  std::vector<double> BamRecord::GetSmartDoubleTag(const std::string& tag) const {
    
    std::vector<double> out;
    std::string tmp;
    
    GetZTag(tag, tmp);
    if (tmp.empty())
      return std::vector<double>();
    
    if (tmp.find(TAG_DELIMITER) != std::string::npos) {
      std::istringstream iss(tmp);
      std::string line;
      while (std::getline(iss, line, CTAG_DELIMITER))
	out.push_back(std::atof(line.c_str())); 
    } else { // single entry
      out.push_back(std::atof(tmp.c_str())); 
    }
    
    assert(out.size());
    return out;
    
  }

  BamRecord::BamRecord(const std::string& name, const std::string& seq, const GenomicRegion * gr, const Cigar& cig) {

    // make sure cigar fits with sequence
    if (cig.NumQueryConsumed() != seq.length())
      throw std::invalid_argument("Sequence string length mismatches cigar consumed query bases");

    // make sure alignment fits
    if (cig.NumReferenceConsumed() != gr->Width())
      throw std::invalid_argument("Alignment position mismatches cigar consumed reference bases");

    init();
    b->core.tid = gr->chr;
    b->core.pos = gr->pos1; //gr->pos1 + 1;
    b->core.qual = 60;
    b->core.flag = 0;
    b->core.n_cigar = cig.size();
    
    // set dumy mate
    b->core.mtid = -1;
    b->core.mpos = -1;
    b->core.isize = 0;
      
    // if alignment is reverse, set it
    if (gr->strand == '-') // just choose this convention to reverse
      b->core.flag |= BAM_FREVERSE;
    
    // allocate all the data
    b->core.l_qname = name.length() + 1;
    b->core.l_qseq = seq.length(); //(seq.length()>>1) + seq.length() % 2; // 4-bit encoding
    b->l_data = b->core.l_qname + (b->core.n_cigar<<2) + ((b->core.l_qseq+1)>>1) + (b->core.l_qseq);
    b.get()->data = (uint8_t*)malloc(b.get()->l_data);
    
    // allocate the qname
    memcpy(b->data, name.c_str(), name.length() + 1);
      
    // allocate the cigar. 32 bits per elem (4 type, 28 length)
    uint32_t * cigr = bam_get_cigar(b);
    for (size_t i = 0; i < cig.size(); ++i)
      cigr[i] = cig[i].raw(); //Length << BAM_CIGAR_SHIFT | BAM_CMATCH;
    
    // allocate the sequence
    uint8_t* m_bases = b->data + b->core.l_qname + (b->core.n_cigar<<2);

      // TODO move this out of bigger loop
      int slen = seq.length();
      for (int i = 0; i < slen; ++i) {
	// bad idea but works for now
	uint8_t base = 15;
	if (seq.at(i) == 'A')
	  base = 1;
	else if (seq.at(i) == 'C')
	  base = 2;
	else if (seq.at(i) == 'G')
	  base = 4;
	else if (seq.at(i) == 'T')
	  base = 8;
	
	m_bases[i >> 1] &= ~(0xF << ((~i & 1) << 2));   ///< zero out previous 4-bit base encoding
	m_bases[i >> 1] |= base << ((~i & 1) << 2);  ///< insert new 4-bit base encoding
	
      }
  }
  

  CigarField::CigarField(char  t, uint32_t len) {
    int op = CigarCharToInt[(int)t];
    if (op < 0)
      throw std::invalid_argument("Cigar type must be one of MIDSHPN=X");      
    data = len << BAM_CIGAR_SHIFT;
    data = data | static_cast<uint32_t>(op);
  }

  std::ostream& operator<<(std::ostream& out, const CigarField& c) { 
    out << bam_cigar_oplen(c.data) << bam_cigar_opchr(c.data); 
    return out; 
  }


  std::ostream& operator<<(std::ostream& out, const Cigar& c) { 
    for (Cigar::const_iterator i = c.begin(); i != c.end(); ++i)
      out << *i;
    return out; 
  }


  Cigar::Cigar(const std::string& cig) {

    //Cigar tc;

    // get the ops (MIDSHPN)
    std::vector<char> ops;
    for (size_t i = 0; i < cig.length(); ++i)
      if (!isdigit(cig.at(i))) {
	ops.push_back(cig.at(i));
      }
    
    std::size_t prev = 0, pos;
    std::vector<std::string> lens;
    while ((pos = cig.find_first_of("MIDSHPNX", prev)) != std::string::npos) {
        if (pos > prev)
	  lens.push_back(cig.substr(prev, pos-prev));
        prev = pos+1;
      }
    if (prev < cig.length())
      lens.push_back(cig.substr(prev, std::string::npos));

    assert(ops.size() == lens.size());
    for (size_t i = 0; i < lens.size(); ++i) {
      add(CigarField(ops[i], std::atoi(lens[i].c_str())));
    }
    
    //return tc;

  }

  bool Cigar::operator==(const Cigar& c) const { 
     if (m_data.size() != c.size())
       return false;
     if (!m_data.size()) // both empty
       return true;
     for (size_t i = 0; i < m_data.size(); ++i)
       if (m_data[i].Type() != c[i].Type() || m_data[i].Length() != c[i].Length())
	 return false;
     return true;
  }


  int BamRecord::OverlappingCoverage(const BamRecord& r) const {
    
    uint32_t* c  = bam_get_cigar(b);
    uint32_t* c2 = bam_get_cigar(r.b);
    
    //uint8_t * cov1 = (uint8_t*)calloc(l > 0 ? l : b->core.l_qseq, sizeof(uint8_t));
    uint8_t * cov1 = (uint8_t*)calloc(GetCigar().NumQueryConsumed(), sizeof(uint8_t));
    size_t pos = 0;
    for (int k = 0; k < b->core.n_cigar; ++k) {
      if (bam_cigar_opchr(c[k]) == 'M')  // is match, so track locale
	for (size_t j = 0; j < bam_cigar_oplen(c[k]); ++j)
	  cov1[pos + j] = 1;
      if (bam_cigar_type(bam_cigar_op(c[k]))&1)  // consumes query, so move position
	pos = pos + bam_cigar_oplen(c[k]);
    }
    
    pos = 0;
    size_t ocov = 0; // overlapping coverage
    for (int k = 0; k < r.b->core.n_cigar; ++k) {
      if (bam_cigar_opchr(c2[k]) == 'M')  // is match, so track local
	for (size_t j = 0; j < bam_cigar_oplen(c2[k]); ++j)
	  if (cov1[pos+j]) // r is covered. Check again this too
	    ++ocov;
      if (bam_cigar_type(bam_cigar_op(c2[k]))&1)  // consumes query, so move position
	pos = pos + bam_cigar_oplen(c2[k]);
    }
    
    free(cov1);
    
    return ocov;
  }

  
}

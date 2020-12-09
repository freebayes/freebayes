#ifndef SEQLIB_BAM_RECORD_H
#define SEQLIB_BAM_RECORD_H

#include <stdint.h>
//#include <cstdint> //+11
#include <vector>
#include <iostream>
#include <sstream>
#include <cassert>
#include <algorithm>

extern "C" {
#include "htslib/htslib/hts.h"
#include "htslib/htslib/sam.h"
#include "htslib/htslib/bgzf.h"
#include "htslib/htslib/kstring.h"
#include "htslib/htslib/faidx.h"

}

#include "SeqLib/SeqLibUtils.h"
#include "SeqLib/GenomicRegion.h"
#include "SeqLib/UnalignedSequence.h"

static const char BASES[16] = {' ', 'A', 'C', ' ',
                               'G', ' ', ' ', ' ', 
                               'T', ' ', ' ', ' ', 
                               ' ', ' ', ' ', 'N'};

static std::string cigar_delimiters = "MIDNSHPX";

static const uint8_t CIGTAB[255] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                                    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                                    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                                    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                                    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                                    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                                    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                                    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                                    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
#define FRORIENTATION 0
#define FFORIENTATION 1
#define RFORIENTATION 2
#define RRORIENTATION 3
#define UDORIENTATION 4

namespace SeqLib {

/** Basic container for a single cigar operation
 *
 * Stores a single cigar element in a compact 32bit form (same as HTSlib).
 */
class CigarField {

  friend class Cigar;

 public:

  /** Construct the cigar op by type (MIDNSHP=X) and length 
   * @param t Cigar op (MIDNSHP=X)
   * @param l Cigar length
   * @exception Throws an invalid_argument if l <= 0 or invalid cigar op
   */
  CigarField(char t, uint32_t l); 

  /** Construct the cigar op from the raw sam.h uint32_t (first 4 bits op, last 28 len) */
  CigarField(uint32_t f) : data(f) {}

  /** Return the raw sam.h uint8_t cigar data */
  inline uint32_t raw() const { return data; }

  /** Print the cigar field (eg 35M) */
  friend std::ostream& operator<<(std::ostream& out, const CigarField& c);

  /** Return the cigar op type (one of MIDNSHPX) as a char */
  inline char Type() const { return bam_cigar_opchr(data); } 

  /** Return the raw sam.h uint8_t cigar type (bam_cigar_op(data)) */
  inline uint8_t RawType() const { return bam_cigar_op(data); } 

  /** Return the length of the cigar op (eg 35M returns 35) */
  inline uint32_t Length() const { return bam_cigar_oplen(data); } 

  /** Returns true if cigar op matches bases on the reference (MDN=X) */
  inline bool ConsumesReference() const { return bam_cigar_type(bam_cigar_op(data))&2;  }

  /** Returuns true cigar op matches bases on the query (MIS=X) */
  inline bool ConsumesQuery() const { return bam_cigar_type(bam_cigar_op(data))&1;  }

  /** Return whether two CigarField objects have same op and len */
  inline bool operator==(const CigarField& c) const { return c.Type() == Type() && c.Length() == Length(); }

  /** Return whether two CigarField objects have different op and/or len */
  inline bool operator!=(const CigarField& c) const { return !(c == *this); } 

 private:

  // first 4 bits hold op, last 28 hold len
  uint32_t data;
  
};

/** CIGAR for a single gapped alignment
 *
 * Constructed as a vector of CigarField objects. 
 */
 class Cigar {

 public:

   /** Construct an empty CIGAR */
   Cigar() {} 

   /** Construct from a CIGAR string 
    * @param cig CIGAR string, e.g. 54M46S
    */
   Cigar(const std::string& cig);

   typedef std::vector<CigarField>::iterator iterator; ///< Iterator for move between CigarField ops
   typedef std::vector<CigarField>::const_iterator const_iterator; ///< Iterator (const) for move between CigarField ops
   iterator begin() { return m_data.begin(); } ///< Iterator (aka std::vector<CigarField>.begin()
   iterator end()   { return m_data.end(); } ///< Iterator (aka std::vector<CigarField>.end()
   const_iterator begin() const { return m_data.begin(); } ///< Iterator (aka std::vector<CigarField>.begin()
   const_iterator end() const   { return m_data.end(); } ///< Iterator (aka std::vector<CigarField>.end() 

   /** Const reference to last cigar op */
   inline const CigarField& back() const { return m_data.back(); }

   /** Reference to last cigar op */
   inline CigarField& back() { return m_data.back(); }

   /** Const reference to first cigar op */
   inline const CigarField& front() const { return m_data.front(); }

   /** Reference to first cigar op */
   inline CigarField& front() { return m_data.front(); }

   /** Returns the number of cigar ops */
   inline size_t size() const { return m_data.size(); }

   /** Returns the i'th cigar op */
   inline CigarField& operator[](size_t i) { return m_data[i]; }

   /** Returns the i'th cigar op (const) */
   const CigarField& operator[](size_t i) const { return m_data[i]; }

   /** Return the sum of all of the lengths for all kinds */
   inline int TotalLength() const {
     int t = 0;
     for (Cigar::const_iterator c = m_data.begin(); c != m_data.end(); ++c)
       //for (auto& c : m_data)
       t += c->Length();
     return t;
   }

   /** Return the number of query-consumed bases */
   inline int NumQueryConsumed() const {
     int out = 0;
     for (Cigar::const_iterator c = m_data.begin(); c != m_data.end(); ++c)
       if (c->ConsumesQuery())
	 out += c->Length();
     return out;
   }

   /** Return the number of reference-consumed bases */
   inline int NumReferenceConsumed() const {
     int out = 0;
     //    for (auto& c : m_data)
     for (Cigar::const_iterator c = m_data.begin(); c != m_data.end(); ++c)
       if (c->ConsumesReference())
	 out += c->Length();
     return out;
   }

   /** Add a new cigar op */
   inline void add(const CigarField& c) { 
     m_data.push_back(c); 
   }

   /** Return whether two Cigar objects are equivalent */
   bool operator==(const Cigar& c) const;
   
   /** Return whether two Cigar objects are not equivalent */
   bool operator!=(const Cigar& c) const { return !(c == *this); }

  /** Print cigar string (eg 35M25S) */
  friend std::ostream& operator<<(std::ostream& out, const Cigar& c);
  
   
 private:
   
   std::vector<CigarField> m_data; // should make this simpler

 };

 typedef SeqHashMap<std::string, size_t> CigarMap;

/** Class to store and interact with a SAM alignment record
 *
 * HTSLibrary reads are stored in the bam1_t struct. Memory allocation
 * is taken care of by bam1_t init, and deallocation by destroy_bam1. This
 * class is a C++ interface that automatically takes care of memory management
 * for these C allocs/deallocs. The only member of BamRecord is a bam1_t object.
 * Alloc/dealloc is taken care of by the constructor and destructor.
 */
class BamRecord {

  friend class BLATWraper;
  friend class BWAWrapper;

 public:

  /** Construct a BamRecord manually from a name, sequence, cigar and location
   * @param name Name of the read
   * @param seq Sequence of the read (compsed of ACTG or N).
   * @param gr Location of the alignment
   * @param cig Cigar alignment
   * @exception Throws an invalid_argument exception if length of seq is not commensurate
   * with number of query-bases consumed in cigar. 
   * @exception Throws an invalid_argument exception if width of gr is not commensurate
   * with number of reference-bases consumed in cigar. 
   */
  BamRecord(const std::string& name, const std::string& seq, const GenomicRegion * gr, const Cigar& cig);
  
  /** Construct an empty BamRecord by calling bam_init1() 
   */
  void init();

  /** Check if a read is empty (not initialized)
   * @return true if read was not initialized with any values
   */
  bool isEmpty() const { return !b; }

  /** Explicitly pass a bam1_t to the BamRecord. 
   *
   * The BamRecord now controls the memory, and will delete at destruction
   * @param a An allocated bam1_t
   */
  void assign(bam1_t* a);

  /** Make a BamRecord with no memory allocated and a null header */
  BamRecord() {}

  /** BamRecord is aligned on reverse strand */
  inline bool ReverseFlag() const { return b ? ((b->core.flag&BAM_FREVERSE) != 0) : false; }

  /** BamRecord has mate aligned on reverse strand */
  inline bool MateReverseFlag() const { return b ? ((b->core.flag&BAM_FMREVERSE) != 0) : false; }

  /** BamRecord has is an interchromosomal alignment */
  inline bool Interchromosomal() const { return b ? b->core.tid != b->core.mtid && PairMappedFlag() : false; }

  /** BamRecord is a duplicate */
  inline bool DuplicateFlag() const { return b ? ((b->core.flag&BAM_FDUP) != 0) : false; }

  /** BamRecord is a secondary alignment */
  inline bool SecondaryFlag() const { return b ? ((b->core.flag&BAM_FSECONDARY) != 0) : false; }

  /** BamRecord is paired */
  inline bool PairedFlag() const { return b ? ((b->core.flag&BAM_FPAIRED) != 0) : false; }

  /** Get the relative pair orientations 
   * 
   * 0 - FR (RFORIENTATION) (lower pos read is Fwd strand, higher is reverse)
   * 1 - FF (FFORIENTATION)
   * 2 - RF (RFORIENTATION)
   * 3 - RR (RRORIENTATION)
   * 4 - Undefined (UDORIENTATION) (unpaired or one/both is unmapped)
   */
  inline int PairOrientation() const {
    if (!PairMappedFlag())
      return UDORIENTATION;
    else if ( (!ReverseFlag() && Position() <= MatePosition() &&  MateReverseFlag() ) || // read 1
	      (ReverseFlag()  && Position() >= MatePosition() && !MateReverseFlag() ) ) // read 2
      return FRORIENTATION;
    else if (!ReverseFlag() && !MateReverseFlag())
      return FFORIENTATION;
    else if (ReverseFlag() && MateReverseFlag())
      return RRORIENTATION;
    else if (   ( ReverseFlag() && Position() < MatePosition() && !MateReverseFlag()) ||
                (!ReverseFlag() && Position() > MatePosition() &&  MateReverseFlag()))
      return RFORIENTATION;
    assert(false);
  }
  
  /** BamRecord is failed QC */
  inline bool QCFailFlag() const { return b ? ((b->core.flag&BAM_FQCFAIL) != 0) : false; }

  /** BamRecord is supplementary alignment */
  inline bool SupplementaryFlag() const { return b ? ((b->core.flag&BAM_FSUPPLEMENTARY) != 0) : false; }

  /** BamRecord is mapped */
  inline bool MappedFlag() const { return b ? ((b->core.flag&BAM_FUNMAP) == 0) : false; }

  /** BamRecord mate is mapped */
  inline bool MateMappedFlag() const { return b ? ((b->core.flag&BAM_FMUNMAP) == 0) : false; }

  /** BamRecord is mapped and mate is mapped and in pair */
  inline bool PairMappedFlag() const { return b ? (!(b->core.flag&BAM_FMUNMAP) && !(b->core.flag&BAM_FUNMAP) && (b->core.flag&BAM_FPAIRED) ) : false; }

  /** BamRecord is mapped in proper pair */
  inline bool ProperPair() const { return b ? (b->core.flag&BAM_FPROPER_PAIR) : false;} 

  /** BamRecord has proper orientation (FR) */
  inline bool ProperOrientation() const { 
    if (!b)
      return false;
    
    // mate on diff chrom gets false
    if (b->core.tid != b->core.mtid)
      return false;

    // if FR return true
    if (b->core.pos < b->core.mpos) {
      return (b->core.flag&BAM_FREVERSE) == 0 && (b->core.flag&BAM_FMREVERSE) != 0 ? true : false;
    } else {
      return (b->core.flag&BAM_FREVERSE) == 0 && (b->core.flag&BAM_FMREVERSE) != 0 ? false : true;
    }
      
  }

  /** Count the total number of N bases in this sequence */
  int32_t CountNBases() const;

  /** Trim the sequence down by removing bases from ends with low quality scores. Stores the
   * trimmed sequence in the GV tag, but does not affect any other part of read.
   * @param qualTrim Minimal quality score, zero-based (eg # == 2)
   * @param startpoint Returns the new starting point for the sequence
   * @param endpoint Return the new ending point for the sequence
   */
  void QualityTrimmedSequence(int32_t qualTrim, int32_t& startpoint, int32_t& endpoint) const;

  /** Retrieve the quality trimmed seqeuence from QT tag if made. Otherwise return normal seq */
  std::string QualitySequence() const;

  /** Get the alignment position */
  inline int32_t Position() const { return b ? b->core.pos : -1; }

  /** Get the alignment position, including soft clips */
  int32_t PositionWithSClips() const;
  
  /** Get the alignment position of mate */
  inline int32_t MatePosition() const { return b ? b->core.mpos: -1; }

  /** Count the number of secondary alignments by looking at XA tag.
   * @note A secondary alignment is an alternative mapping. This may not
   * work for non-BWA aligners that may not place the XA tag.
   */
  int32_t CountBWASecondaryAlignments() const;

  /** Count the number of chimeric alignments by looking at XP and SA tags 
   * @note A secondary alignment is an alternative mapping. This may not
   * work for non-BWA aligners that may not place the XP/SA tags. BWA-MEM 
   * used the XP tag prior to v0.7.5, and SA aftewards.
   */
  int32_t CountBWAChimericAlignments() const;

  /** Get the end of the alignment */
  int32_t PositionEnd() const;

  /** Get the end of the alignment, including soft clips */
  int32_t PositionEndWithSClips() const;

  /** Get the end of the aligment mate pair */
  int32_t PositionEndMate() const;

  /** Get the chromosome ID of the read */
  inline int32_t ChrID() const { return b ? b->core.tid : -1; }
  
  /** Get the chrosome ID of the mate read */
  inline int32_t MateChrID() const { return b ? b->core.mtid : -1; }
  
  /** Get the mapping quality */
  inline int32_t MapQuality() const { return b ? b->core.qual : -1; }

  /** Set the qc fail flag on/off (true -> on) */
  inline void SetQCFail(bool f) { 
    if (f)
      b->core.flag |= BAM_FQCFAIL;
    else
      b->core.flag &= ~BAM_FQCFAIL;
  }

  /** Set the mapping quality */
  inline void SetMapQuality(int32_t m) { if (b) b->core.qual = m; }

  /** Set the chr id */
  inline void SetChrID(int32_t i) { b->core.tid = i; }

  /** Set the chr id of mate */
  inline void SetChrIDMate(int32_t i) { b->core.mtid = i; }
  
  /** Set the position of the mate read */
  inline void SetPositionMate(int32_t i) { b->core.mpos = i; }

  /** Set the pair mapped flag on/off (true -> on) */
  inline void SetPairMappedFlag(bool f) { 
    if (f)
      b->core.flag |= BAM_FPAIRED;
    else
      b->core.flag &= ~BAM_FPAIRED;
  }

  /** Set the mate reverse flag on/off (true -> on) */
  inline void SetMateReverseFlag(bool f) { 
    if (f)
      b->core.flag |= BAM_FMREVERSE;
    else
      b->core.flag &= ~BAM_FMREVERSE;
  }

  /** Get the number of cigar fields */
  inline int32_t CigarSize() const { return b ? b->core.n_cigar : -1; }
  
  /** Check if this read is first in pair */
  inline bool FirstFlag() const { return (b->core.flag&BAM_FREAD1); }
  
  /** Get the qname of this read as a string */
  inline std::string Qname() const { return std::string(bam_get_qname(b)); }
  
  /** Get the qname of this read as a char array */
  inline char* QnameChar() const { return bam_get_qname(b); }
  
  /** Get the full alignment flag for this read */
  inline uint32_t AlignmentFlag() const { return b->core.flag; }
  
  /** Get the insert size for this read */
  inline int32_t InsertSize() const { return b->core.isize; } 

  /** Get the read group, first from qname, then by RG tag 
   * @return empty string if no readgroup found
   */
  inline std::string ParseReadGroup() const {

    // try to get from RG tag first
    std::string RG;
    if (GetZTag("RG", RG))
      return RG;

    // try to get the read group tag from qname second
    std::string qn = Qname();
    size_t posr = qn.find(":", 0);
    return (posr != std::string::npos) ? qn.substr(0, posr) : "NA";
  }

  /** Get the insert size, absolute value, and always taking into account read length */
  inline int32_t FullInsertSize() const {

    if (b->core.tid != b->core.mtid || !PairMappedFlag())
      return 0;

    return std::abs(b->core.pos - b->core.mpos) + GetCigar().NumQueryConsumed();

  }
  
  /** Get the number of query bases of this read (aka length) */
  inline int32_t Length() const { return b->core.l_qseq; }
  
  /** Append a tag with new value, delimited by 'x' */
  void SmartAddTag(const std::string& tag, const std::string& val);
  
  /** Set the query name */
  void SetQname(const std::string& n);

  /** Set the quality scores 
   * @param n String of quality scores or empty string
   * @param offset Offset parameter for encoding (eg 33)
   * @exception Throws an invalid_argument if n is non-empty
   * and different length than sequence
   */
  void SetQualities(const std::string& n, int offset);

  /** Set the sequence name 
   * @param seq Sequence in upper-case (ACTGN) letters. 
   */
  void SetSequence(const std::string& seq);

  /** Set the cigar field explicitly 
   * @param c Cigar operation to set
   * @note Will not check if the cigar ops are consistent with 
   * the length of the sequence.
   */
  void SetCigar(const Cigar& c);

  /** Print a SAM-lite record for this alignment */
  friend std::ostream& operator<<(std::ostream& out, const BamRecord &r);

  /** Return read as a GenomicRegion */
  GenomicRegion AsGenomicRegion() const;

  /** Return mate read as a GenomicRegion */
  GenomicRegion AsGenomicRegionMate() const;

   /** Return the number of "aligned bases" in the same style as BamTools
    *
    * BamTools reports AlignedBases, which for example returns the literal strings (for diff CIGARs):
    * 3S5M - CTG
    * 5M - CTAGC
    * 3M1D3M - ATG-TGA 
    * 3M1I3M - ATGCTGA
    *
    * @return The number of M, D, X, = and I bases
    */
  inline int NumAlignedBases() const {
    int out = 0;
    uint32_t* c = bam_get_cigar(b);
    for (size_t i = 0; i < b->core.n_cigar; i++) 
      if (bam_cigar_opchr(c[i]) == 'M' || 
	  bam_cigar_opchr(c[i]) == 'I' || 
	  bam_cigar_opchr(c[i]) == '=' || 
	  bam_cigar_opchr(c[i]) == 'X' || 
	  bam_cigar_opchr(c[i]) == 'D')
	out += bam_cigar_oplen(c[i]);
    return out;
  }
  

  /** Return the max single insertion size on this cigar */
  inline uint32_t MaxInsertionBases() const {
    uint32_t* c = bam_get_cigar(b);
    uint32_t imax = 0;
    for (size_t i = 0; i < b->core.n_cigar; i++) 
      if (bam_cigar_opchr(c[i]) == 'I')
	imax = std::max(bam_cigar_oplen(c[i]), imax);
    return imax;
  }

  /** Return the max single deletion size on this cigar */
  inline uint32_t MaxDeletionBases() const {
    uint32_t* c = bam_get_cigar(b);
    uint32_t dmax = 0;
    for (size_t i = 0; i < b->core.n_cigar; i++) 
      if (bam_cigar_opchr(c[i]) == 'D')
	dmax = std::max(bam_cigar_oplen(c[i]), dmax);
    return dmax;
  }

  /** Get the number of matched bases in this alignment */
  inline uint32_t NumMatchBases() const {
    uint32_t* c = bam_get_cigar(b);
    uint32_t dmax = 0;
    for (size_t i = 0; i < b->core.n_cigar; i++) 
      if (bam_cigar_opchr(c[i]) == 'M')
	dmax += bam_cigar_oplen(c[i]);
    return dmax;
  }


  /** Retrieve the CIGAR as a more managable Cigar structure */
  Cigar GetCigar() const {
    uint32_t* c = bam_get_cigar(b);
    Cigar cig;
    for (size_t k = 0; k < b->core.n_cigar; ++k) {
      cig.add(CigarField(c[k]));
    }
    return cig;
  }

  /** Retrieve the inverse of the CIGAR as a more managable Cigar structure */
  Cigar GetReverseCigar() const {
    uint32_t* c = bam_get_cigar(b);
    Cigar cig;
    for (int k = b->core.n_cigar - 1; k >= 0; --k) 
      cig.add(CigarField(c[k]));
    return cig;
  }

  /** Remove the sequence, quality and alignment tags. 
   * Make a more compact alignment stucture, without the string data
   */
  void ClearSeqQualAndTags();

  /** Retrieve the sequence of this read as a string (ACTGN) */
  std::string Sequence() const;

  /** Return the mean quality score 
   */
  double MeanPhred() const;

  /** Performa a Smith-Waterman alignment between two strings
   * @param name Name of the query sequence to align
   * @param seq Sequence (ACTGN) of the query string
   * @param ref Sequence (ACTGN) of the reference string
   * @param gr Location of the reference string. The alignment record after Smith-Waterman alignment
   * will be relative to this location.
   */
  BamRecord(const std::string& name, const std::string& seq, const std::string& ref, const GenomicRegion * gr);

  /** Get the quality scores of this read as a string 
   * @param offset Encoding offset for phred quality scores. Default 33
   * @return Qualties scores after converting offset. If first char is empty, returns empty string
   */
  inline std::string Qualities(int offset = 33) const { 
    uint8_t * p = bam_get_qual(b);
    if (!p)
      return std::string();
    //if (!p[0])
    //  return std::string();
    std::string out(b->core.l_qseq, ' ');
    for (int32_t i = 0; i < b->core.l_qseq; ++i) 
      out[i] = (char)(p[i] + offset);
    return out;
  }

  /** Get the start of the alignment on the read, by removing soft-clips
   * Do this in the reverse orientation though.
   */
  inline int32_t AlignmentPositionReverse() const {
    uint32_t* c = bam_get_cigar(b);
    int32_t p = 0;
    for (int32_t i = b->core.n_cigar - 1; i >= 0; --i) {
      if ( (bam_cigar_opchr(c[i]) == 'S') || (bam_cigar_opchr(c[i]) == 'H'))
	p += bam_cigar_oplen(c[i]);
      else // not a clip, so stop counting
	break;
    }
    return p;
  }
  
  /** Get the end of the alignment on the read, by removing soft-clips
   * Do this in the reverse orientation though.
   */
  inline int32_t AlignmentEndPositionReverse() const {
    uint32_t* c = bam_get_cigar(b);
    int32_t p = 0;
    for (size_t i = 0; i < b->core.n_cigar; ++i) { // loop from the end
      if ( (bam_cigar_opchr(c[i]) == 'S') || (bam_cigar_opchr(c[i]) == 'H'))
	p += bam_cigar_oplen(c[i]);
      else // not a clip, so stop counting
	break;
    }
    return (b->core.l_qseq - p);
  }


  /** Get the start of the alignment on the read, by removing soft-clips
   */
  inline int32_t AlignmentPosition() const {
    uint32_t* c = bam_get_cigar(b);
    int32_t p = 0;
    for (size_t i = 0; i < b->core.n_cigar; ++i) {
      if (bam_cigar_opchr(c[i]) == 'S')
	p += bam_cigar_oplen(c[i]);
      else if (bam_cigar_opchr(c[i]) != 'H') 
	break;
    }
    return p;
  }
  
  /** Get the end of the alignment on the read, by removing soft-clips
   */
  inline int32_t AlignmentEndPosition() const {
    uint32_t* c = bam_get_cigar(b);
    int32_t p = 0;
    for (int32_t i = b->core.n_cigar - 1; i >= 0; --i) { // loop from the end
      if ( (bam_cigar_opchr(c[i]) == 'S') || (bam_cigar_opchr(c[i]) == 'H'))
	p += bam_cigar_oplen(c[i]);
      else // not a clip, so stop counting
	break;
    }
    return (b->core.l_qseq - p);
  }

  /** Get the number of soft clipped bases */
  inline int32_t NumSoftClip() const {
      int32_t p = 0;
      uint32_t* c = bam_get_cigar(b);
      for (size_t i = 0; i < b->core.n_cigar; ++i)
	if (bam_cigar_opchr(c[i]) == 'S')
	  p += bam_cigar_oplen(c[i]);
      return p;
    }

  /** Get the number of hard clipped bases */
  inline int32_t NumHardClip() const {
      int32_t p = 0;
      uint32_t* c = bam_get_cigar(b);
      for (size_t i = 0; i < b->core.n_cigar; ++i) 
	if (bam_cigar_opchr(c[i]) == 'H')
	  p += bam_cigar_oplen(c[i]);
      return p;
    }


  /** Get the number of clipped bases (hard clipped and soft clipped) */
  inline int32_t NumClip() const {
    int32_t p = 0;
    uint32_t* c = bam_get_cigar(b);
    for (size_t i = 0; i < b->core.n_cigar; ++i)
      if ( (bam_cigar_opchr(c[i]) == 'S') || (bam_cigar_opchr(c[i]) == 'H') )
	p += bam_cigar_oplen(c[i]);
    return p;
  }
  
  /** Get a string (Z) tag 
   * @param tag Name of the tag. eg "XP"
   * @param s The string to be filled in with the tag information
   * @return Returns true if the tag is present, even if empty. Return false if no tag or not a Z tag.
   */
  bool GetZTag(const std::string& tag, std::string& s) const;
  
  /** Get a string of either Z, f or i type. Useful if tag type not known at compile time.
   * @param tag Name of the tag. eg "XP"
   * @param s The string to be filled in with the tag information
   * @return Returns true if the tag is present and is either Z or i, even if empty. Return false if no tag or not Z or i.
   */  
  bool GetTag(const std::string& tag, std::string& s) const;
  
  /** Get a vector of type int from a Z tag delimited by "^"
   * Smart-tags allow one to store vectors of strings, ints or doubles in the alignment tags, and
   * do not require an additional data structure on top of bseq1_t. 
   * @param tag Name of the tag eg "AL"
   * @return A vector of ints, retrieved from the x delimited Z tag
   * @exception Throws an invalid_argument if cannot convert delimited field val to int
   */
  std::vector<int> GetSmartIntTag(const std::string& tag) const;

  /** Get a vector of type double from a Z tag delimited by "x"
   * Smart-tags allow one to store vectors of string, ints or doubles in the alignment tags, and
   * do not require an additional data structure on top of bseq1_t. 
   * @param tag Name of the tag eg "AL"
   * @return A vector of double elems, retrieved from the "^" delimited Z tag
   * @exception Throws an invalid_argument if cannot convert delimited field val to double
   */
  std::vector<double> GetSmartDoubleTag(const std::string& tag) const;

  /** Get a vector of strings from a Z tag delimited by "^"
   * Smart-tags allow one to store vectors of strings, ints or doubles in the alignment tags, and
   * do not require an additional data structure on top of bseq1_t. 
   * @param tag Name of the tag eg "CN"
   * @return A vector of strngs, retrieved from the x delimited Z tag
   */
  std::vector<std::string> GetSmartStringTag(const std::string& tag) const;

  /** Get an int (i) tag 
   * @param tag Name of the tag. eg "XP"
   * @param t Value to be filled in with the tag value.
   * @return Return true if the tag exists.
   */
  inline bool GetIntTag(const std::string& tag, int32_t& t) const {
    uint8_t* p = bam_aux_get(b.get(),tag.c_str());
    if (!p)
      return false;
    t = bam_aux2i(p);
    int type = *p++;
    if (!(type == 'i' || type == 'C' || type=='S' || type=='s' || type =='I' || type=='c'))
      return false;

    return true;
  }

  /** Get a float (f) tag 
   * @param tag Name of the tag. eg "AS"
   * @param t Value to be filled in with the tag value.
   * @return Return true if the tag exists.
   */
  inline bool GetFloatTag(const std::string& tag, float& t) const {
    uint8_t* p = bam_aux_get(b.get(),tag.c_str());
    if (!p)
      return false;

    t = bam_aux2f(p);
    int type = *p;
    type = *p++;
    if (!(type == 'f' || type == 'd'))
      return false;

    return true;
  }

  /** Add a string (Z) tag
   * @param tag Name of the tag. eg "XP"
   * @param val Value for the tag
   */
  void AddZTag(std::string tag, std::string val);

  /** Add an int (i) tag
   * @param tag Name of the tag. eg "XP"
   * @param val Value for the tag
   */
  inline void AddIntTag(const std::string& tag, int32_t val) {
    bam_aux_append(b.get(), tag.data(), 'i', 4, (uint8_t*)&val);
  }

  /** Set the chr id number 
   * @param id Chromosome id. Typically is 0 for chr1, etc
   */
  inline void SetID(int32_t id) {
    b->core.tid = id;
  }
  
  /** Set the alignment start position
   * @param pos Alignment start position
   */
  inline void SetPosition(int32_t pos) {
    b->core.pos = pos;
  }

  /** Convert CIGAR to a string
   */
  inline std::string CigarString() const {
    std::stringstream cig;
    uint32_t* c = bam_get_cigar(b);
    for (size_t k = 0; k < b->core.n_cigar; ++k)
      cig << bam_cigar_oplen(c[k]) << "MIDNSHP=XB"[c[k]&BAM_CIGAR_MASK];
    return cig.str();
  }
  
  /** Return a human readable chromosome name assuming chr is indexed
   * from 0 (eg id 0 return "1")
   * @note This is a quick convienence function, and is not robust to non-numbered
   * chromosomes (eg chrX becomes 23). For accurate string representation of 
   * any chromosomes, use the full ChrName with BamHeader input.
   */
  inline std::string ChrName() const {
    std::stringstream ss;
    ss << (b->core.tid + 1);

    return ss.str();
    //return std::to_string(b->core.tid + 1); //c++11
  }

  /** Retrieve the human readable chromosome name. 
   * @param h Dictionary for chr name lookup. If it is empty, assumes this is chr1 based reference.
   * @exception Throws an out_of_range exception if chr id is not in dictionary
   * @return Empty string if chr id < 0, otherwise chromosome name from dictionary.
   */
  inline std::string ChrName(const SeqLib::BamHeader& h) const {
    if (b->core.tid < 0)
      return std::string();

    if (!h.isEmpty())
      return h.IDtoName(b->core.tid);

    // c++98    
    std::stringstream ss;
    ss << b->core.tid;
    
    // no header, assume zero based
    return ss.str(); //std::to_string(b->core.tid + 1);
    
  }

  /** Return a short description (chr:pos) of this read */
  inline std::string Brief() const {
    //if (!h)
    // c++11
    //  return(std::to_string(b->core.tid + 1) + ":" + AddCommas<int32_t>(b->core.pos) + "(" + ((b->core.flag&BAM_FREVERSE) != 0 ? "+" : "-") + ")");
    // c++98
    std::stringstream ss;
    ss << (b->core.tid + 1) << ":" << AddCommas(b->core.pos) << "(" << ((b->core.flag&BAM_FREVERSE) != 0 ? "+" : "-") << ")";
    return ss.str();
    //else
    // return(std::string(h->target_name[b->core.tid]) + ":" + AddCommas<int32_t>(b->core.pos) + "(" + ((b->core.flag&BAM_FREVERSE) != 0 ? "+" : "-") + ")");      
  }

  /** Return a short description (chr:pos) of this read's mate */
  inline std::string BriefMate() const {
    //if (!h)
    // c++11
    // return(std::to_string(b->core.mtid + 1) + ":" + AddCommas<int32_t>(b->core.mpos) + "(" + ((b->core.flag&BAM_FMREVERSE) != 0 ? "+" : "-") + ")");
    std::stringstream ss;
    ss << (b->core.mtid + 1) << ":" << AddCommas(b->core.mpos) << "(" << ((b->core.flag&BAM_FMREVERSE) != 0 ? "+" : "-") << ")";
    return ss.str();
    //else
    //  return(std::string(h->target_name[b->core.mtid]) + ":" + AddCommas<int32_t>(b->core.mpos) + "(" + ((b->core.flag&BAM_FMREVERSE) != 0 ? "+" : "-") + ")");      
  }

  /** Strip a particular alignment tag 
   * @param tag Tag to remove
   */
  inline void RemoveTag(const char* tag) {
    uint8_t* p = bam_aux_get(b.get(), tag);
    if (p)
      bam_aux_del(b.get(), p);
  }

  /** Strip all of the alignment tags */
  inline void RemoveAllTags() {
    size_t keep = (b->core.n_cigar<<2) + b->core.l_qname + ((b->core.l_qseq + 1)>>1) + b->core.l_qseq;
    b->data = (uint8_t*)realloc(b->data, keep); // free the end, which has aux data
    b->l_data = keep;
    b->m_data = b->l_data;
  }

  /** Return the raw pointer */
  inline bam1_t* raw() const { return b.get(); }

  /** Return the number of bases on the query that are covered by a match (M) on both reads 
   * This is for tracking overlapping coverage on the reads, regardless of their alignment locations.
   * For instance, two reads with 101M will have overlapping coverage of 101, regardless of alignment location.
   * A read with 50S50M and 50M50S will have 0 overlapping coverage.
   */
  int OverlappingCoverage(const BamRecord& r) const;
  
  /** Return the shared pointer */
  SeqPointer<bam1_t> shared_pointer() const { return b; }

  protected:
  
  SeqPointer<bam1_t> b; // bam1_t shared pointer

};

 typedef std::vector<BamRecord> BamRecordVector; ///< Store a vector of alignment records
 
 typedef std::vector<BamRecordVector> BamRecordClusterVector; ///< Store a vector of alignment vectors

 /** @brief Sort methods for alignment records
  */
 namespace BamRecordSort {

   /** @brief Sort by read position 
    */
   struct ByReadPosition
   {
     bool operator()( const BamRecord& lx, const BamRecord& rx ) const {
       return (lx.ChrID() < rx.ChrID()) || (lx.ChrID() == rx.ChrID() && lx.Position() < rx.Position());
     }
   };

   /** @brief Sort by mate position 
    */
   struct ByMatePosition
   {
     bool operator()( const BamRecord& lx, const BamRecord& rx ) const {
       return (lx.MateChrID() < rx.MateChrID()) || (lx.MateChrID() == rx.MateChrID() && lx.MatePosition() < rx.MatePosition());
     }
   };

}

}
#endif

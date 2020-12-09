#ifndef SEQLIB_GENOMIC_REGION_H__
#define SEQLIB_GENOMIC_REGION_H__

#include <vector>
#include <iostream>
#include <stdint.h>
#include <utility>
#include <list>
#include <cstring>

#include "SeqLib/SeqLibCommon.h"
#include "SeqLib/SeqLibUtils.h"
#include "SeqLib/BamHeader.h"

namespace SeqLib {

  /** @brief Container for an interval on the genome 
   */
class GenomicRegion {

  template<typename T> friend class GenomicRegionCollection;
  
 public:

  /** Construct an "empty" GenomicRegion at (chr -1), pos 0, width = 1
   */
 GenomicRegion() : chr(-1), pos1(0), pos2(0), strand('*') {};

  /** Construct a GenomicRegion from another 
   * @param gr A GenomicRegion to copy
   */
 GenomicRegion(const GenomicRegion& gr) : chr(gr.chr), pos1(gr.pos1), pos2(gr.pos2), strand(gr.strand) {}

  /** Construct a GenomicRegion at a specific start and end location 
   * @param t_chr Chromosome id  (chr1 = 0, etc)
   * @param t_pos1 Start position
   * @param t_pos2 End position. Must be >= start position.
   * @param t_strand +, -, or * (default is *)
   * @exception throws an invalid_argument exception if pos2 < pos1
   * @exception throws an invalid_argument exception if char not one of +, - , *
  */
  GenomicRegion(int32_t t_chr, int32_t t_pos1, int32_t t_pos2, char t_strand = '*');

  /** Construct a GenomicRegion from a set of strings 
   * @param tchr Chromosome name
   * @param tpos1 Position 1
   * @param tpos2 Position 2
   * @param hdr Header to be used as sequence dictionary to convert chromosome name to id
   * @exception Throws an invalid_argument if cannot convert string to int
   * @exception Throws an out_of_range if number if greater than int32_t max
   * @note If an empty BamHeader is provided, will try to guess chromosome id.
   * eg "1" -> 0, "X" -> 22, "chr2" -> 1.
   */
  GenomicRegion(const std::string& tchr, const std::string& tpos1, const std::string& tpos2, const BamHeader& hdr);

  /** Construct a GenomicRegion from a samtools style region string.
   *
   * This calls the samtools-like parser, which accepts in form "chr7:10,000-11,100".
   * Note that this requires that a BamHeader be provided as well 
   * to convert the text representation of the chr to the id number.
   * @param reg Samtools-style string (e.g. "1:1,000,000-2,000,000") or single chr
   * @param hdr Pointer to BAM header that will be used to convert chr string to ref id
   * @exception throws an invalid_argument exception if cannot parse correctly
   */
  GenomicRegion(const std::string& reg, const BamHeader& hdr);

  /** Return a string representation of just the first base-pair 
   * @param h BamHeader used as a lookup table for id to string (i.e. chr name)
   * e.g. 1:10,000
   */
  std::string PointString(const BamHeader& h) const;

  // Randomize the position of this GenomicRegion on the genome
  // 
  // Creates a GenomicRegion with pos1 = pos2. Simulates a random value
  // with val <= genome_size_XY and then converts to GenomicRegion
  // @note Seed is set before-hand at any time with srand
  //
  //void Random();

  /** Check if the GenomicRegion is empty (aka chr -1 and pos1=pos2=0)   */
  bool IsEmpty() const;

  /** Find the absolute distance between start of two GenomicRegion objects 
   * 
   * If chr1 != chr2, then -1 is returned
   * @param gr GenomicRegion object to compare with
   */
  int32_t DistanceBetweenStarts(const GenomicRegion &gr) const;

  /** Find the absolute distance between ends of two GenomicRegion objects 
   * 
   * If chr1 != chr2, then -1 is returned
   * @param gr GenomicRegion object to compare with
   */
  int32_t DistanceBetweenEnds(const GenomicRegion &gr) const;

  /** Returns true if a.chr < b.chr or a.pos1 < a.pos1 if on same chrome, or if a.pos2 < b.pos2 if same chrom and same pos1 */
  bool operator < (const GenomicRegion& b) const;

  /** Returns true if a.chr > b.chr or a.pos1 > a.pos1 if on same chrome, or if a.pos2 > b.pos2 if same chrom and same pos1 */
  bool operator > (const GenomicRegion& b) const;

  /** Returns true if chr, pos1, pos2. No strand consideration */
  bool operator==(const GenomicRegion& b) const;

  /** Returns opposite of == */
  bool operator!=(const GenomicRegion& b) const;

  /** Returns true if < or == . No strand consideration */
  bool operator<=(const GenomicRegion &b) const;

  /** Returns true if > or == . No strand consideration */
  bool operator>=(const GenomicRegion &b) const;
  
  /** Check if the GenomicRegion has a complete or partial overlap
   * If the argument contains the calling object, returns 3
   * If the argument is contained in the calling object, returns 2
   * If the argument overlaps partially the calling object, returns 1
   * If the argument and calling object do not overlap, returns 0
   * @param gr GenomicRegion to compare against
   */
  int GetOverlap(const GenomicRegion& gr) const;

  /** Print with chr ID bumped up by one to make eg ID 0
   * print as "1"
   */
  friend std::ostream& operator<<(std::ostream& out, const GenomicRegion& gr);

  /** Print the chromosome with the correct chromosome name as store in the header 
   * @param h BamHeader that stores the lookup between chr id and chr name string
   */
  std::string ToString(const BamHeader& h) const;

  /** Extract the chromosome name as a string 
   * @param h BamHeader to serve as sequence dictionary
   * @exception throws an out_of_range exception if ref id >= h->n_targets
   */
  std::string ChrName(const BamHeader& h) const;

  /** Pad the object to make larger or smaller
   * @param pad Amount to pad by.
   * @exception throws an out_of_bounds if for pad < -width/2
   */
  void Pad(int32_t pad);

  /** Return the width (inclusive)
   * @note Width is inclusive, so that if pos1=1 and pos2=2, width is 2
   */
  int Width() const;

  int32_t chr; ///< Chromosome ID

  int32_t pos1; ///< Start position

  int32_t pos2; ///< End Position

  char strand; ///< Strand. Should be one of *, -, +

 private:

  // Convert a chromosome number to a string using default ordering (1-Y)
  // Assumes a 1-based ordering (1, ...), not zero-based.
  // e.g. chrToString(10) return "11"
  // @param ref Reference ID to convert
  // @exception throws an invalid_argument exception if ref < 0
  std::string chrToString(int32_t ref) const;


};

typedef std::vector<GenomicRegion> GenomicRegionVector;

}


#endif

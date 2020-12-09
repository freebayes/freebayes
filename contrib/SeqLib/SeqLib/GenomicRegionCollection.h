#ifndef SWAP_GENOMIC_REGION_COLLECTION_H
#define SWAP_GENOMIC_REGION_COLLECTION_H

#include <vector>
#include <string>
#include <cstdlib>
#include <list>

#include "SeqLib/IntervalTree.h"
#include "SeqLib/GenomicRegionCollection.h"
#include "SeqLib/BamRecord.h"

namespace SeqLib {

  /** Simple structure to store overlap results 
   */
  typedef std::pair<size_t, size_t> OverlapResult;

/** Class to store vector of intervals on the genome */
typedef TInterval<int32_t> GenomicInterval;
typedef SeqHashMap<int, std::vector<GenomicInterval> > GenomicIntervalMap;
typedef TIntervalTree<int32_t> GenomicIntervalTree;
typedef SeqHashMap<int, GenomicIntervalTree> GenomicIntervalTreeMap;
typedef std::vector<GenomicInterval> GenomicIntervalVector;

  /** @brief Template class to store / query a collection of genomic intervals
   *
   * Can hold a collection of GenomicRegion objects, or any object whose
   * class is a child of GenomicRegion. Contains an implementation of an
   * interval tree (as provided by Erik Garrison) for fast interval queries.
   */
template<typename T=GenomicRegion>
class GenomicRegionCollection {

 public:

  /** Construct an empty GenomicRegionCollection 
   */
 GenomicRegionCollection();

 ~GenomicRegionCollection();
 
  /** Construct from a plain vector of GenomicRegion objects
   */
  GenomicRegionCollection(std::vector<T>& vec);

  /** Construct from a single GenomicRegion
   */
  GenomicRegionCollection(const T& gr);

  /** Construct from a vector of reads
   * 
   * @note See BamRecord::AsGenomicRegion 
   */
  GenomicRegionCollection(const BamRecordVector& brv);

  /** Construct a GenomicRegionCollection with overlapping intervals
   * 
   * @param width Desired bin width
   * @param ovlp Amount that the bins should overlap
   * @param gr GenomicRegion to divide into smaller overlapping bins
   */
  GenomicRegionCollection(int width, int ovlp, const T &gr);

 /** Construct a tiled set of intervals across a genome
  *
  * @param width Width of each interval tile
  * @param ovlp Amount of overlap between neighboring tiles
  * @param h Set of chromosomes and their lengths to build the tile on
  */
 GenomicRegionCollection(int width, int ovlp, const HeaderSequenceVector& h);

 // Read in a MuTect call-stats file and adds to GenomicRegionCollection object.
   //
   // Reads a MuTect call-stats file and imports only
   // lines with KEEP marked. 
   // @param file Path to call-stats file
   // @param pad Amount to pad intervals by
   // @return True if file was succesfully read
   //
   //bool ReadMuTect(const std::string &file, const SeqLib::BamHeader& hdr);

  /** Read in a BED file and adds to GenomicRegionCollection object
   * @param file Path to BED file
   * @param hdr Dictionary for converting chromosome strings in BED file to chr indicies
   * @return True if file was succesfully read
   */
   bool ReadBED(const std::string &file, const SeqLib::BamHeader& hdr);

  /** Read in a VCF file and adds to GenomicRegionCollection object
   * @param file Path to VCF file. All elements will be width = 1 (just read start point)
   * @param hdr Dictionary for converting chromosome strings in BED file to chr indicies
   */
  bool ReadVCF(const std::string &file, const SeqLib::BamHeader& hdr);

  /** Shuffle the order of the intervals */
 void Shuffle();

  /** Read in a text file (can be gzipped) and construct a GenomicRegionCollection
   *
   * This function will automatically detect which file type is being input:
   * -- ends in .vcf -> readVCFfile
   * -- ends in .bed -> readBEDfile
   * -- contains ':' -> Assumes single samtools-style region (eg 1:100-100)
   * The values are appended to existing vector of GenomicRegion objects
   * @param file Text file to read and store intervals
   * @param hdr BamHeader to serve as dictionary for chromosomes
   */
   GenomicRegionCollection(const std::string &file, const BamHeader& hdr);

  /** Create the set of interval trees (one per chromosome) 
   *
   * A GenomicIntervalTreeMap is an unordered_map of GenomicIntervalTrees for 
   * each chromosome. A GenomicIntervalTree is an interval tree on the ranges
   * defined by the genomic interval, with cargo set at the same GenomicRegion object.
   */
  void CreateTreeMap();
  
  /** Reduces the GenomicRegion objects to minimal set by merging overlapping intervals
   * @note This will merge intervals that touch. eg [4,6] and [6,8]
   * @note This will also call CreateTreeMap() at end to re-create the interval tree
   */
  void MergeOverlappingIntervals();

  /** Return the number of GenomicRegions stored 
   */
  size_t size() const { return m_grv->size(); }

  /** Add a new GenomicRegion (or child of) to end
   */
 void add(const T& g) { m_grv->push_back(g); /*createTreeMap();*/ }

  /** Is this object empty?
   */
  bool IsEmpty() const { return !m_grv->size(); }

  /** Clear out all of the GenomicRegion objects
   */
  void clear() { m_grv->clear(); 
		 m_tree->clear(); 
		 idx = 0;
  }

 /** Get the number of trees (eg number of chromosomes, each with own tree */
 int NumTree() const { return m_tree->size(); }

 /** Get the IDs of all intervals that overlap with a query range
  *
  * The IDs are created during CreateTreeMap, and are the position of the 
  * the individual intervals from the tree, in genomic order. e.g the first
  * interval on chromosome 1 gets 0, the next one gets 1, etc.
  * The returned IDs can then be used as lookups with [], as long as the 
  * collection is not altered in between.
  * @param gr Query range to check overlaps against
  * @param ignore_strand Should strandedness be ignore when doing overlaps 
  * @return A vector of IDs of intervals in this collection that overlap with gr
  */
 template<class K>
 std::vector<int> FindOverlappedIntervals(const K& gr, bool ignore_strand) const;

 /** Get a const pointer to the genomic interval tree map */
 const GenomicIntervalTreeMap* GetTree() const { return m_tree.get(); }

  /** Retrieve a GenomicRegion at given index. 
   * 
   * Note that this does not move the idx iterator, which is 
   * used to loop through all the regions. Throws an exception
   * if the index is out of bounds.
   * @return GenomicRegion pointed to by index i
   */
  const T& at(size_t i) const;

  /** Find overlaps between this vector and input GenomicRegion.
   *
   * Requires that the GenomicIntervalTreeMap have been created first
   * @param gr Region to test
   * @return Number of overlapping elements in this GenomicRegionCollection
   */
 size_t CountOverlaps(const T &gr) const;

 /** Test if two intervals overlap the same element in the collection
  */
 template<class K>
 bool OverlapSameInterval(const K &gr1, const K &gr2) const;

 /** Count the number of intervals in the collection contained in this range */
 size_t CountContained(const T &gr);

 /** Return the overlaps between the collection and the query collection
  * @param subject Subject collection of intervals
  * @param query_id Indices of the queries that have an overlap. Will be same size as output and subject_id and in same order
  * @param subject_id Indices of the subject that have an overlap. Will be same size as output and query_id and in same order
  * @param ignore_strand If true, won't exclude overlap if on different strand
  * @return A collection of overlapping intervals from this collection, trimmed to be contained
  * @exception Throws a logic_error if this tree is non-empty, but the interval tree has not been made with 
  * CreateTreeMap
  * inside the query collection
  */
 template<class K>
 GenomicRegionCollection<GenomicRegion> FindOverlaps(const GenomicRegionCollection<K> &subject, std::vector<int32_t>& query_id, std::vector<int32_t>& subject_id, bool ignore_strand) const;

 /** Return the overlaps between the collection and the query interval
  * @param gr Query region 
  * @param ignore_strand If true, won't exclude overlap if on different strand
  * @return A collection of overlapping intervals from this collection, trimmed to be contained
  * inside gr
  */
 template<class K>
 GenomicRegionCollection<GenomicRegion> FindOverlaps(const K& gr, bool ignore_strand) const;

 /** Return the number of bases in query that overlap this collection 
  * @param gr Query GenomicRegion (or child of)
  * @param ignore_strand If true, won't exclude overlap if on different strand
  * @return Number of bases in query that overlap with region in collection
  */
 template<class K>
 size_t FindOverlapWidth(const K& gr, bool ignore_strand) const;

 /** Return the total amount spanned by this collection */
 int TotalWidth() const; 

 /** Increase the left and right ends of each contained GenomicRegion by 
  * the pad value.
  * @param v Amount to pad each end by. Result is increase in width by 2*pad.
  * @note See GenomicRegion::Pad
  */
 void Pad(int v);

 /** Set the i'th GenomicRegion */
 T& operator[](size_t i) { return m_grv->at(i); }
 
 /** Retreive the i'th GenomicRegion */
 const T& operator[](size_t i) const { return m_grv->at(i); }
 
  /** Add two GenomicRegionCollection objects together */
  void Concat(const GenomicRegionCollection<T>& g);

  /** Output the GenomicRegionCollection to a BED format
   *
   * @param h Header to convert id to chromosome name
   * @return BED formated string reprsentation 
   */
  std::string AsBEDString(const BamHeader& h) const;

 /** Coordinate sort the interval collection */
  void CoordinateSort();

 /** Expand all the elements so they are sorted and become adjacent 
  * by stretching them to the right up to max 
  * @param max Element furthest to the right will be stretched to max. If set to 0, will not stretch furthest right element.
  * @exception Throws an out_of_range if furthest right position is > max
  */
  void SortAndStretchRight(int max);

 /** Expand all the elements so they are sorted and become adjacent 
  * by stretching them to the left down to min.
  * @param min Element furthest to the left will be stretched to min. If set to < 0, will not stretch furthest right element.
  * @exception Throws an out_of_range if furthest left is < min
  */
 void SortAndStretchLeft(int min);

  /** Rewind the element pointer to the first GenomicRegion */
  void Rewind() { idx = 0; }

 /** Return elements as an STL vector of GenomicRegion objects */
 GenomicRegionVector AsGenomicRegionVector() const;
 
 /** Iterator to first element of the region collection */
 typename std::vector<T>::iterator begin() { return m_grv->begin(); } 

 /** Iterator to end of the region collection */ 
 typename std::vector<T>::iterator end() { return m_grv->end(); } 

 /** Const iterator to first element of the region collection */  
 typename std::vector<T>::const_iterator begin() const { return m_grv->begin(); } 
 
 /** Const iterator to end of the region collection */ 
 typename std::vector<T>::const_iterator end() const { return m_grv->end(); } 

  /** Shortcut to FindOverlaps that just returns the intersecting regions
   * without keeping track of the query / subject ids
   * @param subject Collection of regions to intersect this object with
   * @param ignore_strand Ignore strand considerations when performing intersection
   * @return Intersecting regions between subject and query
   */
  template <class K>
  GenomicRegionCollection<GenomicRegion> Intersection(const GenomicRegionCollection<K>& subject, bool ignore_strand) const;
 
 protected:

 bool m_sorted;
 
 // always construct this object any time m_grv is modifed
 SeqPointer<GenomicIntervalTreeMap> m_tree;
 
 // hold the genomic regions
 SeqPointer<std::vector<T> > m_grv; 
 
 // index for current GenomicRegion
 size_t idx;

 // open the memory
 void allocate_grc();

};

typedef GenomicRegionCollection<GenomicRegion> GRC;

}

#include "SeqLib/GenomicRegionCollection.cpp"

#endif

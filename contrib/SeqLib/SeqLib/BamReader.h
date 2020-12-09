#ifndef SEQLIB_BAM_POLYREADER_H
#define SEQLIB_BAM_POLYREADER_H

#include <cassert>
#include "SeqLib/ReadFilter.h"
#include "SeqLib/BamWalker.h"
#include "SeqLib/ThreadPool.h"

// forward declare this from hts.c
extern "C" {
int hts_useek(htsFile *file, long uoffset, int where);
}


namespace SeqLib {

  class BamReader;

  typedef SeqPointer<hts_idx_t> SharedIndex; ///< Shared pointer to the HTSlib index struct

  typedef SeqPointer<htsFile> SharedHTSFile; ///< Shared pointer to the HTSlib file pointer
 
  // store file accessors for single BAM
  class _Bam {

    friend class BamReader;

  public:

  _Bam(const std::string& m) : m_region_idx(0), m_in(m), empty(true), mark_for_closure(false)  {}

  _Bam() : m_region_idx(0), empty(true), mark_for_closure(false) {}

    //! Return the header for this BAM
    const BamHeader& GetHeader() const {
      return m_hdr;
    }
    
    ~_Bam() {}

    std::string GetFileName() const { return m_in; }

    // point index to this region of bam
    bool SetRegion(const GenomicRegion& gp);

    // which region are we on
    size_t m_region_idx;

  private:

    // do the read loading
    // the return value here is just passed along from sam_read1
    int32_t load_read(BamRecord& r);

    void set_pool(ThreadPool t) {
      if (t.IsOpen() && fp) // probably dont need this, it can handle null
	hts_set_opt(fp.get(),  HTS_OPT_THREAD_POOL, &t.p); //t.p is htsThreadPool
    }

    void reset() {
      empty = true;
      mark_for_closure = false;
      m_region_idx = 0;
    }

    // close this bam
    bool close() {
      if (!fp)
	return false;
      fp.reset();
      idx.reset();
      hts_itr.reset();

      //fp = nullptr; // calls destructor actually
      //idx = nullptr;
      //hts_itr = nullptr;

      empty = true;
      mark_for_closure = false;
      m_region_idx = 0;

      return true;
    }

    // set a pre-loaded index (save on loading each time)
    //void set_index(SharedIndex& i) { idx = i; }

    // set a pre-loaded htsfile (save on loading each time)
    //void set_file(SharedHTSFile& i) { fp = i; }
    
    // set a pre-loaded index and make a deep copy
    //void deep_set_index();

    GRC* m_region; // local copy of region

    SharedHTSFile fp;     // BAM file pointer
    SharedIndex idx;  // bam index
    SeqPointer<hts_itr_t> hts_itr; // iterator to index location
    std::string m_in;                   // file name
    BamHeader m_hdr;                    // the BAM header

    // the next read "slotted" for this BAM
    BamRecord next_read;

    // the next read "slot" is empty
    bool empty;
    
    // if set to true, then won't even attempt to lookup read
    bool mark_for_closure;
    
    // open the file pointer
    bool open_BAM_for_reading(SeqLib::ThreadPool t);

    // hold the reference for CRAM reading
    std::string m_cram_reference;

  };

  typedef SeqHashMap<std::string, _Bam> _BamMap;
  
/** Stream in reads from multiple BAM/SAM/CRAM or stdin */
class BamReader {

 public:

  /** Construct an empty BamReader */
  BamReader();

  /** Destroy a BamReader and close all connections to the BAMs 
   * 
   * Calling the destructor will take care of all of the C-style dealloc
   * calls required within HTSlib to close a BAM or SAM file. 
   */
  ~BamReader() { }

  /** Explicitly set a reference genome to be used to decode CRAM file.
   * If no reference is specified, will automatically load from
   * file pointed to in CRAM header using the SQ tags. 
   * @note This function is useful if the reference path pointed
   * to by the UR field of SQ is not on your system, and you would
   * like to explicitly provide one.
   * @param ref Path to an index reference genome
   */
  void SetCramReference(const std::string& ref);

  /** Set a part of the BAM to walk.
   *
   * This will set the BAM pointer to the given region.
   * @note This clears all other regions and resets the index
   * pointer to this location
   * @param gp Location to point the BAM to
   * @return true if the region is found in the index
   */
  bool SetRegion(const GenomicRegion& gp);

  /** Assign this BamReader a thread pool
   * 
   * The thread pool with stay with this object, but
   * will not be created or destroyed. This must be done
   * separately, which allows for multiple readers/writers
   * to be connected to one thread pool
   * @return false if the thread pool has not been opened 
   */
  bool SetThreadPool(ThreadPool p);
  
  /** Set up multiple regions. Overwrites current regions. 
   * 
   * This will set the BAM pointer to the first element of the
   * input list.
   * @note This clears all other regions and resets the index
   * pointer to the first element of grc
   * @param grc Set of location to point BAM to
   * @return true if the regions are found in the index
   */
  bool SetMultipleRegions(const GRC& grc);

  /** Return if the reader has opened the first file */
  bool IsOpen() const { if (m_bams.size()) return m_bams.begin()->second.fp.get() != NULL; return false; }

  /*
  Set pre-loaded raw htslib index 
    Provide the reader with an index structure that is already loaded.
    This is useful if there are multiple newly created BamReader objects
    that use the same index (e.g. make a BAM index in a loop)
    @note This does not make a copy, so ops on this index are shared with
    every other object that controls it.
    @param i Pointer to an HTSlib index
    @param f Name of the file to set index for
    @return True if the file f is controlled by this object
  bool SetPreloadedIndex(const std::string& f, SharedIndex& i);

  Return a shared pointer to the raw htsFile object
    @exception Throws runtime_error if the requested file has not been opened already with Open
    @param f File to retrieve the htsFile from.
   
  SharedHTSFile GetHTSFile (const std::string& f) const;

  Return a shared pointer to the raw htsFile object from the first BAM
    @exception Throws runtime_error if no files have been opened already with Open
    @param f File to retrieve the htsFile from.
  SharedHTSFile GetHTSFile () const;

  Set a pre-loaded raw index, to the first BAM
   @note see SetPreloadedIndex(const std::string& f, SharedIndex& i)
  bool SetPreloadedIndex(SharedIndex& i);
  */

  /** Return if the reader has opened the file
   * @param f Name of file to check
   */
  bool IsOpen(const std::string& f) const { 
    SeqHashMap<std::string, _Bam>::const_iterator ff = m_bams.find(f);
    if (ff == m_bams.end())
      return false;
    return ff->second.fp.get() != NULL; 
  }

  /** Close all of the BAMs */
  bool Close();

  /** Close a particular BAM/CRAM/SAM
   * @param f Particular file to close
   * @return True if BAM is found and is closable (eg no already closed)
   */
  bool Close(const std::string& f);

  /** Reset the given BAM/SAM/CRAM to the begining, but keep the loaded indicies and file-pointers 
   * @param f Name of file to reset
   * @return Returns false if this BAM is not found in object
   * @note Unlike Reset(), this version will NOT reset the regions, since other BAMs may still be
   * using them.
   */
  bool Reset(const std::string& f);

  /** Return a vector of all of the BAM/SAM/CRAMs in this reader */
  std::vector<std::string> ListFiles() const {
    std::vector<std::string> out;
    for (_BamMap::const_iterator i = m_bams.begin(); i != m_bams.end(); ++i)
      out.push_back(i->first);
    return out;
  }

  /** Create a string representation of 
   * all of the regions to walk
   */
  std::string PrintRegions() const;

  /** Print out some basic info about this reader */
  friend std::ostream& operator<<(std::ostream& out, const BamReader& b);

  /** Open a BAM/SAM/CRAM/STDIN file for streaming in 
   * @param bam Path to a SAM/CRAM/BAM file, or "-" for stdin
   * @return True if open was successful
   */
  bool Open(const std::string& bam);

  /** Open a set of BAM/SAM/CRAM/STDIN files for streaming in 
   * @param bams Path to a vector fo SAM/CRAM/BAM files, or "-" for stdin
   * @return True if open was successful
   */
  bool Open(const std::vector<std::string>& bams);

  /** Retrieve the next read from the available input streams.
   * @note Will chose the read with the lowest left-alignment position
   * from the available streams.
   * @param r Read to fill with data
   * @return true if the next read is available
   */
  bool GetNextRecord(BamRecord &r);

  /** Reset all the regions, but keep the loaded indicies and file-pointers */
  void Reset();

  /** Return a copy of the header to the first file 
   * @note The object returned is a copy of the BamHeader, but 
   * this does not actually copy the actual header contents. Header contents
   * are stored in a shared_ptr, and so the new returned BamHeader
   * have a copy of the shared_ptr that will point to the originally alloced 
   * raw header data.
   */
  BamHeader Header() const;

  /** Return a concatenation of all the headers */
  std::string HeaderConcat() const;

 protected:

  GRC m_region; ///< Regions to access

  _BamMap m_bams; ///< store the htslib file pointers etc to BAM files

 private:
  // hold the reference for CRAM reading
  std::string m_cram_reference;

  // for multicore reading/writing
  ThreadPool pool;

};


}
#endif 



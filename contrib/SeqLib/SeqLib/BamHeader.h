#ifndef SEQLIB_BAM_HEADER_H__
#define SEQLIB_BAM_HEADER_H__

#include "htslib/htslib/hts.h"
#include "htslib/htslib/sam.h"
#include "htslib/htslib/bgzf.h"
#include "htslib/htslib/kstring.h"

#include "SeqLib/SeqLibUtils.h"
#include <string>
#include <vector>

namespace SeqLib {

  /** Store a reference chromosome and its length
   * @note This parallels the data found in SQ tag of BAM header
   */
  struct HeaderSequence {
    
    /** Make a new header sequence 
     * @param n Name of the chromosome
     * @param l Length of the chromosome
     */
  HeaderSequence(const std::string& n, uint32_t l) : Name(n), Length(l) {}

    std::string Name; ///< Name of the sequence (eg "1")
    uint32_t Length; ///< Length of the sequence (eg LN:191469)
  };

  typedef std::vector<HeaderSequence> HeaderSequenceVector;

  /** Store a header to a BAM file 
   *
   * Stores a BAM header, which also acts as a dictionary of 
   * reference sequences, with names and lengths.
   */
  class BamHeader {

  public:

    /** Initializes a new empty BamHeader with no data
     * 
     * @note No memory is allocated here
     */
    BamHeader() {};
    
    /** Construct a new header from ref sequences and lengths 
     *
     */
    BamHeader(const HeaderSequenceVector& hsv);

    /** Initialize a BamHeader from a string containing
     * a BAM header in human-readable form (e.g. PG ... )
     * @param hdr Text of a BAM header, with newlines separating lines
     */
    BamHeader(const std::string& hdr);

    /** Create a new BamHeader from a raw HTSlib header.
     * 
     * @note This will make a copy of the input header
     */
    BamHeader(const bam_hdr_t * hdr);
    
    /** Return the number of sequences store in this dictionary
     * Returns 0 if header is unitialized.
     */
    int NumSequences() const;

    /** Return the length of the sequence */
    int GetSequenceLength(int id) const;

    /** Return the length of the sequence */
    int GetSequenceLength(const std::string& id) const;

    /** Return if the header has been opened  */
    bool IsOpen() const { return h.get() != NULL; }

    /** Return the full text of the header */
    std::string AsString() const;

    /** Convert a numeric sequence ID to a name
     * 
     * @exception Throws an out_of_range if ID is >= then number of 
     * targets in dictionary, or if header is uninitialized..
     * @exception Throws an invalid_argument if ID is < 0;
     */
    std::string IDtoName(int id) const;

    /** Check if the header has been initialized
     */
    bool isEmpty() const { return h.get() == NULL; }

    /** Return the raw bam_hdr_t */
    const bam_hdr_t* get() const { return h.get(); }

    /** Return the raw bam_hdr_t */
    bam_hdr_t* get_() const { return h.get(); }

    /** Get the numeric ID associated with a sequence name.
     * @param name Name of the sequence
     * @return ID of named sequence, or -1 if not in dictionary
     */
    int Name2ID(const std::string& name) const;

    /** Return the reference sequences as vector of HeaderSequence objects */
    HeaderSequenceVector GetHeaderSequenceVector() const;

  private:

    // adapted from sam.c - bam_nam2id
    int bam_name2id_2(const bam_hdr_t *h, const char *ref) const;

    SeqPointer<bam_hdr_t> h;

    // make the name 2 id map (to be used by Name2ID)
    // replaces part of bam_name2id that makes the hash table
    void ConstructName2IDTable();

    // hash table for name to id
    SeqPointer<SeqHashMap<std::string, int> > n2i;

    // adapted from sam_hdr_read
    bam_hdr_t* sam_hdr_read2(const std::string& hdr) const;

  };
  
}


#endif

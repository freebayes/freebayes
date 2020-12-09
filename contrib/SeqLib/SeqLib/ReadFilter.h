#ifndef SEQLIB_READ_FILTER_H
#define SEQLIB_READ_FILTER_H

#define AHO_CORASICK 1

#include <string>
#include <vector>
#include <climits>

#include "json/json.h"

#include "SeqLib/GenomicRegionCollection.h"
#include "SeqLib/BamRecord.h"

#ifdef HAVE_C11
#include "SeqLib/aho_corasick.hpp"
#endif


#define MINIRULES_MATE_LINKED 1
#define MINIRULES_MATE_LINKED_EXCLUDE 2
#define MINIRULES_REGION 3
#define MINIRULES_REGION_EXCLUDE 4

namespace SeqLib {

  typedef SeqHashSet<std::string> StringSet;

  namespace Filter {

#ifdef HAVE_C11
  /** Tool for using the Aho-Corasick method for substring queries of 
   * using large dictionaries 
   * @note Trie construction / searching implemented by https://github.com/blockchaindev/aho_corasick
   */
  struct AhoCorasick {
    
    /** Allocate a new empty trie */
    AhoCorasick() { 
      aho_trie = SeqPointer<aho_corasick::trie>(new aho_corasick::trie()); 
      inv = false;
      count = 0;
    } 

    /** Deallocate the trie */
    ~AhoCorasick() { }

    /** Add a motif to the trie 
     * @note Trie construction is lazy. Won't build trie until 
     * first query. Therefore first query is slow, the rest are
     * O(n) where (n) is length of query string.
     */
    void AddMotif(const std::string& m) { 
      aho_trie->insert(m);
    } 
    
    /** Add a set of motifs to the trie from a file 
     * @param f File storing the motifs (new line separated)
     * @exception Throws a runtime_error if file cannot be opened
     */
    void TrieFromFile(const std::string& f);

    /** Query if a string is in the trie 
     * @param t Text to query
     * @return Returns number of substrings in tree that are in t
     */
    int QueryText(const std::string& t) const;

    SeqPointer<aho_corasick::trie> aho_trie; ///< The trie for the Aho-Corasick search
    
    std::string file; ///< Name of the file holding the motifs

    bool inv; ///< Is this an inverted dictinary (ie exclude hits)
    
    int count; ///< Number of motifs in dictionary
    
  };
#endif

/** Stores a rule for a single alignment flag.
 *
 * Rules for alignment flags can be one of three states:
 * - NA - All flag values are valid
 * - Off - Flag is valid if OFF
 * - On - Flag is valid if ON
 */
class Flag {
  
 public:

  /** Construct a new Flag with NA rule 
   */
  Flag() : on(false), off(false), na(true) {}
    
  /** Set the flag to NA (pass alignment regardless of flag value) */
  void setNA() { on = false; off = false; na = true; } 
  
  /** Set the flag to ON (require flag ON to pass) */
  void setOn() { on = true; off = false; na = false; } 

  /** Set the flag to OFF (require flag OFF to pass) */
  void setOff() { on = false; off = true; na = false; } 

  /** Return if the Flag filter is NA */
  bool isNA()  const { return na; } 

  /** Return if the Flag filter is ON */
  bool isOn()  const { return on; } 

  /** Return if the Flag filter is OFF */
  bool isOff() const { return off; } 

  /** Parse the Flag rule from a JSON entry */
  bool parseJson(const Json::Value& value, const std::string& name);

 private: 

  bool on;
  bool off; 
  bool na;

};

/** Filter numeric values on whether they fall in/out of a range of values (eg mapping quality). 
 *
 * Can optionally invert the range to make rule the complement of the range
 * (eg insert-size NOT in [300,600]
 */
class Range {

public:
  /** Construct a default range with everything accepted
   */
  Range() : m_min(0), m_max(0), m_inverted(false), m_every(true) {}

  /** Construct a Range from min to max, inclusive
   * @param min Minimum for range 
   * @param max Maximum for range 
   * @param inverted Declare if this should be an inverted range (do NOT accept vals in range)
   */
  Range(int min, int max, bool inverted) : m_min(min), m_max(max), m_inverted(inverted), m_every(false) {}
  
  /** Given a query value, determine if the value passes this Range
   * @param val Query value (e.g. mapping quality)
   * @return true if the value passes this Range rule
   */
  bool isValid(int val) {
    if (m_every)
      return true;
    if (!m_inverted)
      return (val >= m_min && val <= m_max);
    else
      return (val < m_min || val > m_max);
  }

  /** Parse a JSON value 
   * @param value 
   * @param name
   */
  void parseJson(const Json::Value& value, const std::string& name);

  /** Print the contents of this Range */
  friend std::ostream& operator<<(std::ostream &out, const Range &r);

  /** Return if this range accepts all values */
  bool isEvery() const { return m_every; }
  
  /** Return the lower bound of the range */
  int lowerBound() const { return m_min; }

  /** Return the upper bound of the range */
  int upperBound() const { return m_max; }

  /** Return true if the range is inverted (e.g. do NOT accept i in [min,max] */
  bool isInverted() const { return m_inverted; }
  
private:
  
  int m_min;
  int m_max;
  bool m_inverted;
  bool m_every;

};

/** Stores a set of Flag objects for filtering alignment flags
 *
 * An alignment can be queried against a FlagRule to check if it 
 * satisfies the requirements for its alignment flag.
 */
class FlagRule {

 public:

  FlagRule() {
    dup  = Flag();
    supp       = Flag();
    qcfail     = Flag();
    hardclip   = Flag();
    fwd_strand = Flag();
    rev_strand = Flag();
    mate_fwd_strand = Flag();
    mate_rev_strand = Flag();
    mapped          = Flag();
    mate_mapped     = Flag();
    ff = Flag();
    fr = Flag();
    rf = Flag();
    rr = Flag();
    ic = Flag();
    paired = Flag();
    m_all_on_flag = 0;
    m_all_off_flag = 0;
    m_any_on_flag = 0;
    m_any_off_flag = 0;
    every = false;
  }
  
  Flag dup; ///< Filter for duplicated flag 
  Flag supp; ///< Flag for supplementary flag 
  Flag qcfail; ///< Flag for qcfail flag
  Flag hardclip; ///< Flag for presence of hardclip in cigar string
  Flag fwd_strand; ///< Flag for forward strand alignment
  Flag rev_strand; ///< Flag for reverse strand alignment
  Flag mate_fwd_strand; ///< Flag for forward strand alignment for mate 
  Flag mate_rev_strand; ///< Flag for reverse strand alignment for mate 
  Flag mapped; ///< Flag for mapped alignment
  Flag mate_mapped; ///< Flag for mate-mapped alignment
  Flag ff; ///< Flag for both reads on forward strand
  Flag fr; ///< Flag for lower (by position) read on forward strand, higher on reverse
  Flag rf; ///< Flag for lower (by position) read on reverse strand, higher on forward
  Flag rr; ///< Flag for both reads on reverse strand 
  Flag ic; ///< Flag for read and mate aligned to different chromosomes
  Flag paired; ///< Flag for read is part of pair

  /** Parse a FlagRule from a JSON object
   */
  void parseJson(const Json::Value& value);

  /** Set rule to pass all alignment flags that have any bit in f on
   * @param f Alignment flags to be on for alignment record to pass
   */  
  void setAnyOnFlag(uint32_t f) { m_any_on_flag = f;   every = (every && f == 0); } 
  //  NOTE: every = (every && f == 0) means to set every to true only if 
  //  input flag is zero and every was already true

  /** Set rule to pass all alignment flags that have any bit in f off
   * @param f Alignment flags to be off for alignment record to pass
   */  
  void setAnyOffFlag(uint32_t f) { m_any_off_flag = f; every = (every && f == 0); } 

  /** Set rule to reject all alignment flags that have any bit in f off
   * @param f Alignment flags to be on for alignment record to pass
   */
  void setAllOnFlag(uint32_t f) { m_all_on_flag = f;   every = (every && f == 0); } 

  /** Set rule to reject all alignment flags that have any bit in f on
   * @param f Alignment flags to be off for alignment record to pass
   */
  void setAllOffFlag(uint32_t f) { m_all_off_flag = f; every = (every && f == 0); } 

  /** Return whether a read passes the alignment flag rules in this object 
   * @param r Alignment record to query
   * @return true if record passes rules
   */
  bool isValid(const BamRecord &r);

  /** Print the flag rule */
  friend std::ostream& operator<<(std::ostream &out, const FlagRule &fr);

  /** Return if every this object will pass all records provided to it */
  bool isEvery() const { return every; }

private:

  bool every; // does this pass all flags? 
  
  uint32_t m_all_on_flag;  // if read has all of these, keep
  uint32_t m_all_off_flag; // if read has all of these, fail

  uint32_t m_any_on_flag; // if read has any of these, keep
  uint32_t m_any_off_flag;// if read has any of these, fail

  int parse_json_int(const Json::Value& v);

};

/** Stores a full rule (Flag + Range + motif etc)
 *
 * An alignment can be queried with an AbstractRule object
 * to check if it passes that rule.
 */
class AbstractRule {

  friend class ReadFilter;
  friend class ReadFilterCollection;

 public:

  /** Create empty rule with default to accept all */
 AbstractRule() : m_count(0), subsam_frac(1), subsam_seed(999) { }

  /** Destroy the filter */
  ~AbstractRule() {}

  /** Add a list of motifs that will be search as sub-strings
   * of the read sequence
   * @param f Path to new-line separted file of motifs
   * @param inverted If true, the reads that have a matching motif will fail isValid
   */
  void addMotifRule(const std::string& f, bool inverted);

  /** Query a read against this rule. If the
   * read passes this rule, return true.
   * @param r An aligned sequencing read to query against filter
   */
  bool isValid(const BamRecord &r);

  /** Supply the rule parameters with a JSON
   * @param value JSON object created by parsing a string
   */
  void parseJson(const Json::Value& value);

  /** Print some basic information about this filter
   */
  friend std::ostream& operator<<(std::ostream &out, const AbstractRule &fr);

  /** Return if this rule accepts all reads
   */
  bool isEvery() const;

  /** Set the rate to subsample (default 1 = no subsampling) 
   * @param s A rate between 0 and 1
   */
  void SetSubsampleRate(double s) { subsam_frac = s; };

  /** Supply a name for this rule 
   * @param s ID to be associated with this rule
   */
  void SetRuleID(const std::string& s) { id = s; };

  /** Specify a read-group for this filter.
   * Reads that do not belong to this read group
   * will not pass isValid
   * @param rg Read group to be matched against RG:Z:readgroup
   */
  void SetReadGroup(const std::string& rg) { read_group = rg; }

  FlagRule fr; ///< FlagRule specifying the alignment flag filter

  Range isize; ///< Range object for insert-size filter
  Range mapq; ///< Range object for mapping quality filter
  Range len; ///< Range object for length filter
  Range phred; ///< Range object for base-quality filter
  Range clip; ///< Range object for number of clipped bases filter
  Range nm; ///< Range object for NM (num mismatch) filter
  Range nbases; ///< Range object for number of "N" bases filer
  Range ins; ///< Range object for max CIGAR insertion size filter
  Range del; ///< Range object for max CIGAR deletion size filter
  Range xp; ///< Range object for number of secondary alignments

 private:

  void parseSeqLine(const Json::Value& value);
  
  // read group 
  std::string read_group;

  // how many reads pass this rule?
  size_t m_count;

  // the aho-corasick trie
#ifdef HAVE_C11
  AhoCorasick aho;
#endif

  // id for this rule
  std::string id;

  // fraction reads to subsample
  double subsam_frac;

  // data
  uint32_t subsam_seed; // random seed for subsampling

  void parseSubLine(const Json::Value& value);

};

class ReadFilterCollection;

/** 
 * A set of AbstractRules on a region united by logi rules.
 *
 * ReadFilter stores an arbitrarily complex collection of AbstractRules
 * (e.g. (Mapped && Clipped) || (Unmapped)).
 */
class ReadFilter {
  
  friend class ReadFilterCollection;

  public:

  /** Construct an empty filter that passes all reads */
  ReadFilter() : excluder(false), m_applies_to_mate(false), m_count(0) {}

  /** Destroy the filter */
  ~ReadFilter();

  // ReadFilter(const ReadFilter& rf);

  // Make a ReadFilter with an all exclude or include rule
  // @param Samtools style string, BED file or VCF 
  // @param reg_type The type of rule this will be
  // @param h BAM header that defines available chromosomes
  ///
  //ReadFilter(const CommandLineRegion& c, const BamHeader& hdr);

  /** Return whether a read passes this filter
   * @param r A read to query
   * @note If this is an excluder rule, then this
   * returns false if the read passes the filter
   */
  bool isValid(const BamRecord &r);

  /** Add a rule to this filter. A read must pass all 
   * of the rules contained in this filter to pass 
   * @param ar A rule (eg MAPQ > 30) that the read must satisfy to pass this filter.
   */
  void AddRule(const AbstractRule& ar);
   
  /** Provide the region covered by this read filter
   * @param g Region that this filter applies to
   */
  void setRegions(const GRC& g);

  /** Add additional regions to the filtered region 
   * @param g Additional regions to be included in filter
   */
  void addRegions(const GRC& g);
  
  /** Check if a read is overlapping the region defined by this filter 
   * @param r Read to query whether it overlaps (even partially) the region.
   * @note If this is a mate-linked region, then the read will overlap
   * if its mate overlaps as well.
   */
  bool isReadOverlappingRegion(const BamRecord &r) const;

  /** Print basic information about this filter */
  friend std::ostream& operator<<(std::ostream& out, const ReadFilter &mr);

  /** Return the number of rules in this filter */
  size_t size() const {
    return m_abstract_rules.size();
  }

  /** Set as an excluder region 
   * An excluder region is such that if a read satisfies
   * this rule, then it will fail isValid, rather than pass
   */
  void SetExcluder(bool e) { excluder = e; }

  /** Set as a mate linked region */
  void SetMateLinked(bool e) { m_applies_to_mate = e; }
  
 private:

  GRC m_grv; // the interval tree with the regions this rule applies to. Empty is whole-genome

  std::string id; // set a unique id for this filter
 
  bool excluder; // this filter is such that if read passes, it gets excluded
 
  std::string m_region_file;

  std::vector<AbstractRule> m_abstract_rules; // hold all of the rules

  // rule applies to mate too
  bool m_applies_to_mate;

  // how many reads pass this MiniRule
  size_t m_count;

};

/** A full set of rules across any number of regions
 *
 * Stores the entire set of ReadFilter, each defined on a unique interval.
 * A single ReadFilterCollection object is sufficient to store any combination of rules,
 * and is the highest in the rule hierarchy. (ReadFilterCollection stores ReadFilter 
 * stores AbstractRules stores FlagRule/Ranges).
 */
class ReadFilterCollection {

 public: 
  
  /** Construct an empty ReadFilterCollection 
   * that will pass all reads.
   */
 ReadFilterCollection() : m_count(0), m_count_seen(0) {}

  /** Create a new filter collection directly from a JSON 
   * @param script A JSON file or directly as JSON formatted string
   * @param h BamHeader to convert chr sequence to id
   * @exception invalid_argument if cannot parse JSON
   */
  ReadFilterCollection(const std::string& script, const SeqLib::BamHeader& h);

  /** Add a new rule to the collection. 
   * If a read passes this rule, it will be included,
   * even if it fails the other filters. Or, if this filter
   * has the excluder tag, then if a read passes this filter
   * then it will be excluded, regardless of the other filters
   */
  void AddReadFilter(const ReadFilter& rf);

  /** Provide a global rule set (applies to each filter)
   * @param rule A filter specified in JSON format
   */
  void addGlobalRule(const std::string& rule);

  /** Query a read to see if it passes any one of the
   * filters contained in this collection */
  bool isValid(const BamRecord &r);
  
  /** Print some basic information about this object */
  friend std::ostream& operator<<(std::ostream& out, const ReadFilterCollection &mr);

  /** Return a GenomicRegionCollection of all
   * of the regions specified by the filters.
   * @note This returns the raw regions. It may be useful
   * to run mergeOverlappingIntervals on the output to see
   * the minimal covered regions.
   */
  GRC getAllRegions() const;

  /** Return the number of filters in this collection */
  size_t size() const { return m_regions.size(); } 

  /** Return the total number of rules in this collection.
   * Filters are composed of collections of rules, and this
   * returns the total number of rules (e.g. MAPQ > 30) across
   * all of the filters
   */
  size_t numRules() const {
    size_t num = 0;
    for (std::vector<ReadFilter>::const_iterator it = m_regions.begin(); it != m_regions.end(); ++it)
      num += it->size();
    return num;
  }

  /** Check that the filter collection includes something. 
   * If not, then give it a universal includer
   */
  void CheckHasIncluder();

  // Return the a tab-delimited tally of which filters were satisfied.
   // Includes the header:
   // total_seen_count total_passed_count region region_passed_count rule rule_passed_count
   //
  //std::string EmitCounts() const;

 private:  

  // the global rule that all other rules are inherited from
  AbstractRule rule_all;

  size_t m_count; // passed
  size_t m_count_seen; // tested

  // store all of the individual filters
  std::vector<ReadFilter> m_regions;

  bool ParseFilterObject(const std::string& filterName, const Json::Value& filterObject);

};

}

}
#endif

#include "SeqLib/GenomicRegionCollection.h"

#include <iostream>
#include <sstream>
#include <fstream>
#include <cassert>
#include <set>
#include <stdexcept>
#include <algorithm>
#include <zlib.h>

#define GZBUFFER 65472

//#define DEBUG_OVERLAPS 1

namespace SeqLib {

  template<class T>
  GenomicRegionCollection<T>::GenomicRegionCollection(int width, int ovlp, const HeaderSequenceVector& h) {

    idx = 0;
    allocate_grc();

    // undefined otherwise
    if (width <= ovlp)
      throw std::invalid_argument("Width should be > ovlp");

    size_t chr = 0;
    for (HeaderSequenceVector::const_iterator i = h.begin(); i != h.end(); ++i) {
      
      T gr;
      gr.chr = chr;
      gr.pos1 = 0;
      gr.pos2 = i->Length;
      ++chr;

      if (width >= gr.Width()) {
	m_grv->push_back(gr);
	continue;
      }

      int32_t start = gr.pos1;
      int32_t end = gr.pos1 + width;
      
      // region is smaller than width
      if ( end > gr.pos2 ) {
	std::cerr << "GenomicRegionCollection constructor: GenomicRegion is smaller than bin width" << std::endl;
	return; 
      }
      
      // loop through the sizes until done
      while (end <= gr.pos2) {
	T tg;
	tg.chr = gr.chr;
	tg.pos1 = start;
	tg.pos2 = end;
	m_grv->push_back(tg);
	end += width - ovlp; // make the new one
	start += width - ovlp;
      }
      assert(m_grv->size() > 0);


    }
  }

  template<class T>
  void GenomicRegionCollection<T>::CoordinateSort() {
    
    if (m_grv) {
      std::sort(m_grv->begin(), m_grv->end());
      m_sorted = true;
    }
  }
  
  template<class T>
  void GenomicRegionCollection<T>::Shuffle() {
    std::random_shuffle ( m_grv->begin(), m_grv->end() );
  }

  template<class T>
  void GenomicRegionCollection<T>::SortAndStretchRight(int max) {

    if (!m_grv->size())
      return;
    
    CoordinateSort();

    if (max > 0 && max < m_grv->back().pos2) 
      throw std::out_of_range("GenomicRegionCollection::SortAndStrech Can't stretch to max, as we are already past max.");

    for (size_t i = 0; i < m_grv->size() - 1; ++i) 
      m_grv->at(i).pos2 = m_grv->at(i+1).pos1 - 1;

    if (max > 0)
      m_grv->back().pos2 = max;

  }

  template<class T>
  void GenomicRegionCollection<T>::SortAndStretchLeft(int min) {

    if (!m_grv->size())
      return;
    
    CoordinateSort();

    if (min >= 0 && min < m_grv->begin()->pos1) 
      throw std::out_of_range("GenomicRegionCollection::SortAndStrechLeft - Can't stretch to min, as we are already below min");

    if (min >= 0)
      m_grv->at(0).pos1 = min;

    for (size_t i = 1; i < m_grv->size(); ++i) 
      m_grv->at(i).pos1 = m_grv->at(i-1).pos2 + 1;

  }

template<class T>
bool GenomicRegionCollection<T>::ReadBED(const std::string & file, const BamHeader& hdr) {

  m_sorted = false;
  idx = 0;

  gzFile fp = NULL;
  fp = strcmp(file.c_str(), "-")? gzopen(file.c_str(), "r") : gzdopen(fileno(stdin), "r");

  if (file.empty() || !fp) {
    std::cerr << "BED file not readable: " << file << std::endl;
    return false;
  }

  // http://www.lemoda.net/c/gzfile-read/
  while (1) {

    int err;                    
    char buffer[GZBUFFER];
    gzgets(fp, buffer, GZBUFFER);
    int bytes_read = strlen(buffer);

    // get one line
    if (bytes_read < GZBUFFER - 1) {
      if (gzeof (fp)) break;
      else {
	const char * error_string;
	error_string = gzerror (fp, &err);
	if (err) {
	  fprintf (stderr, "Error: %s.\n", error_string);
	  exit (EXIT_FAILURE);
	}
      }
    }

    // prepare to loop through each field of BED line
    //size_t counter = 0;
    std::string chr, pos1, pos2;
    std::string line(buffer);
    std::istringstream iss_line(line);
    std::string val;
    if (line.find("#") != std::string::npos) 
      continue;
    
    // read first three BED columns
    iss_line >> chr >> pos1 >> pos2;

    // construct the GenomicRegion
    T gr(chr, pos1, pos2, hdr);
    
    if (gr.chr >= 0)
      m_grv->push_back(gr);
  }

  return true;
}

template<class T>
bool GenomicRegionCollection<T>::ReadVCF(const std::string & file, const BamHeader& hdr) {

  m_sorted = false;
  idx = 0;

  gzFile fp = NULL;
  fp = strcmp(file.c_str(), "-")? gzopen(file.c_str(), "r") : gzdopen(fileno(stdin), "r");
  
  if (file.empty() || !fp) {
    std::cerr << "VCF file not readable: " << file << std::endl;
    return false;
  }

  // http://www.lemoda.net/c/gzfile-read/
  while (1) {

    int err;                    
    char buffer[GZBUFFER];
    gzgets(fp, buffer, GZBUFFER);
    int bytes_read = strlen(buffer);

    // get one line
    if (bytes_read < GZBUFFER - 1) {
      if (gzeof (fp)) break;
      else {
	const char * error_string;
	error_string = gzerror (fp, &err);
	if (err) {
	  fprintf (stderr, "Error: %s.\n", error_string);
	  exit (EXIT_FAILURE);
	}
      }
    }

    // prepare to loop through each field of BED line
    std::string chr, pos;
    std::string line(buffer);
    std::istringstream iss_line(line);
    std::string val;
    if (line.empty() || line.at(0) == '#')
      continue;

    // read first two columnes
    iss_line >> chr >> pos;

    // construct the GenomicRegion
    T gr;
    try {
      gr = T(chr, pos, pos, hdr);
    } catch (...) {
      std::cerr << "...Could not parse pos: " << pos << std::endl << std::endl
		<< "...on line " << line << std::endl;
      
    }
    if (gr.chr >= 0) 
      m_grv->push_back(gr);
  }

  return true;
}

template<class T>
GenomicRegionCollection<T>::GenomicRegionCollection(const std::string &file, const BamHeader& hdr) {

  allocate_grc();

  idx = 0;

  // check if it's samtools-style file
  if (file.find(":") != std::string::npos) {
    m_sorted = true; // only one, so sorted
    m_grv->push_back(T(file, hdr));
    return;
  }

  // BED file
  if (file.find(".bed") != std::string::npos)
    ReadBED(file, hdr);
  // VCF file
  else if (file.find(".vcf") != std::string::npos) 
    ReadVCF(file, hdr);
  else // default is BED file
    ReadBED(file, hdr);    

}

// reduce a set of GenomicRegions into the minium overlapping set (same as GenomicRanges "reduce")
template <class T>
void GenomicRegionCollection<T>::MergeOverlappingIntervals() {

  // make the list
  std::list<T> intervals(m_grv->begin(), m_grv->end());

  intervals.sort();
  typename std::list<T>::iterator inext(intervals.begin());
  ++inext;
  for (typename std::list<T>::iterator i(intervals.begin()), iend(intervals.end()); inext != iend;) {
    if((i->pos2 >= inext->pos1) && (i->chr == inext->chr)) // change >= to > to not overlap touching intervals (eg [4,5][5,6])
      {
	if(i->pos2 >= inext->pos2) intervals.erase(inext++);
	else if(i->pos2 < inext->pos2)
	  { i->pos2 = inext->pos2; intervals.erase(inext++); }
      }
    else { ++i; ++inext; }
  }

  // move it over to a grv
  m_grv->clear(); // clear the old data 

  // c++11
  //std::vector<T> v{ std::make_move_iterator(std::begin(intervals)), 
  //    std::make_move_iterator(std::end(intervals)) };
  //m_grv->insert(m_grv->end(), v.begin(), v.end());

  // non c++11
  //std::vector<T> v;
  // v.push_back(std::make_move_iterator(std::begin(intervals)));
  //v.push_back(std::make_move_iterator(std::end(intervals)));
  //std::vector<T> v{ std::make_move_iterator(std::begin(intervals)), 
  //    std::make_move_iterator(std::end(intervals)) };
  //m_grv->insert(m_grv->end(), v.begin(), v.end());
  //m_grv->reserve(intervals.size());
  //m_grv->append(intervals.begin(), intervals.end());
  m_grv->insert(m_grv->end(), intervals.begin(), intervals.end());

  // clear the old interval tree
  m_tree->clear();
}

template <class T>
GenomicRegionVector GenomicRegionCollection<T>::AsGenomicRegionVector() const { 
  GenomicRegionVector gg;
  for (typename std::vector<T>::const_iterator i = m_grv->begin(); i != m_grv->end(); ++i)
    gg.push_back(GenomicRegion(i->chr, i->pos1, i->pos2, i->strand));
  return gg; 
} 
  
template <class T>
void GenomicRegionCollection<T>::CreateTreeMap() {

  if (!m_grv->size())
    return;

  // sort the genomic intervals
  if (!m_sorted)
    CoordinateSort();

  // loop through and make the intervals for each chromosome
  GenomicIntervalMap map;
  for (size_t i = 0; i < m_grv->size(); ++i) {
    map[m_grv->at(i).chr].push_back(GenomicInterval(m_grv->at(i).pos1, m_grv->at(i).pos2, i));
  }

  // for each chr, make the tree from the intervals
  //for (auto it : map) {
  for (GenomicIntervalMap::iterator it = map.begin(); it != map.end(); ++it) {
    GenomicIntervalTreeMap::iterator ff = m_tree->find(it->first);
    if (ff != m_tree->end())
      ff->second = GenomicIntervalTree(it->second);
    else
      m_tree->insert(std::pair<int, GenomicIntervalTree>(it->first, GenomicIntervalTree(it->second)));
    //old //m_tree[it.first] = GenomicIntervalTree(it.second);
  }

}

template<class T>
int GenomicRegionCollection<T>::TotalWidth() const { 
  int wid = 0; 
  for (typename std::vector<T>::const_iterator i = m_grv->begin(); i != m_grv->end(); ++i)
    //  for (auto& i : *m_grv) 
    wid += i->Width(); 
  return wid; 
}

// divide a region into pieces of width and overlaps
template<class T>
GenomicRegionCollection<T>::GenomicRegionCollection(int width, int ovlp, const T &gr) {

  idx = 0;
  allocate_grc();

  // undefined otherwise
  if (width <= ovlp)
    throw std::invalid_argument("Width should be > ovlp");
  if (width >= gr.Width()) {
    m_grv->push_back(gr);
    return;
  }

  int32_t start = gr.pos1;
  int32_t end = gr.pos1 + width;

  // region is smaller than width
  if ( end > gr.pos2 ) {
    std::cerr << "GenomicRegionCollection constructor: GenomicRegion is smaller than bin width" << std::endl;
    return; 
  }

  // loop through the sizes until done
  while (end <= gr.pos2) {
    m_grv->push_back(T(gr.chr, start, end));
    end += width - ovlp; // make the new one
    start += width - ovlp;
  }
  assert(m_grv->size() > 0);
  
  // finish the last one if we need to
  if (m_grv->back().pos2 != gr.pos2) {
    start = m_grv->back().pos2 - ovlp; //width;
    end = gr.pos2;
    m_grv->push_back(T(gr.chr, start, end));
  }

  m_sorted = true;

}

template<class T>
size_t GenomicRegionCollection<T>::CountOverlaps(const T &gr) const {

  if (m_tree->size() == 0 && m_grv->size() != 0) 
    {
      std::cerr << "!!!!!! WARNING: Trying to find overlaps on empty tree. Need to run this->createTreeMap() somewhere " << std::endl;
      return 0;
    }

  GenomicIntervalVector giv;

  GenomicIntervalTreeMap::const_iterator ff = m_tree->find(gr.chr);
  if (ff == m_tree->end())
    return 0;
  ff->second.findOverlapping(gr.pos1, gr.pos2, giv);
  return (giv.size());
}

  template<class T>
  template<class K>
  bool GenomicRegionCollection<T>::OverlapSameInterval(const K &gr1, const K &gr2) const {
    
    // events on diff chr do not overlap same bin
    if (gr1.chr != gr2.chr)
      return false;
    
    if (m_tree->size() == 0 && m_grv->size() != 0) {
      std::cerr << "!!!!!! WARNING: Trying to find overlaps on empty tree. Need to run this->createTreeMap() somewhere " << std::endl;
      return false;
    }
    
    GenomicIntervalTreeMap::const_iterator ff1 = m_tree->find(gr1.chr);
    GenomicIntervalTreeMap::const_iterator ff2 = m_tree->find(gr2.chr);
  if (ff1 == m_tree->end() || ff2 == m_tree->end())
    return false;

  // do the interval tree query
   GenomicIntervalVector giv1, giv2;
   ff1->second.findOverlapping(gr1.pos1, gr1.pos2, giv1);
   ff2->second.findOverlapping(gr2.pos1, gr2.pos2, giv2);

   if (!giv1.size() || !giv2.size())
     return false;

   // each one only overlapped one element
   if (giv1.size() == 1 && giv2.size() == 1)
     return (giv1[0].value == giv2[0].value);

  // make a set of the possible starts
  SeqHashSet<int> vals;
  //  for (auto& i : giv1)
  for (GenomicIntervalVector::iterator i = giv1.begin(); i != giv1.end(); ++i)
    vals.insert(i->value);
  
  // loop the other side and see if they mix
  for (GenomicIntervalVector::iterator i = giv2.begin(); i != giv2.end(); ++i)
    if (vals.count(i->value))
      return true;

  return false;
  
  }

template<class T>
std::string GenomicRegionCollection<T>::AsBEDString(const BamHeader& h) const {
  
  if (m_grv->size() ==  0)
    return std::string(); 

  std::stringstream ss;
  //for (auto& i : *m_grv)
  for (typename std::vector<T>::const_iterator i = m_grv->begin(); i != m_grv->end(); ++i)
    ss << i->ChrName(h) << "\t" << i->pos1 << "\t" << i->pos2 << "\t" << i->strand << std::endl;

  return ss.str();

}

template<class T>
void GenomicRegionCollection<T>::Concat(const GenomicRegionCollection<T>& g)
{
  if (!g.size())
    return;
  m_sorted = false;
  m_grv->insert(m_grv->end(), g.m_grv->begin(), g.m_grv->end());
}

template<class T>
GenomicRegionCollection<T>::GenomicRegionCollection() {
  idx = 0;
  allocate_grc();
}

template<class T>
GenomicRegionCollection<T>::~GenomicRegionCollection() {
}


template<class T>
void GenomicRegionCollection<T>::allocate_grc() {
  m_sorted = false;
  m_grv =  SeqPointer<std::vector<T> >(new std::vector<T>()) ;
  m_tree = SeqPointer<GenomicIntervalTreeMap>(new GenomicIntervalTreeMap()) ;
}

template<class T>
GenomicRegionCollection<T>::GenomicRegionCollection(const BamRecordVector& brv) {
  idx = 0;

  allocate_grc();

  //for (auto& i : brv) 
  for (BamRecordVector::const_iterator i = brv.begin(); i != brv.end(); ++i) 
    m_grv->push_back(GenomicRegion(i->ChrID(), i->Position(), i->PositionEnd()));

}

template<class T>
const T& GenomicRegionCollection<T>::at(size_t i) const
{ 
  if (i >= m_grv->size()) 
    throw 20;
  return m_grv->at(i); 
}  


// this is query
template<class T>
template<class K>
std::vector<int> GenomicRegionCollection<T>::FindOverlappedIntervals(const K& gr, bool ignore_strand) const {  

  if (m_tree->size() == 0 && m_grv->size() != 0) 
    throw std::logic_error("Need to run CreateTreeMap to make the interval tree before doing range queries");
  
  // which chr (if any) are common between query and subject
  GenomicIntervalTreeMap::const_iterator ff = m_tree->find(gr.chr);

  std::vector<int> output;  

  //must as least share a chromosome  
  if (ff == m_tree->end())
    return output;

  // get the subject hits
  GenomicIntervalVector giv;
  ff->second.findOverlapping(gr.pos1, gr.pos2, giv);

  for (GenomicIntervalVector::const_iterator i = giv.begin(); i != giv.end(); ++i)
    if (ignore_strand || m_grv->at(i->value).strand == gr.strand) 
      output.push_back(i->value);

  return output;

}

template<class T>
template<class K>
size_t GenomicRegionCollection<T>::FindOverlapWidth(const K& gr, bool ignore_strand) const {
  
  SeqLib::GRC out = FindOverlaps<K>(gr, ignore_strand);
  if (!out.size())
    return 0;

  // make sure merged down
  out.MergeOverlappingIntervals();
  
  size_t val = 0;
  for (size_t i = 0; i < out.size(); ++i)
    val += out[i].Width();

  return val;
}

// this is query
template<class T>
template<class K>
GenomicRegionCollection<GenomicRegion> GenomicRegionCollection<T>::FindOverlaps(const K& gr, bool ignore_strand) const
{  

  GenomicRegionCollection<GenomicRegion> output;

  if (m_tree->size() == 0 && m_grv->size() != 0) 
    throw std::logic_error("Need to run CreateTreeMap to make the interval tree before doing range queries");
  
  // which chr (if any) are common between query and subject
  GenomicIntervalTreeMap::const_iterator ff = m_tree->find(gr.chr);

  //must as least share a chromosome  
  if (ff == m_tree->end())
    return output;

  // get the subject hits
  GenomicIntervalVector giv;
  ff->second.findOverlapping(gr.pos1, gr.pos2, giv);
  
#ifdef DEBUG_OVERLAPS
  std::cerr << "ff->second.intervals.size() " << ff->second.intervals.size() << std::endl;
  for (auto& k : ff->second.intervals)
    std::cerr << " intervals " << k.start << " to " << k.stop << " value " << k.value << std::endl;
  std::cerr << "GIV NUMBER OF HITS " << giv.size() << " for query " << gr << std::endl;
#endif

  // loop through the hits and define the GenomicRegion
  for (GenomicIntervalVector::const_iterator j = giv.begin(); j != giv.end(); ++j) {
    //for (auto& j : giv) { // giv points to positions on subject
    if (ignore_strand || (m_grv->at(j->value).strand == gr.strand) ) {
#ifdef DEBUG_OVERLAPS
	std::cerr << "find overlaps hit " << j->start << " " << j->stop << " -- " << j->value << std::endl;
#endif
	output.add(GenomicRegion(gr.chr, std::max(static_cast<int32_t>(j->start), gr.pos1), std::min(static_cast<int32_t>(j->stop), gr.pos2)));
      }
  }

  return output;
  
}

  // this is query
  template<class T>
  template<class K>
GenomicRegionCollection<GenomicRegion> GenomicRegionCollection<T>::FindOverlaps(const GenomicRegionCollection<K>& subject, std::vector<int32_t>& query_id, std::vector<int32_t>& subject_id, bool ignore_strand) const
{  

  GenomicRegionCollection<GenomicRegion> output;
  if (subject.NumTree() == 0 && subject.size() != 0) {
    std::cerr << "!!!!!! findOverlaps: WARNING: Trying to find overlaps on empty tree. Need to run this->createTreeMap() somewhere " << std::endl;
    return output;
  }

  // we loop through query, so want it to be smaller
  if (subject.size() < m_grv->size() && m_grv->size() - subject.size() > 20) 
    std::cerr << "findOverlaps warning: Suggest switching query and subject for efficiency." << std::endl;

#ifdef DEBUG_OVERLAPS
  std::cerr << "OVERLAP SUBJECT: " << std::endl;
  for (auto& i : subject)
    std::cerr << i << std::endl;
#endif

  // loop through the query GRanges (this) and overlap with subject
  for (size_t i = 0; i < m_grv->size(); ++i) 
    {
      // which chr (if any) are common between query and subject
      GenomicIntervalTreeMap::const_iterator ff = subject.GetTree()->find(m_grv->at(i).chr);

      GenomicIntervalVector giv;

#ifdef DEBUG_OVERLAPS
      std::cerr << "TRYING OVERLAP ON QUERY " << m_grv->at(i) << std::endl;
#endif
      //must as least share a chromosome
      if (ff != m_tree->end())
	{
	  // get the subject hits
	  ff->second.findOverlapping(m_grv->at(i).pos1, m_grv->at(i).pos2, giv);

#ifdef DEBUG_OVERLAPS
	  std::cerr << "ff->second.intervals.size() " << ff->second.intervals.size() << std::endl;
	  for (auto& k : ff->second.intervals)
	    std::cerr << " intervals " << k.start << " to " << k.stop << " value " << k.value << std::endl;
	  std::cerr << "GIV NUMBER OF HITS " << giv.size() << " for query " << m_grv->at(i) << std::endl;
#endif
	  // loop through the hits and define the GenomicRegion
	  for (GenomicIntervalVector::const_iterator j = giv.begin(); j != giv.end(); ++j) {
	    //for (auto& j : giv) { // giv points to positions on subject
	    if (ignore_strand || (subject.at(j->value).strand == m_grv->at(i).strand) ) {
	      query_id.push_back(i);
	      subject_id.push_back(j->value);
#ifdef DEBUG_OVERLAPS
	      std::cerr << "find overlaps hit " << j->start << " " << j->stop << " -- " << j->value << std::endl;
#endif
	      output.add(GenomicRegion(m_grv->at(i).chr, std::max(static_cast<int32_t>(j->start), m_grv->at(i).pos1), std::min(static_cast<int32_t>(j->stop), m_grv->at(i).pos2)));
	    }
	  }
	}
    }

  return output;
  
}


template<class T>
GenomicRegionCollection<T>::GenomicRegionCollection(const T& gr)
{
  m_sorted = true;
  idx = 0;
  allocate_grc();
  m_grv->push_back(gr);
}

template<class T>
template<class K>
GRC GenomicRegionCollection<T>::Intersection(const GenomicRegionCollection<K>& subject, bool ignore_strand) const
{
  std::vector<int32_t> sub, que;
  GRC out;
  if (subject.size() > this->size()) // do most efficient ordering
    out = this->FindOverlaps<K>(subject, que, sub, ignore_strand);
  else
    out = subject.FindOverlaps(*this, que, sub, ignore_strand);    
  return out;
}

template<class T>
void GenomicRegionCollection<T>::Pad(int v)
{
  //for (auto& i : *m_grv)
  for (typename std::vector<T>::iterator i = m_grv->begin(); i != m_grv->end(); ++i)
    i->Pad(v);
}

}


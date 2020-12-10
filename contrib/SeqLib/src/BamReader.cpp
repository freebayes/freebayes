#include "SeqLib/BamReader.h"

//#define DEBUG_WALKER 1

namespace SeqLib {

// set the bam region
bool _Bam::SetRegion(const GenomicRegion& gp) {

  // mark it "open" again, may be new reads here
  mark_for_closure = false;
    
  //HTS set region 
  if ( (fp->format.format == 4 || fp->format.format == 6) && !idx)  // BAM (4) or CRAM (6)
    idx = SharedIndex(sam_index_load(fp.get(), m_in.c_str()), idx_delete());
  
  if (!idx) {
    if (m_in != "-")
      std::cerr << "Failed to load index for " << m_in << ". Rebuild samtools index" << std::endl;
    else
      std::cerr << "Random access with SetRegion not available for STDIN reading (no index file)" << std::endl;
    return false;
  }
  
  if (gp.chr >= m_hdr.NumSequences()) {
    std::cerr << "Failed to set region on " << gp << ". Chr ID is bigger than n_targets=" << m_hdr.NumSequences() << std::endl;
    return false;
  }
  
  // should work for BAM or CRAM
  hts_itr = SeqPointer<hts_itr_t>(sam_itr_queryi(idx.get(), gp.chr, gp.pos1, gp.pos2), hts_itr_delete());
  
  if (!hts_itr) {
    std::cerr << "Error: Failed to set region: " << gp << std::endl; 
    return false;
  }
  
  return true;
}

void BamReader::Reset() {
  for (_BamMap::iterator b = m_bams.begin(); b != m_bams.end(); ++b) 
     b->second.reset();
  m_region = GRC();
}

  bool BamReader::Reset(const std::string& f) {
    
    // cant reset what we don't have
    if (!m_bams.count(f))
      return false;
    m_bams[f].reset();
    return true;
}

  bool BamReader::Close() {
    
    bool success = true;
  for (_BamMap::iterator b = m_bams.begin(); b != m_bams.end(); ++b) 
      success = success && b->second.close();
    return success;
  }

  bool BamReader::Close(const std::string& f) {
    
    // cant close what we don't have
    if (!m_bams.count(f)) 
      return false;

    return m_bams[f].close();
  }

  /*  SharedHTSFile BamReader::GetHTSFile () const {
    if (!m_bams.size())
      throw std::runtime_error("No BAMs have been opened yet");
    return m_bams.begin()->second.fp;
  }

  SharedHTSFile BamReader::GetHTSFile(const std::string& f) const {
    _BamMap::const_iterator ff = m_bams.find(f);
    if (ff == m_bams.end())
      throw std::runtime_error("File " + f + " has not been opened yet");
    return ff->second.fp;
  }
  

  bool BamReader::SetPreloadedIndex(SharedIndex& i) {
    if (!m_bams.size())
      return false;
    m_bams.begin()->second.set_index(i);
    return true;
  }

  bool BamReader::SetPreloadedIndex(const std::string& f, SharedIndex& i) {
    if (!m_bams.count(f))
      return false;
    m_bams[f].set_index(i);
    return true;
  }

  */

  bool BamReader::SetRegion(const GenomicRegion& g) {
    m_region.clear();
    m_region.add(g);
    
    bool success = true;
    if (m_region.size()) {
      for (_BamMap::iterator b = m_bams.begin(); b != m_bams.end(); ++b) {
	b->second.m_region = &m_region;
	b->second.m_region_idx = 0; // set to the begining
	success = success && b->second.SetRegion(m_region[0]);
    }
    return success;
  }

  return false;
  
}

  bool BamReader::SetMultipleRegions(const GRC& grc) 
{
  if (grc.size() == 0) {
    std::cerr << "Warning: Trying to set an empty bam region"  << std::endl;
    return false;
  }
  
  m_region = grc;

  // go through and start all the BAMs at the first region
  bool success = true;
  if (m_region.size()) {
    for (_BamMap::iterator b = m_bams.begin(); b != m_bams.end(); ++b) {
      b->second.m_region = &m_region;
      b->second.m_region_idx = 0; // set to the begining
      success = success && b->second.SetRegion(m_region[0]);
    }
    return success;
  }
  
  return false;
}

  bool BamReader::Open(const std::string& bam) {

    // dont open same bam twice
    if (m_bams.count(bam))
      return false;
    
    _Bam new_bam(bam);
    if (!m_cram_reference.empty()) new_bam.m_cram_reference = m_cram_reference;
    new_bam.m_region = &m_region;
    bool success = new_bam.open_BAM_for_reading(pool);
    m_bams.insert(std::pair<std::string, _Bam>(bam, new_bam));
    return success;
  }

  bool BamReader::Open(const std::vector<std::string>& bams) {
    
    bool pass = true;
    for (std::vector<std::string>::const_iterator i = bams.begin(); i != bams.end(); ++i)
      pass = pass && Open(*i);
    return pass;
  }
  
BamReader::BamReader() {}

  std::string BamReader::HeaderConcat() const {
    std::stringstream ss;
    for (_BamMap::const_iterator i = m_bams.begin(); i != m_bams.end(); ++i) 
      ss << i->second.m_hdr.AsString();
    return ss.str();

  }

  bool BamReader::SetThreadPool(ThreadPool p) {
    if (!p.IsOpen())
      return false;
    pool = p;
    for (_BamMap::iterator b = m_bams.begin(); b != m_bams.end(); ++b)
      b->second.set_pool(p);
    return true;
  }
  
  BamHeader BamReader::Header() const { 
    if (m_bams.size()) 
      return m_bams.begin()->second.m_hdr; 
    return BamHeader(); 
  }

  bool _Bam::open_BAM_for_reading(SeqLib::ThreadPool t) {

    // HTS open the reader
    fp = SharedHTSFile(hts_open(m_in.c_str(), "r"), htsFile_delete()); 

    // connect the thread pool (may already be done, but its ok
    set_pool(t);

    // open cram reference
    if (!m_cram_reference.empty()) {
      int ret = hts_set_fai_filename(fp.get(), m_cram_reference.c_str());
      if (ret < 0) 
	throw std::invalid_argument("Could not read reference genome " + m_cram_reference + " for CRAM opt");
    }
   
    // check if opening failed
    if (!fp) 
      return false; 
    
    // read the header and create a BamHeader
    bam_hdr_t * hdr = sam_hdr_read(fp.get());
    m_hdr = BamHeader(hdr); // calls BamHeader(bam_hdr_t), makes a copy

    // deallocate the memory we just made
    if (hdr)
      bam_hdr_destroy(hdr);
    
    // if BAM header opening failed, return false
    if (!m_hdr.get()) 
      return false;
    
    // everything worked
    return true;
    
  }

  void BamReader::SetCramReference(const std::string& ref) {
    m_cram_reference = ref;
    for (_BamMap::iterator b = m_bams.begin(); b != m_bams.end(); ++b)
      b->second.m_cram_reference = ref;
  }

bool BamReader::GetNextRecord(BamRecord& r) {

  // shortcut if we have only a single bam
  if (m_bams.size() == 1) {
    
    if (m_bams.begin()->second.fp.get() == NULL || m_bams.begin()->second.mark_for_closure) // cant read if not opened
      return false;
    
    // try and get the next read
    int32_t status = m_bams.begin()->second.load_read(r);
    if (status >= 0)
      return true;
    if (status == -1) {
      // didn't find anything, clear it
      m_bams.begin()->second.mark_for_closure = true;
      return false;
    }
    
    // run time error
    std::stringstream ss;
    ss << "sam_read1 return status: " << status << " file: " << m_bams.begin()->first;
    throw std::runtime_error(ss.str());
    return false;
  }

  // loop the files and load the next read
  // for the one that was emptied last
  for (_BamMap::iterator bam = m_bams.begin(); bam != m_bams.end(); ++bam) {

    _Bam *tb = &(bam->second);

    // if marked, then don't even try on this BAM
    if (tb->mark_for_closure)
      continue;

    // skip un-opened BAMs
    if (tb->fp.get() == NULL) 
      continue;

    // if next read is not marked as empty, skip to next
    if (!tb->empty)
      continue; 
    
    // load the next read
    int32_t status = tb->load_read(r);
    if (status == -1) {
      // can't load, so mark for closing
      tb->empty = true;
      tb->mark_for_closure = true; // no more reads in this BAM
      continue; 
    } else if (status < 0) { // error sent back from sam_read1
      // run time error
      std::stringstream ss;
      ss << "sam_read1 return status: " << status << " file: " << bam->first;
      throw std::runtime_error(ss.str());
    }
    
  }

  // for multiple bams, choose the one to return
  // sort based on chr and left-most alignment pos. Same as samtools
  int min_chr = INT_MAX;
  int min_pos = INT_MAX;
  _BamMap::iterator hit; 
  bool found = false; // did we find a valid read

  for (_BamMap::iterator bam = m_bams.begin(); 
       bam != m_bams.end(); ++bam) {

    // dont check if already marked for removal or doesnt need new read
    if (bam->second.empty || bam->second.mark_for_closure) 
      continue;

    found = true;
    if (bam->second.next_read.ChrID() < min_chr || // check if read in this BAM is lowest
	(bam->second.next_read.Position() <  min_pos && bam->second.next_read.ChrID() == min_chr)) {
      min_pos = bam->second.next_read.Position();
      min_chr = bam->second.next_read.ChrID();
      hit = bam; // read is lowest, so mark this BAM as having the hit
    }
  }

  // mark the one we just found as empty
  if (found) {
    r = hit->second.next_read; // read is lowest, so assign
    hit->second.empty = true;  // mark as empty, so we fill this slot again
  }
  
  return found;
}
  
std::string BamReader::PrintRegions() const {

  std::stringstream ss;
  //for (GRC::const_iterator r = m_region.begin(); r != m_region.end(); ++r)
  //  ss << *r << std::endl;
  return(ss.str());

}

  int32_t _Bam::load_read(BamRecord& r) {

  // allocated the memory
  bam1_t* b = bam_init1(); 
  int32_t valid = -1; // start with EOF return code

  if (hts_itr.get() == NULL) {
    valid = sam_read1(fp.get(), m_hdr.get_(), b);    

    if (valid < 0) { 
      
#ifdef DEBUG_WALKER
      std::cerr << "ended reading on null hts_itr" << std::endl;
#endif
      //goto endloop;
      bam_destroy1(b);
      return valid;
    }
  } else {
    
    //changed to sam from hts_itr_next
    // move to next region of bam
    valid = sam_itr_next(fp.get(), hts_itr.get(), b);
  }
  
  if (valid < 0) { // read still not found
    do {
      
#ifdef DEBUG_WALKER
      std::cerr << "Failed read, trying next region. Moving counter to " << m_region_idx << " of " << m_region.size() << " FP: "  << fp_htsfile << " hts_itr " << std::endl;
#endif
      // try next region, return if no others to try
      ++m_region_idx; // increment to next region
      if (m_region_idx >= m_region->size()) {
	bam_destroy1(b);
	return valid;
      }
	//goto endloop;
      
      // next region exists, try it
      SetRegion(m_region->at(m_region_idx));
      valid = sam_itr_next(fp.get(), hts_itr.get(), b);
    } while (valid <= 0); // keep trying regions until works
  }
  
  // if we got here, then we found a read in this BAM
  empty = false;
  next_read.assign(b); // assign the shared_ptr for the bam1_t
  r = next_read;

  return valid;
}

std::ostream& operator<<(std::ostream& out, const BamReader& b)
{
  for(_BamMap::const_iterator bam = b.m_bams.begin(); bam != b.m_bams.end(); ++bam)
    out << ":" << bam->second.GetFileName() << std::endl; 

  if (b.m_region.size() && b.m_region.size() < 20) {
    out << " ------- BamReader Regions ----------" << std::endl;;
    //for (GRC::const_iterator r = b.m_region.begin(); r != b.m_region.end(); ++r)
    //  out << *i << std::endl;
  } 
  else if (b.m_region.size() >= 20) {
    int wid = 0;
    //for (GRC::const_iterator r = b.m_region.begin(); r != b.m_region.end(); ++r)
    //  wid += r->Width();
    out << " ------- BamReader Regions ----------" << std::endl;;
    out << " -- " << b.m_region.size() << " regions covering " << AddCommas(wid) << " bp of sequence"  << std::endl;
  }
  else 
    out << " - BamReader - Walking whole genome -" << std::endl;

  out <<   " ------------------------------------";
  return out;
}

}

#include "SeqLib/BamWalker.h"
#include "SeqLib/BamWriter.h"

#include <stdexcept>

//#define DEBUG_WALKER 1

namespace SeqLib {

  void BamWriter::SetHeader(const SeqLib::BamHeader& h) {
    hdr = h;
  }

  bool BamWriter::WriteHeader() const {
    
    if (hdr.isEmpty()) {
      std::cerr << "BamWriter::WriteHeader - No header supplied. Provide with SetWriteHeader" << std::endl;
      return false;
    }

    if (!fop) {
      std::cerr << "BamWriter::WriteHeader - Output not open for writing. Open with Open()" << std::endl;
      return false;
    }
    
    if (sam_hdr_write(fop.get(), hdr.get()) < 0) {
      std::cerr << "Cannot write header. sam_hdr_write exited with < 0" << std::endl;
      return false;
    }

    return true;
    
  }
  
  bool BamWriter::Close() {

    if (!fop)
      return false;

    fop.reset(); //tr1 compatible
    //fop = NULL; // this clears shared_ptr, calls sam_close (c++11)

    return true;
  }

bool BamWriter::BuildIndex() const {
  
  // throw an error if BAM is not already closed
  if (fop) {
    std::cerr << "Trying to index open BAM. Close first with Close()" << std::endl;
    return false;
  }

  if (m_out.empty()) {
    std::cerr << "Trying to make index, but no BAM specified" << std::endl;
    return false;
  }
  
  // call to htslib to build bai index
  if (sam_index_build(m_out.c_str(), 0) < 0) { // 0 is "min_shift", which is 0 for bai index
    std::cerr << "Failed to create index";
    return false;
  }

  return true;

}

  bool BamWriter::Open(const std::string& f) {

    // don't reopen
    if (fop)
      return false;

    m_out = f;

    // hts open the writer
    fop = SeqPointer<htsFile>(hts_open(m_out.c_str(), output_format.c_str()), htsFile_delete());

    // open the thread pool. It's OK if already connected before opening
    SetThreadPool(pool);

    if (!fop) {
      return false;
      //throw std::runtime_error("BamWriter::Open - Cannot open output file: " + f);
    }

    return true;
  }

  BamWriter::BamWriter(int o) {

    switch(o) {
    case BAM :  output_format = "wb"; break;
    case CRAM : output_format = "wc"; break;
    case SAM :  output_format = "w"; break;
    default : throw std::invalid_argument("Invalid writer type");
    }

  }
  

bool BamWriter::WriteRecord(const BamRecord &r)
{
  if (!fop) {
    return false;
  } else {
    if (sam_write1(fop.get(), hdr.get(), r.raw()) < 0)
      return false;
  }

  return true;
}

std::ostream& operator<<(std::ostream& out, const BamWriter& b)
{
  if (b.fop)
    out << "Write format: " << b.fop->format.format;
  out << " Write file " << b.m_out; 
  return out;
}

bool BamWriter::SetThreadPool(ThreadPool p) {
  if (!p.IsOpen()) 
    return false;
  pool = p;
  if (fop.get())
    hts_set_opt(fop.get(),  HTS_OPT_THREAD_POOL, &pool.p);
  return true;
}
  //does not return false if file not found
bool BamWriter::SetCramReference(const std::string& ref) {

  if (!fop)
    return false;

  // need to open reference for CRAM writing 
  char* fn_list = samfaipath(ref.c_str()); // eg ref = my.fa  returns my.fa.fai
  if (fn_list) {

    // in theory hts_set_fai_filename should give back < 0
    // if fn_list not there, but it doesnt
    if (!read_access_test(std::string(fn_list)))
      return false;
	
    int status = hts_set_fai_filename(fop.get(), fn_list); 
    if (status < 0) {
      fprintf(stderr, "Failed to use reference \"%s\".\n", fn_list);
      return false;
    }
  } else {
    std::cerr << "Failed to get the reference for CRAM compression" << std::endl;
    return false;
  }

  return true;
}

}

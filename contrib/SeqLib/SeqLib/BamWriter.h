#ifndef SEQLIB_BAM_WRITER_H
#define SEQLIB_BAM_WRITER_H

#include <cassert>
#include "SeqLib/BamRecord.h"
#include "SeqLib/ThreadPool.h"

namespace SeqLib {

  const int BAM = 4;
  const int SAM = 3;
  const int CRAM = 6;

/** Walk along a BAM or along BAM regions and stream in/out reads
 */
class BamWriter  {

 public:

  /** Construct an empty BamWriter to write BAM */
 BamWriter() : output_format("wb") {}

  /** Construct an empty BamWriter and specify output format 
   * @param o One of SeqLib::BAM, SeqLib::CRAM, SeqLib::SAM
   * @exception Throws an invalid_argument if not one of accepted values
   */
  BamWriter(int o);

  /** Destroy a BamWriter and close all connections to the BAM 
   * 
   * Calling the destructor will take care of all of the C-style dealloc
   * calls required within HTSlib to close a BAM or SAM file. 
   */
  ~BamWriter() {}

  /** Write the BAM header 
   * @return False if cannot write header
   */
  bool WriteHeader() const;

  /** Assign this BamWriter a thread pool
   * 
   * The thread pool with stay with this object, but
   * will not be created or destroyed. This must be done
   * separately, which allows for multiple readers/writers
   * to be connected to one thread pool
   * @return false if the thread pool has not been opened
   */
  bool SetThreadPool(ThreadPool p);

  /** Provide a header to this writer 
   * @param h Header for this writer. Copies contents
   */
  void SetHeader(const SeqLib::BamHeader& h);

  /** Close a file explitily. This is required before indexing with makeIndex.
   * @note If not called, BAM will close properly on object destruction
   * @return False if BAM already closed or was never opened
   */
  bool Close();

  /** Create the index file for the output bam in BAI format.
   *
   * This will make a call to HTSlib bam_index_build for the output file. 
   * @return Returns false if sam_index_build exits with < 0 status
   */
  bool BuildIndex() const;

  /** Print out some basic info about this writer */
  friend std::ostream& operator<<(std::ostream& out, const BamWriter& b);

  /** Open a BAM file for streaming out.
   * @param f Path to the output BAM/SAM/CRAM or "-" for stdout
   * @return False if cannot openf for writing
   */
  bool Open(const std::string& f);
  
  /** Return if the writer has opened the file */
  bool IsOpen() const { return fop.get() != NULL; }

  /** Write an alignment to the output BAM file 
   * @param r The BamRecord to save
   * @return False if cannot write alignment
   * @exception Throws a runtime_error if cannot write alignment
   */
  bool WriteRecord(const BamRecord &r);

  /** Explicitly set a reference genome to be used to decode CRAM file.
   * If no reference is specified, will automatically load from
   * file pointed to in CRAM header using the SQ tags. 
   * @note This function is useful if the reference path pointed
   * to by the UR field of SQ is not on your system, and you would
   * like to explicitly provide one.
   * @param ref Path to an index reference genome
   * @return Returns true if reference loaded.
   */
  bool SetCramReference(const std::string& ref);

  /** Return the BAM header */
  BamHeader Header() const { return hdr; };

 private:

  // path to output file
  std::string m_out; 

  // output format
  std::string output_format; 
  
  // hts
  SeqPointer<htsFile> fop;

  // header
  SeqLib::BamHeader hdr;

  // for multicore reading/writing
  ThreadPool pool;
  
};


}
#endif 



#ifndef SEQLIB_THREAD_POOL_H
#define SEQLIB_THREAD_POOL_H

#include <stdexcept>
#include "SeqLib/BamWalker.h"
#include "htslib/thread_pool.h"

namespace SeqLib{

class ThreadPool {

 public: 
  
 ThreadPool() : nthreads(1) { p.pool = NULL; }

 ThreadPool(int n) : nthreads(1) {
    p.pool = NULL;
    if (n < 1)
      throw std::invalid_argument( "n threads must be > 0");
    if (!(p.pool = hts_tpool_init(n))) 
      throw std::runtime_error( "Error creating thread pool");
  }

  bool IsOpen() { return p.pool != NULL; }

  htsThreadPool p;
  size_t nthreads;

};

}
#endif

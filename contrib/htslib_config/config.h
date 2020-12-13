// freebayes htslib configuration
//
// This file is used when building a local htslib. The local version has
// LZMA CRAM support. 

// #define HAVE_LIBBZ2 1   // may add BZ2 support for CRAM later
#define HAVE_LIBLZMA 1
#ifndef __APPLE__
#define HAVE_LZMA_H 1
#endif
#define HAVE_DRAND48 1
// #define HAVE_LIBCURL 1  


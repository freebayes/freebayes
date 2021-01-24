#ifndef FREEBAYES_LARGEFILESUPPORT_H
#define FREEBAYES_LARGEFILESUPPORT_H

#ifdef WIN32
#define ftell64(a)     _ftelli64(a)
#define fseek64(a,b,c) _fseeki64(a,b,c)
typedef __int64 off_type;
#elif defined(__CYGWIN__)
#define ftell64(a)     ftell(a)
#define fseek64(a,b,c) fseek(a,b,c)
typedef __int64_t off_type;
#elif defined(__APPLE__)
#define ftell64(a)     ftello(a)
#define fseek64(a,b,c) fseeko(a,b,c)
typedef off_t off_type;
#else
#define ftell64(a)     ftello(a)
#define fseek64(a,b,c) fseeko(a,b,c)
typedef __off64_t off_type;
#endif

#endif

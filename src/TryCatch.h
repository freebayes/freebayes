#ifndef FREEBAYES_TRYCATCH_H
#define FREEBAYES_TRYCATCH_H

#include <stdexcept> // out_of_range exception
#include <stdlib.h> // abort

// macros which improve our error handling
#ifndef TRY
#define TRY try
#endif
#ifndef CATCH
#define CATCH \
    catch (const std::out_of_range& outOfRange) { \
        cerr << "exception: " << outOfRange.what() \
        << " at line " << __LINE__ \
        << " in file " << __FILE__ << endl; \
        abort(); \
    }
#endif

#endif

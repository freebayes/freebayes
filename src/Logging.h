#ifndef FREEBAYES_LOGGING_H
#define FREEBAYES_LOGGING_H

// helper debugging macros to improve code readability
#define DEBUG(msg) \
    if (parameters.debug) { cerr << msg << endl; }

// lower-priority messages, enabled with "make debug"
#ifdef VERBOSE_DEBUG
#define DEBUG2(msg) \
    if (parameters.debug2) { cerr << msg << endl; }
#else
#define DEBUG2(msg)
#endif

// must-see error messages
#define ERROR(msg) \
    cerr << "ERROR(freebayes): " << msg << endl;

// must-see warning messages
#define WARNING(msg) \
    cerr << "WARNING(freebayes): " << msg << endl;

#endif

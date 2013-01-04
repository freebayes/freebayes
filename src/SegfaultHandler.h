#ifndef SEGFAULTHANDLER_H
#define SEGFAULTHANDLER_H

#include <stdio.h>
#include <execinfo.h>
#include <signal.h>
#include <stdlib.h>

void segfaultHandler(int sig);

#endif

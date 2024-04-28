# Contrib

This directory contains (stripped down) editions of projects that are not really actively maintained.

The goal is to have a clean repository of files that are in (active) use for freebayes.

## fastahack

libfastahack-dev is in Debian. This is actually a vcflib dependency that can be removed when vcflib moves fastahack out of its libvcflib.so (Variant.h)

## SmithWaterman

libsmithwaterman-dev is in Debian. This is actually a vcflib dependency that can be removed when vcflib moves SW out of its libvcflib.so (Variant.h)

## vcflib

libvcflib-dev is in Debian.

## tabixpp

libtabixpp-dev is in Debian.

## htslib

htslib is in Debian. To do a local build we have some prepared include files in contrib/htslib

## SeqLib

SeqLib sources are included for Guix. Should upstream that.

A subset of SeqLib was included that was patched for ARM64 using simde
in ssw.c:

    #define SIMDE_ENABLE_NATIVE_ALIASES
    #include <simde/x86/sse2.h>

I think the default simde is fine now and we no longer uses this submodule.

## ttmath

A few big double operations are included. The code gives warnings, but the newer version of ttmath does not solve that. We should look at using another library instead.

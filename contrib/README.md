# Contrib

This directory contains (stripped down) editions of projects that are
no longer actively maintained.

The goal is to have a clean repository of files that are in (active)
use for freebayes.

## Htslib

Freebayes has a submodule for htslib CRAM/BAM reading

## SeqLib

A subset of SeqLib is included that is patched for ARM64 using simde
in ssw.c:

    #define SIMDE_ENABLE_NATIVE_ALIASES
    #include <simde/x86/sse2.h>

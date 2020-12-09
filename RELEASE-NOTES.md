## ChangeLog v1.3.3 (2020????)

This is a maintenance release of Freebayes:

+ Freebayes compiles for ARM64
+ Freebayes supports CRAM and we added a test for CRAM and a hint in
  the README
+ Added Meson+Ninja build system for faster builds
+ Moved unused files from ./src to ./contrib/freebayes
+ Moved stripped SeqLib source into ./contrib/SeqLib as SeqLib is no
  longer maintained
+ Made htslib a primary dependency of freebayes
+ Updated htslib to the latest version

Some minor stuff:

+ Added regression test so we track some changes through git repo

## Earlier release notes

See github [releases](https://github.com/ekg/freebayes/releases)

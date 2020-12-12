## ChangeLog v1.3.3 (2020????)

This is a maintenance release of Freebayes:

+ Freebayes compiles for ARM64
+ Freebayes supports CRAM and we added a test for CRAM and a hint in
  the README
+ Added Meson+Ninja build system for faster builds
+ Merged: Set CRAM reference for the reader #632 (thanks @SebastianHollizeck)
+ Merged: Fix help messages for python tools #635 (thanks @ibebio)
+ Fixed: Please port to Python3 #569 which removes Debian patch (thanks @tillea)
+ Moved unused files from ./src to ./contrib/freebayes
+ Moved stripped SeqLib source into ./contrib/SeqLib as SeqLib is no
  longer maintained
+ Made htslib a primary dependency of freebayes
+ Updated htslib to the latest version
+ Fixed clang build, requires a patch on htslib sent upstream with https://github.com/samtools/htslib/pull/1190
+ Applied https://github.com/walaj/SeqLib/pull/53 (thanks @jmarshall)
+ Freebayes builds against stock htslib and htslib sources are removed
  from the tree! Moved them back in again for Debian builds and codecs
  support
+ Updated travis-CI; passes new meson builds for gcc

Some minor stuff:

+ Added regression test so we track some changes through git repo
+ Added github agent to deal with stale issues on the issue tracker

## Earlier release notes

See github [releases](https://github.com/freebayes/freebayes/releases)
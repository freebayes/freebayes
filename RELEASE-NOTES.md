For contributions
see
[contributors](https://github.com/freebayes/freebayes/graphs/contributors)
and
[commits](https://github.com/freebayes/freebayes/commits/master).


## ChangeLog v1.3.5 (20210210)

This is a maintenance and bug-fix release of Freebayes:

+ Namespace Fasta* classes which has a build conflict with fastahack
  https://github.com/freebayes/freebayes/commit/00594a1db84855f99c9fe2c23305b86fa542dff2
  (thanks @jnumm)

## ChangeLog v1.3.4 (20210129)

This is a maintenance and bug-fix release of Freebayes:

+ 9 issues [closed](https://github.com/freebayes/freebayes/milestone/1?closed=1)
+ Added support for --trim-complex-tail with #139 and 97735089a1bbf658862fb16f4514c0ad93195e0a (thanks @tsibley)
+ Fixed meson to not build dependencies if they are found to be installed. This should facilitate Debian builds (thanks @jnumm)
+ Cleanup and small fixes (thanks @jnumm)
+ Tested running combinations of builds with installed libs and local builds (htslib, vcflib, tabixpp --- only seqlib needs confirmation)
+ Updated scripts to use python3 (thanks @tillea)

## ChangeLog v1.3.3 (20201213)

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
+ Freebayes builds against stock htslib for Debian etc. Unfortunately
  Debian does not include htslib headers, nor does it handle LZMA
  codecs for CRAM, so the source tree is pulled in for local
  builds. See #664.
+ Updated travis-CI: passes new meson builds for gcc
+ Added github-CI: passes

Some minor stuff:

+ Added regression test so we track some changes through git repo
+ Added github agent to deal with stale issues on the issue tracker

## Earlier release notes

See github [releases](https://github.com/freebayes/freebayes/releases)

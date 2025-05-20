# Meson build definition
#
#     meson --buildtype [debug|release] [--default-library static] build/
#     meson test -C build /
#
# To compile with local git module sources use
#
#     meson build -Dprefer_system_deps=false
#
# to build static binary for release
#
#      meson setup build/ -Dstatic=true -Dprefer_system_deps=false --buildtype release
#
# or for debugging
#
#      meson setup build/ --buildtype debug --wipe -Dprefer_system_deps=false -Dstatic=true
#
# See also meson_options.txt and contrib/README.md
#
# We use guix.scm for development -- see the header of the guix.scm file

project('freebayes', ['cpp', 'c'],
        version : '1.3.10',
        default_options : ['warning_level=1', 'cpp_std=c++17', 'optimization=3'])
add_global_arguments('-Wno-maybe-uninitialized', language : 'cpp') # for ttmath warnings

static_build = get_option('static')

# used by Debian: add_global_arguments('-D_GLIBCXX_ASSERTIONS=1', language : 'cpp')

cc = meson.get_compiler('cpp')

lzma_dep = dependency('liblzma', static: static_build)
thread_dep = dependency('threads', static: static_build)
wfa2lib_dep = cc.find_library('wfa2') # need to link for vcflib code
zlib_dep = dependency('zlib', static: static_build)
htslib_dep = dependency('htslib', static: static_build, required: false)
vcflib_dep = cc.find_library('vcflib', required: false, static: static_build)
seqlib_dep = dependency('libseqlib', static: static_build, required: false)
tabixpp_dep = cc.find_library('tabixpp', required: false, static: static_build)

# later versions of vcflib have Variant.hpp in a system subdirectory
fs = import('fs')
if fs.is_dir('/usr/include/vcflib')
  message('Using /usr/include/vcflib dir')
  vcflib_inc = include_directories('/usr/include/vcflib')
else
  vcflib_inc = files()
endif

if get_option('prefer_system_deps') # defaults to true, otherwise use git modules and local sources
  fastahack_dep = cc.find_library('libfastahack', required: true) # also link in for vcflib
  fastahack_src = files()
  smithwaterman_dep = cc.find_library('libsmithwaterman')
  if fs.is_dir('/usr/include/fastahack')
    message('Using /usr/include/fastahack dir')
    fastahack_inc = include_directories('/usr/include/fastahack')
  else
    fastahack_inc = files()
  endif
else
  # uses the minimal local git submodules or checkout trees in ./contrib/
  # see `git submodule`
  # htslib and vcflib have to come from the underlying distro
  fastahack_dep = dependency('', required : false) # from this tree
  smithwaterman_dep = dependency('', required : false) # from this tree

  fastahack_inc = include_directories('contrib/fastahack')
  fastahack_src = files(
      'contrib/fastahack/Fasta.cpp',
      'contrib/fastahack/split.cpp',
  )

  smithwaterman_inc = include_directories('contrib/smithwaterman')
  smithwaterman_src = files(
      'contrib/smithwaterman/SmithWatermanGotoh.cpp',
      'contrib/smithwaterman/disorder.cpp',
      'contrib/smithwaterman/Repeats.cpp',
      'contrib/smithwaterman/LeftAlign.cpp',
      'contrib/smithwaterman/IndelAllele.cpp',
  )

endif # build git modules

# for setting a warning_level on the external code in custom_* targets below
warn_quiet = ['warning_level=0']

if not seqlib_dep.found() # Seqlib is not in Guix (yet) but it is in Debian
  seqlib_inc = include_directories(
      'contrib/SeqLib',
  )
  seqlib_src = files(
      'contrib/SeqLib/src/BamReader.cpp',
      'contrib/SeqLib/src/BamRecord.cpp',
      'contrib/SeqLib/src/BamHeader.cpp',
      'contrib/SeqLib/src/BamWriter.cpp',
      'contrib/SeqLib/src/GenomicRegion.cpp',
      'contrib/SeqLib/src/ssw_cpp.cpp',
      'contrib/SeqLib/src/ssw.c',
  )

  seqlib_lib = static_library('custom_seqlib',
                              seqlib_src,
                              include_directories : [seqlib_inc],
                              override_options : warn_quiet)
  seqlib_dep = declare_dependency(link_with : seqlib_lib,
                                  include_directories : seqlib_inc)
endif

#
# Sources
#
freebayes_common_src = files(
    'src/Allele.cpp',
    'src/AlleleParser.cpp',
    'src/BedReader.cpp',
    'src/Bias.cpp',
    'src/CNV.cpp',
    'src/Contamination.cpp',
    'src/DataLikelihood.cpp',
    'src/Dirichlet.cpp',
    'src/Ewens.cpp',
    'src/FBFasta.cpp',
    'src/Genotype.cpp',
    'src/IndelAllele.cpp',
    'src/LeftAlign.cpp',
    'src/Marginals.cpp',
    'src/Multinomial.cpp',
    'src/NonCall.cpp',
    'src/Parameters.cpp',
    'src/Result.cpp',
    'src/ResultData.cpp',
    'src/Sample.cpp',
    'src/SegfaultHandler.cpp',
    'src/Utility.cpp',
    )
freebayes_src = files('src/freebayes.cpp')
bamleftalign_src = files('src/bamleftalign.cpp')

# Include paths
incdir = include_directories(
    'src',
    'contrib/ttmath',
    )

extra_cpp_args = cc.get_supported_arguments(
    '-fpermissive',
    '-Wno-reorder',
    '-Wno-sign-compare',
    '-Wno-unused-variable',
    '-Wno-unused-but-set-variable',
    )

freebayes_lib = static_library(
    'freebayes_common',
    [freebayes_common_src, fastahack_src],
    # pic: true,
    include_directories : [incdir, vcflib_inc, fastahack_inc],
    cpp_args : extra_cpp_args,
    dependencies : [vcflib_dep, zlib_dep, lzma_dep, thread_dep, htslib_dep, tabixpp_dep,
                    wfa2lib_dep, seqlib_dep],
    install : false,
    )

if static_build
  link_arguments = '-static'
else
  link_arguments = []
endif

executable('freebayes',
           [freebayes_src],
           include_directories : [incdir, vcflib_inc, fastahack_inc],
           cpp_args : extra_cpp_args,
           link_args: link_arguments,
           dependencies: [zlib_dep,
                          lzma_dep,
                          wfa2lib_dep,
                          vcflib_dep,
                          thread_dep,
                          htslib_dep,
                          tabixpp_dep,
                          seqlib_dep,
                          fastahack_dep, # for vcflib
                          smithwaterman_dep # for vcflib
                          ],
           link_with : [freebayes_lib],
           install: true
          )

executable('bamleftalign',
           [bamleftalign_src],
           include_directories : [incdir, vcflib_inc, fastahack_inc],
           cpp_args : extra_cpp_args,
           link_args: link_arguments,
           dependencies: [zlib_dep,
                          lzma_dep,
                          thread_dep,
                          htslib_dep,
                          vcflib_dep,
                          seqlib_dep,
                          fastahack_dep, # for vcflib
                          smithwaterman_dep, # for vcflib
                          wfa2lib_dep],
           link_with : freebayes_lib,
           install: true
          )

testdir = meson.current_source_dir()+'/test'

prove = find_program('prove')
test('T00', prove, args : ['-e','bash','-v','t/00_region_and_target_handling.t'], workdir: testdir)
test('T01', prove, args : ['-e','bash','-v','t/01_call_variants.t'], workdir : testdir )
test('T01b', prove, args : ['-e','bash','-v','t/01b_call_variants.t'], workdir : testdir )
test('T02', prove, args : ['-e','bash','-v','t/02_multi_bam.t'], workdir : testdir )
test('T03', prove, args : ['-e','bash','-v','t/03_reference_bases.t'], workdir: testdir )

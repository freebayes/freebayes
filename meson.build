# Meson build definition
#
#     meson --buildtype [debug|release] [--default-library static] build/
#     meson test -C build /
#
# Meson builds with locally installed htslib, vcflib, tabixpp and seqlib.
# If one is missing it should pick up the sources from git submodules.
#
# to compile htslib with local sources use
#
#     meson build -Dprefer_system_deps=false
#
# to build static binary for release
#
#     meson build -Dstatic=true -Dprefer_system_deps=false --buildtype release

project('freebayes', ['cpp', 'c'],
        version : '1.3.7-pre1',
        default_options : ['warning_level=1', 'cpp_std=c++17', 'optimization=3'])
static_build = get_option('static')

add_global_arguments('-D_GLIBCXX_ASSERTIONS=1', language : 'cpp')
cc = meson.get_compiler('cpp')

# fastahack_dep = cc.find_library('libfastahack')
lzma_dep = dependency('liblzma', static: static_build)
seqlib_dep = dependency('libseqlib', static: static_build, required: false)
# smithwaterman_dep = cc.find_library('libsmithwaterman')
tabixpp_dep = cc.find_library('tabixpp', required: false, static: static_build)
thread_dep = dependency('threads', static: static_build)
vcflib_dep = cc.find_library('libvcflib')
wfa2lib_dep = cc.find_library('libwfa2')
zlib_dep = dependency('zlib', static: static_build)

# to compile htslib use
#   meson build -Dprefer_system_deps=false
# otherwise it tries to pick up the local system libs
if get_option('prefer_system_deps')
  htslib_dep = dependency('htslib', static: static_build, required: false)
else
  # uses the local git submodules
  htslib_dep = dependency('', required : false)
  tabixpp_dep = dependency('', required : false)
  vcflib_dep = dependency('', required : false)
  seqlib_dep = dependency('', required : false)
endif

# for setting a warning_level on the external code in custom_* targets below
warn_quiet = ['warning_level=0']

# We pull in the include files because they are currently not packaged in Debian, note that we expect
# the share libs nevertheless
fastahack_inc = include_directories('contrib/fastahack')
fastahack_src = files(
    'contrib/fastahack/Fasta.cpp',
    )

multichoose_inc = include_directories('contrib/multichoose')
multichoose_src = files(
    )

smithwaterman_inc = include_directories('contrib/smithwaterman')
smithwaterman_src = files(
    'contrib/smithwaterman/SmithWatermanGotoh.cpp',
    'contrib/smithwaterman/disorder.cpp',
    'contrib/smithwaterman/Repeats.cpp',
    'contrib/smithwaterman/LeftAlign.cpp',
    'contrib/smithwaterman/IndelAllele.cpp',
    )

# optionally build from source. As we unbundled just check out the vcflib repo in contrib.
# note that the build version adds support for lzma and bzip2 in CRAM
if not htslib_dep.found()
  htslib_inc = include_directories(
    'contrib/vcflib/contrib/tabixpp/htslib/',
    'contrib/htslib_config',
    'contrib/vcflib/contrib/tabixpp/htslib/cram',
    )
  htslib_src = files(
    'contrib/vcflib/contrib/tabixpp/htslib/tbx.c',      # used by vcflib
    'contrib/vcflib/contrib/tabixpp/htslib/bgzf.c',
    'contrib/vcflib/contrib/tabixpp/htslib/sam.c',
    'contrib/vcflib/contrib/tabixpp/htslib/hts.c',
    'contrib/vcflib/contrib/tabixpp/htslib/hfile.c',
    'contrib/vcflib/contrib/tabixpp/htslib/textutils.c',
    'contrib/vcflib/contrib/tabixpp/htslib/thread_pool.c',
    'contrib/vcflib/contrib/tabixpp/htslib/region.c',
    'contrib/vcflib/contrib/tabixpp/htslib/hts_os.c',
    'contrib/vcflib/contrib/tabixpp/htslib/cram/cram_index.c',
    'contrib/vcflib/contrib/tabixpp/htslib/cram/cram_io.c',
    'contrib/vcflib/contrib/tabixpp/htslib/header.c',
    'contrib/vcflib/contrib/tabixpp/htslib/cram/pooled_alloc.c',
    'contrib/vcflib/contrib/tabixpp/htslib/cram/string_alloc.c',
    'contrib/vcflib/contrib/tabixpp/htslib/md5.c',
    'contrib/vcflib/contrib/tabixpp/htslib/cram/cram_encode.c',
    'contrib/vcflib/contrib/tabixpp/htslib/cram/cram_stats.c',
    'contrib/vcflib/contrib/tabixpp/htslib/cram/cram_codecs.c',
    'contrib/vcflib/contrib/tabixpp/htslib/cram/cram_decode.c',
    'contrib/vcflib/contrib/tabixpp/htslib/kstring.c',
    'contrib/vcflib/contrib/tabixpp/htslib/cram/cram_external.c',
    'contrib/vcflib/contrib/tabixpp/htslib/cram/mFILE.c',
    'contrib/vcflib/contrib/tabixpp/htslib/faidx.c',
    'contrib/vcflib/contrib/tabixpp/htslib/hts_expr.c',
    'contrib/vcflib/contrib/tabixpp/htslib/htscodecs/htscodecs/rANS_static.c',
    'contrib/vcflib/contrib/tabixpp/htslib/htscodecs/htscodecs/arith_dynamic.c',
    'contrib/vcflib/contrib/tabixpp/htslib/htscodecs/htscodecs/fqzcomp_qual.c',
    'contrib/vcflib/contrib/tabixpp/htslib/htscodecs/htscodecs/htscodecs.c',
    'contrib/vcflib/contrib/tabixpp/htslib/htscodecs/htscodecs/rANS_static.c',
    'contrib/vcflib/contrib/tabixpp/htslib/htscodecs/htscodecs/tokenise_name3.c',
    'contrib/vcflib/contrib/tabixpp/htslib/htscodecs/htscodecs/utils.c',
    'contrib/vcflib/contrib/tabixpp/htslib/cram/open_trace_file.c',
    'contrib/vcflib/contrib/tabixpp/htslib/multipart.c',
    'contrib/vcflib/contrib/tabixpp/htslib/htscodecs/htscodecs/rANS_static4x16pr.c',
    'contrib/vcflib/contrib/tabixpp/htslib/htscodecs/htscodecs/rANS_static32x16pr.c',
    'contrib/vcflib/contrib/tabixpp/htslib/htscodecs/htscodecs/pack.c',
    'contrib/vcflib/contrib/tabixpp/htslib/htscodecs/htscodecs/rle.c',
)
    htslib_lib = static_library('custom_htslib',
                                htslib_src,
                                include_directories : htslib_inc,
                                override_options : warn_quiet)
    htslib_dep = declare_dependency(link_with : htslib_lib,
                                    include_directories : htslib_inc)
else
    htslib_inc = []
endif


if not tabixpp_dep.found()
    tabixpp_inc = include_directories('contrib/vcflib/contrib/tabixpp/')
    tabixpp_src = files('contrib/vcflib/contrib/tabixpp/tabix.cpp')
    tabixpp_lib = static_library('custom_tabixpp',
                                 tabixpp_src,
                                 include_directories : [tabixpp_inc, htslib_inc],
                                 override_options : warn_quiet)
    tabixpp_dep = declare_dependency(link_with : tabixpp_lib,
                                     include_directories : tabixpp_inc)
else
    tabixpp_inc = []
endif

vcflib_inc = include_directories()

# optionally build from source. As we unbundled just check out the vcflib repo in contrib.
if not vcflib_dep.found()
  vcflib_inc = include_directories(
    'contrib/vcflib/src',
    'contrib/vcflib/contrib/smithwaterman',
    'contrib/vcflib/contrib/multichoose',
    'contrib/vcflib/contrib/filevercmp',
    'contrib/vcflib/contrib/fastahack',
    )
  vcflib_src = files(
    'contrib/vcflib/src/Variant.cpp',
    'contrib/vcflib/src/allele.cpp',
    'contrib/vcflib/src/cigar.cpp',
    'contrib/vcflib/src/split.cpp',
    )
  vcflib_lib = static_library('custom_vcflib',
                              vcflib_src,
                              include_directories : [htslib_inc, tabixpp_inc, vcflib_inc],
                              override_options : warn_quiet)
  vcflib_dep = declare_dependency(link_with : vcflib_lib,
                                  include_directories : vcflib_inc)
endif


if not seqlib_dep.found()

  seqlib_inc = include_directories(
    'contrib',
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
                              include_directories : [htslib_inc, vcflib_inc, seqlib_inc],
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
    freebayes_common_src,
    include_directories : [multichoose_inc, fastahack_inc, smithwaterman_inc, incdir],
    cpp_args : extra_cpp_args,
    dependencies : [zlib_dep, lzma_dep, thread_dep, htslib_dep, tabixpp_dep,
                    vcflib_dep, wfa2lib_dep, seqlib_dep],
    install : false,
    )

if static_build
  link_arguments = '-static'
else
  link_arguments = []
endif

executable('freebayes',
           [freebayes_src,multichoose_src,fastahack_src,smithwaterman_src],
           include_directories : [multichoose_inc, fastahack_inc, smithwaterman_inc, incdir],
           cpp_args : extra_cpp_args,
           link_args: link_arguments,
           dependencies: [zlib_dep,
                          lzma_dep,
                          thread_dep,
                          htslib_dep,
                          tabixpp_dep,
                          vcflib_dep,
                          seqlib_dep,
                          # fastahack_dep,
                          # smithwaterman_dep,
                          wfa2lib_dep],
           link_with : freebayes_lib,
           install: true
          )

executable('bamleftalign',
           [bamleftalign_src,fastahack_src,smithwaterman_src],
           include_directories : [multichoose_inc, fastahack_inc, smithwaterman_inc, incdir],
           cpp_args : extra_cpp_args,
           link_args: link_arguments,
           dependencies: [zlib_dep,
                          lzma_dep,
                          thread_dep,
                          htslib_dep,
                          vcflib_dep,
                          seqlib_dep,
                          # fastahack_dep,
                          # smithwaterman_dep,
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

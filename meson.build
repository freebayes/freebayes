# Meson build definition
#
# meson --buildtype [debug|release] [--default-library static] build/
# meson test -C build /

project('freebayes', ['cpp', 'c'],
        version : '1.3.3',
        default_options : ['warning_level=3', 'cpp_std=c++14'])

zlib_dep = dependency('zlib')
lzma_dep = dependency('liblzma')
bzip2_dep = dependency('bz2lib', required: false)
htslib_dep = dependency('htslib', required : false)
thread_dep = dependency('threads')

if htslib_dep.found()
  htslib_includes = ''
  htslib_src = []
else
  # build from source - necessary on some Debian systems
  # note that the build version adds support for lzma and bzip2 in CRAM
  htslib_includes = [
    'contrib/htslib',
    'contrib/htslib_config',
    'contrib/htslib/cram',
    ]
  htslib_src = [
    'contrib/htslib/tbx.c',      # used by vcflib
    'contrib/htslib/bgzf.c',
    'contrib/htslib/sam.c',
    'contrib/htslib/hts.c',
    'contrib/htslib/hfile.c',
    'contrib/htslib/hfile_net.c',
    'contrib/htslib/knetfile.c',
    'contrib/htslib/textutils.c',
    'contrib/htslib/thread_pool.c',
    'contrib/htslib/region.c',
    'contrib/htslib/hts_os.c',
    'contrib/htslib/cram/cram_index.c',
    'contrib/htslib/cram/cram_io.c',
    # the older 'contrib/htslib/cram/sam_header.c', was renamed to
    'contrib/htslib/header.c',
    'contrib/htslib/cram/pooled_alloc.c',
    'contrib/htslib/cram/string_alloc.c',
    'contrib/htslib/md5.c',
    'contrib/htslib/cram/cram_encode.c',
    'contrib/htslib/cram/cram_stats.c',
    'contrib/htslib/cram/cram_codecs.c',
    'contrib/htslib/cram/cram_decode.c',
    'contrib/htslib/cram/cram_samtools.c',
    'contrib/htslib/kstring.c',
    'contrib/htslib/cram/cram_external.c',
    # 'contrib/htslib/cram/files.c',
    'contrib/htslib/cram/mFILE.c',
    'contrib/htslib/faidx.c',
    'contrib/htslib/cram/rANS_static.c',
    'contrib/htslib/cram/open_trace_file.c',
    'contrib/htslib/multipart.c',
]
endif



#
# Sources
#
freebayes_src = [
    'src/freebayes.cpp',
    'src/Allele.cpp',
    'src/AlleleParser.cpp',
    'src/BedReader.cpp',
    'src/BGZF.cpp',
    'src/Bias.cpp',
    'src/CNV.cpp',
    'src/Contamination.cpp',
    'src/DataLikelihood.cpp',
    'src/Dirichlet.cpp',
    'src/Ewens.cpp',
    'src/Fasta.cpp',
    'src/freebayes.cpp',
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
    'src/split.cpp',
    'src/Utility.cpp'
]



seqlib_src = [
    'contrib/SeqLib/src/BamReader.cpp',
    'contrib/SeqLib/src/BamRecord.cpp',
    'contrib/SeqLib/src/BamHeader.cpp',
    'contrib/SeqLib/src/GenomicRegion.cpp',
    'contrib/SeqLib/src/ssw_cpp.cpp',
    'contrib/SeqLib/src/ssw.c',
]

vcflib_src = [
    'vcflib/tabixpp/tabix.cpp',
    'vcflib/src/Variant.cpp',
    'vcflib/smithwaterman/SmithWatermanGotoh.cpp',
    'vcflib/smithwaterman/disorder.cpp',
    'vcflib/smithwaterman/Repeats.cpp',
    'vcflib/smithwaterman/LeftAlign.cpp',
    'vcflib/smithwaterman/IndelAllele.cpp',
]

bamleftalign_src = [
    'src/bamleftalign.cpp',
    'src/IndelAllele.cpp',
    'contrib/SeqLib/src/BamWriter.cpp',
    'vcflib/fastahack/Fasta.cpp',
    'vcflib/smithwaterman/LeftAlign.cpp',
    'vcflib/smithwaterman/IndelAllele.cpp',
    'vcflib/src/split.cpp',
    'src/LeftAlign.cpp',
]

# Include paths
incdir = include_directories(
    'src',
    htslib_includes,
    'bamtools/src/api',
    'bamtools/src',
    'ttmath',
    'contrib',
    'contrib/SeqLib',
    'vcflib/src',
    'vcflib/tabixpp',
    'vcflib/smithwaterman',
    'vcflib/multichoose',
    'vcflib/filevercmp')

c_args = ['-fpermissive','-w']
cpp_args = ['-fpermissive','-w','-Wc++14-compat']

executable('freebayes',
           freebayes_src,
           seqlib_src,
           vcflib_src,
           htslib_src,
           # htslib_src,
           include_directories : incdir,
           cpp_args : cpp_args,
           c_args :   c_args,
           dependencies: [zlib_dep, lzma_dep, htslib_dep, thread_dep],
           install: true
          )

executable('bamleftalign',
           bamleftalign_src,
           seqlib_src,
           vcflib_src,
           htslib_src,
           # htslib_src,
           include_directories : incdir,
           cpp_args : cpp_args,
           c_args :   c_args,
           dependencies: [zlib_dep, lzma_dep, htslib_dep, thread_dep],
           install: true
          )

testdir = meson.current_source_dir()+'/test'

prove = find_program('prove')
test('T00', prove, args : ['-e','bash','-v','t/00_region_and_target_handling.t'], workdir: testdir)
test('T01', prove, args : ['-e','bash','-v','t/01_call_variants.t'], workdir : testdir )
test('T01b', prove, args : ['-e','bash','-v','t/01b_call_variants.t'], workdir : testdir )
test('T02', prove, args : ['-e','bash','-v','t/02_multi_bam.t'], workdir : testdir )
test('T03', prove, args : ['-e','bash','-v','t/03_reference_bases.t'], workdir: testdir )

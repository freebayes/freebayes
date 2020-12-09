project('freebayes', ['cpp', 'c'],
        version : '1.3.3',
        default_options : ['warning_level=3', 'cpp_std=c++14'])

zlib_dep = dependency('zlib', required: false)

thread_dep = dependency('threads')

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
    'bamtools/src/api',
    'bamtools/src',
    'ttmath',
    'contrib',
    'contrib/SeqLib',
    'contrib/htslib_config',
    'contrib/htslib/cram',
    'contrib/htslib/htslib',
    'contrib/htslib',
    'vcflib/src',
    'vcflib/tabixpp',
    'vcflib/smithwaterman',
    'vcflib/multichoose',
    'vcflib/filevercmp')

# Parsing pkg-config by hand because it does not get honoured
pkgconfig_prog = find_program('pkg-config')
r = run_command(pkgconfig_prog.path(), '--libs', 'zlib')

executable('freebayes',
           freebayes_src,
           seqlib_src,
           vcflib_src,
           htslib_src,
           include_directories : incdir,
           cpp_args : ['-fpermissive','-w'],
           c_args : ['-fpermissive','-w'],
           link_args : '-Wl,'+','.join(r.stdout().strip().split()),
           dependencies: [zlib_dep, thread_dep],
           install: true
          )

executable('bamleftalign',
           bamleftalign_src,
           seqlib_src,
           vcflib_src,
           htslib_src,
           include_directories : incdir,
           cpp_args : ['-w','-fpermissive','-Wc++14-compat'],
           c_args : ['-w','-fpermissive'],
           link_args : '-Wl,'+','.join(r.stdout().strip().split()),
           dependencies: [zlib_dep, thread_dep],
           install: true
          )

prove = find_program('prove')
test('T00', prove, args : ['-e','bash','-v','../test/t/00_region_and_target_handling.t'])
test('T01', prove, args : ['-e','bash','-v','../test/t/01_call_variants.t'])
test('T02', prove, args : ['-e','bash','-v','../test/t/02_multi_bam.t'])
test('T03', prove, args : ['-e','bash','-v','../test/t/03_reference_bases.t'])

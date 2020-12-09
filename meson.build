project('freebayes', ['cpp', 'c'],
        version : '1.0',
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


vcflib_src = [
    'vcflib/tabixpp/tabix.cpp',
    'vcflib/src/Variant.cpp',
    'vcflib/smithwaterman/SmithWatermanGotoh.cpp',
    'vcflib/smithwaterman/disorder.cpp',
    'vcflib/smithwaterman/Repeats.cpp',
    'vcflib/smithwaterman/LeftAlign.cpp',
    'vcflib/smithwaterman/IndelAllele.cpp',
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
    'htslib/tbx.c',      # used by vcflib
    'htslib/bgzf.c',
    'htslib/sam.c',
    'htslib/hts.c',
    'htslib/hfile.c',
    'htslib/hfile_net.c',
    'htslib/knetfile.c',
    'htslib/textutils.c',
    'htslib/thread_pool.c',
    'htslib/region.c',
    'htslib/hts_os.c',
    'htslib/cram/cram_index.c',
    'htslib/cram/cram_io.c',
    # the older 'htslib/cram/sam_header.c', was renamed to
    'htslib/header.c',
    'htslib/cram/pooled_alloc.c',
    'htslib/cram/string_alloc.c',
    'htslib/md5.c',
    'htslib/cram/cram_encode.c',
    'htslib/cram/cram_stats.c',
    'htslib/cram/cram_codecs.c',
    'htslib/cram/cram_decode.c',
    'htslib/cram/cram_samtools.c',
    'htslib/kstring.c',
    'htslib/cram/cram_external.c',
    # 'htslib/cram/files.c',
    'htslib/cram/mFILE.c',
    'htslib/faidx.c',
    'htslib/cram/rANS_static.c',
    'htslib/cram/open_trace_file.c',
    'htslib/multipart.c',
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

# Include paths - there are some duplicates in there
incdir = include_directories(
    'src',
    'bamtools/src/api',
    'bamtools/src',
    'ttmath',
    'contrib/SeqLib',
    'htslib/cram',
    'htslib/htslib',
    'htslib',
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
           link_args : '-Wl,'+','.join(r.stdout().strip().split()),
           dependencies: [zlib_dep, thread_dep]
          )

# executable('bamleftalign',
#            bamleftalign_src,
#            seqlib_src,
#            vcflib_src,
#            htslib_src,
#            include_directories : incdir,
#            cpp_args : '-fpermissive',
#            link_args : '-Wl,'+','.join(r.stdout().strip().split()),
#            dependencies: [zlib_dep, thread_dep]
#           )

prove = find_program('prove')
test('T00', prove, args : ['-e','bash','-v','../test/t/00_region_and_target_handling.t'])
test('T01', prove, args : ['-e','bash','-v','../test/t/01_call_variants.t'])
test('T02', prove, args : ['-e','bash','-v','../test/t/02_multi_bam.t'])
test('T03', prove, args : ['-e','bash','-v','../test/t/03_reference_bases.t'])

project('freebayes', language: ['cpp', 'c'],
        version : '1.0',
        default_options : ['warning_level=3', 'cpp_std=c++14'])

zlib_dep = dependency(
        'zlib', required : false)
# conf.set10('HAVE_LIBZ', libz.found())

thread_dep = dependency('threads')

#
# Sources
#
freebayes_src = [

    'src/freebayes.cpp',
    'src/Allele.cpp',
    'src/AlleleParser.cpp',
    # 'src/alleles.cpp',
    # 'src/bamfiltertech.cpp',
    # 'src/bamleftalign.cpp',
    'src/BedReader.cpp',
    'src/BGZF.cpp',
    'src/Bias.cpp',
    'src/CNV.cpp',
    'src/Contamination.cpp',
    'src/DataLikelihood.cpp',
    'src/Dirichlet.cpp',
    # 'src/dummy.cpp',
    'src/Ewens.cpp',
    'src/Fasta.cpp',
    'src/freebayes.cpp',
    'src/Genotype.cpp',
    # 'src/GenotypePriors.cpp',
    'src/IndelAllele.cpp',
    'src/LeftAlign.cpp',
    # 'src/levenshtein.cpp',
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
    'SeqLib/htslib/tbx.c',      # used by vcflib
    'vcflib/tabixpp/tabix.cpp',
    'vcflib/src/Variant.cpp',
    'vcflib/smithwaterman/SmithWatermanGotoh.cpp',
    'vcflib/smithwaterman/disorder.cpp',
    'vcflib/smithwaterman/Repeats.cpp',
    'vcflib/smithwaterman/LeftAlign.cpp',
    'vcflib/smithwaterman/IndelAllele.cpp',
]

seqlib_src = [
    'SeqLib/src/BamReader.cpp',
    'SeqLib/src/BamRecord.cpp',
    'SeqLib/src/BamHeader.cpp',
    'SeqLib/src/GenomicRegion.cpp',
    'SeqLib/htslib/bgzf.c',
    'SeqLib/htslib/sam.c',
    'SeqLib/htslib/hts.c',
    'SeqLib/htslib/hfile.c',
    'SeqLib/htslib/hfile_net.c',
    'SeqLib/htslib/knetfile.c',
    'SeqLib/htslib/textutils.c',
    'SeqLib/htslib/thread_pool.c',
    'SeqLib/htslib/cram/cram_index.c',
    'SeqLib/htslib/cram/cram_io.c',
    'SeqLib/htslib/cram/sam_header.c',
    'SeqLib/htslib/cram/pooled_alloc.c',
    'SeqLib/htslib/cram/string_alloc.c',
    'SeqLib/htslib/md5.c',
    'SeqLib/htslib/cram/cram_encode.c',
    'SeqLib/htslib/cram/cram_stats.c',
    'SeqLib/htslib/cram/cram_codecs.c',
    'SeqLib/htslib/cram/cram_decode.c',
    'SeqLib/htslib/cram/cram_samtools.c',
    'SeqLib/htslib/kstring.c',
    'SeqLib/htslib/cram/cram_external.c',
    'SeqLib/htslib/cram/files.c',
    'SeqLib/htslib/cram/mFILE.c',
    'SeqLib/htslib/faidx.c',
    'SeqLib/htslib/cram/rANS_static.c',
    'SeqLib/htslib/cram/open_trace_file.c',
    'SeqLib/htslib/multipart.c',
    'SeqLib/src/ssw_cpp.cpp',
    'SeqLib/src/ssw.c',

]

# Include paths - there are some duplicates in there
incdir = include_directories(
    'src',
    'bamtools/src/api',
    'bamtools/src',
    'ttmath',
    'SeqLib',
    'SeqLib/htslib/cram',
    'SeqLib/htslib',
    'vcflib/src',
    'vcflib/tabixpp',
    'vcflib/smithwaterman',
    'vcflib/multichoose',
    'vcflib/filevercmp')

executable('freebayes',
           freebayes_src,
           seqlib_src,
           vcflib_src,
           include_directories : incdir,
           cpp_args : '-Wfatal-errors',
           link_args : '-Wl,-L/gnu/store/899vqdi2wxf50lzip42wfxqb588zri5p-profile/lib,-lz',
           dependencies: [zlib_dep, thread_dep]
          )

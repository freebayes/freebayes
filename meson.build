project('freebayes', ['cpp', 'c'],
        version : '1.3.3',
        default_options : ['warning_level=3', 'cpp_std=c++14'])

zlib_dep = dependency('zlib')
htslib_dep = dependency('htslib', required : false)
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
    'vcflib/src',
    'vcflib/tabixpp',
    'vcflib/smithwaterman',
    'vcflib/multichoose',
    'vcflib/filevercmp')

executable('freebayes',
           freebayes_src,
           seqlib_src,
           vcflib_src,
           # htslib_src,
           include_directories : incdir,
           cpp_args : ['-fpermissive','-w'],
           c_args : ['-fpermissive','-w'],
           dependencies: [zlib_dep, htslib_dep, thread_dep],
           install: true
          )

executable('bamleftalign',
           bamleftalign_src,
           seqlib_src,
           vcflib_src,
           # htslib_src,
           include_directories : incdir,
           cpp_args : ['-w','-fpermissive','-Wc++14-compat'],
           c_args : ['-w','-fpermissive'],
           dependencies: [zlib_dep, htslib_dep, thread_dep],
           install: true
          )

testdir = meson.current_source_dir()+'/test'

prove = find_program('prove')
test('T00', prove, args : ['-e','bash','-v','t/00_region_and_target_handling.t'], workdir: testdir)
# test('T01', prove, args : ['-e','bash','-v','../test/t/01_call_variants.t'])
test('T02', prove, args : ['-e','bash','-v','t/02_multi_bam.t'], workdir : testdir )
test('T03', prove, args : ['-e','bash','-v','t/03_reference_bases.t'], workdir: testdir )

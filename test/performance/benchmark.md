# Benchmarking

In this file we keep track of some benchmarking stats over time.

Download the test BAM file from https://www.internationalgenome.org/

    wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00100/alignment/HG00100.chrom20.ILLUMINA.bwa.GBR.low_coverage.20130415.bam

And a reference genome. E.g.

    wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr20.fa.gz

and unpack. Note you need to change the first 'chr20' tag to '20' for this to work. And index with

    samtools faidx chr20.fa

Create BAM index

    samtools index HG00100.chrom20.ILLUMINA.bwa.GBR.low_coverage.20130415.bam

## AMD Ryzen 7 3700X 8-Core Processor

```sh
time ./freebayes-static-1.3.10-b0d8efd -f chr20.fa HG00100.chrom20.ILLUMINA.bwa.GBR.low_coverage.20130415.bam --region 20:60031-229081 > fb-test-1.3.0-1.vcf
real    0m30.446s
user    0m30.274s
sys     0m0.144s
```

similar to the fastest

```sh
time ./freebayes-1.3.6-linux-amd64-static -f chr20.fa HG00100.chrom20.ILLUMINA.bwa.GBR.low_coverage.20130415.bam --region 20:60031-229081 > fb-test-1.3.0-1.vcf

real    0m29.773s
user    0m29.491s
sys     0m0.184s
```

a native compile (znver2) is the same speed, so we ship the generic version:

```
freebayes 1.3.10

  User defined options
    buildtype         : debugoptimized
    prefix            : /gnu/store/7h6d10wbkwzjfc1ili1rvypisbbpbvrl-freebayes-static-git-1.3.10-HEAD.b0d8efd
    c_link_args       : -Wl,-rpath=/gnu/store/7h6d10wbkwzjfc1ili1rvypisbbpbvrl-freebayes-static-git-1.3.10-HEAD.b0d8efd/lib
    cpp_link_args     : -Wl,-rpath=/gnu/store/7h6d10wbkwzjfc1ili1rvypisbbpbvrl-freebayes-static-git-1.3.10-HEAD.b0d8efd/lib
    prefer_system_deps: false
    static            : true
```

## Penguin2 56x Intel(R) Xeon(R) CPU E5-2683 v3 @ 2.00GHz, 256Gb

Note we used in the older runs:

    head -6000 chr20.fa > chr20-6K.fa

First test an older 1.3.0 release:

```sh
time ./freebayes-v1.3.0-1 -f chr20-6K.fa HG00100.chrom20.ILLUMINA.bwa.GBR.low_coverage.20130415.bam > fb-test-1.3.0-1.vcf
real    1m29.208s
user    1m28.696s
sys     0m0.508s
```

First test an older 1.2.0 release:

```sh
time ./.conda/pkgs/freebayes-1.2.0-htslib1.7_0/bin/freebayes -f chr20-6K.fa HG00100.chrom20.ILLUMINA.bwa.GBR.low_coverage.20130415.bam > fb-test-1.2.0-conda.vcf
real    4m31.958s
user    4m31.364s
sys     0m0.584s
```

Found 49265 variants. The output for 1.2.0 is slightly different, so stick with later versions (after 1.3.0) of freebayes!

Latest static build finds 49265 variants:

```sh
time ./freebayes-1.3.6 -f chr20-6K.fa HG00100.chrom20.ILLUMINA.bwa.GBR.low_coverage.20130415.bam > fb-test-1.3.6.vcf
real    1m10.926s
user    1m10.472s
sys     0m0.444s
```

1.3.0 and 1.3.6 give the same results.  The speed difference between 1.3.0 to 1.3.6 is due to htslib improvements.

The clang (9.0.1) + llvm (13.0.0) build has a similar runtime.

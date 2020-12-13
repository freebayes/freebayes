#!/usr/bin/env bash

cd ../test # for meson

BASH_TAP_ROOT=bash-tap
. bash-tap/bash-tap-bootstrap

root=$(dirname $0)/../..

PATH=../build:$root/build:$root/../build:$root/bin:$PATH
PATH=../scripts:$PATH # for freebayes-parallel

plan tests 10


is $(freebayes -f tiny/q.fa tiny/NA12878.chr22.tiny.bam | grep -v "^#" | wc -l) $(freebayes-parallel tiny/q.regions 2 -f tiny/q.fa tiny/NA12878.chr22.tiny.bam | grep -v "^#" | wc -l) "running in parallel makes no difference"

is $(freebayes -f 'tiny/q with spaces.fa' tiny/NA12878.chr22.tiny.bam | grep -v "^#" | wc -l) $(freebayes-parallel 'tiny/q with spaces.regions' 2 -f 'tiny/q\ with\ spaces.fa' tiny/NA12878.chr22.tiny.bam | grep -v "^#" | wc -l) "freebayes handles spaces in file names"

# check input can hand colons in name like the HLA contigs in GRCh38
is $(freebayes -f tiny/hla.fa -@ tiny/hla.vcf.gz -r HLA-DRB1*16:02:01:1-10000 tiny/NA12878.chr22.tiny.hla.bam | grep -v "^#" | cut -f1,2 | grep -P "(\t500$|\t11000$|\t1000$)" | wc -l) 2 "freebayes handles region and variant input even with : in contig names"

is $(freebayes -f splice/1:883884-887618.fa splice/1:883884-887618.bam -F 0.05 -C 1 | grep ^1 | wc -l) 3 "freebayes can handle spliced reads"

is $(freebayes -f tiny/q.fa -F 0.2 tiny/NA12878.chr22.tiny.bam --gvcf | grep '<\*>' | wc -l) 20 "freebayes produces the expected number of lines of gVCF output"

is $(freebayes -f tiny/q.fa -F 0.2 tiny/NA12878.chr22.tiny.bam --gvcf --gvcf-chunk 50 | grep '<\*>' | wc -l) 250 "freebayes produces the expected number of lines of gVCF output"

is $(freebayes -f tiny/q.fa -F 0.2 tiny/NA12878.chr22.tiny.bam --gvcf --gvcf-dont-use-chunk true | grep '<\*>' | wc -l) 12250 "freebayes produces the expected number of lines of gVCF output"

samtools view -h tiny/NA12878.chr22.tiny.bam | sed s/NA12878D_HiSeqX_R1.fastq.gz/222.NA12878D_HiSeqX_R1.fastq.gz/ | sed s/SM:1/SM:2/ >tiny/x.sam
is $(freebayes -f tiny/q.fa -F 0.2 tiny/NA12878.chr22.tiny.bam tiny/x.sam -A <(echo 1 8; echo 2 13) | grep 'AN=21' | wc -l) 19 "the CNV map may be used to specify per-sample copy numbers"
rm -f tiny/x.sam

is $(freebayes -f tiny/q.fa --skip-coverage 30 tiny/NA12878.chr22.tiny.bam | grep -v '^#' | wc -l) 22 "freebayes makes the expected number of calls when capping coverage"

# is $(freebayes -f tiny/q.fa -g 30 tiny/NA12878.chr22.tiny.bam | vcf2tsv | cut -f 8 | tail -n+2 | awk '$1 <= 30 { print }' | wc -l) 22 "all coverage capped calls are below the coverage threshold"

> cnv-map.bed
for i in {1..10}; do
    echo -e "q\t$((i * 2000))\t$((i * 2000 + 2000))\t1\t$((2 + i))" >> cnv-map.bed
done
is $(freebayes -f tiny/q.fa -F 0.2 tiny/NA12878.chr22.tiny.bam -A cnv-map.bed | grep -v "^#" | cut -f2,10 | cut -f1 -d":" |
    awk 'BEGIN{errors = 0} {cn = int($1 / 2000) + 2; obs = gsub(/\/|\|/, "", $2) + 1; if (cn != obs) errors += 1 } END {OFS=":"; print errors,NR}') "0:19" \
    "freebayes correctly uses CNV map with multiple entries"
rm cnv-map.bed

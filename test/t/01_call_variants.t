#!/usr/bin/env bash

BASH_TAP_ROOT=bash-tap
. ./bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for freebayes
PATH=../scripts:$PATH # for freebayes-parallel
PATH=../vcflib/bin:$PATH # for vcf binaries used by freebayes-parallel

plan tests 19

is $(echo "$(comm -12 <(cat tiny/NA12878.chr22.tiny.giab.vcf | grep -v "^#" | cut -f 2 | sort) <(freebayes -f tiny/q.fa tiny/NA12878.chr22.tiny.bam | grep -v "^#" | cut -f 2 | sort) | wc -l) >= 13" | bc) 1 "variant calling recovers most of the GiAB variants in a test region"

by_region=$((for region in \
    q:180-191 \
    q:1002-1013 \
    q:1811-1825 \
    q:1911-1922 \
    q:2344-2355 \
    q:3257-3268 \
    q:4443-4454 \
    q:5003-5014 \
    q:5074-5085 \
    q:5089-5100 \
    q:5632-5646 \
    q:6412-6423 \
    q:8840-8851 \
    q:9245-9265 \
    q:9785-9796 \
    q:10526-10537 \
    q:11255-11266 \
    q:11530-11541 \
    q:12119-12130;
do
    freebayes -f tiny/q.fa tiny/NA12878.chr22.tiny.bam -r $region | grep -v "^#"
done) |wc -l)

at_once=$(freebayes -f tiny/q.fa tiny/NA12878.chr22.tiny.bam | grep -v "^#" | wc -l)

is $by_region $at_once "freebayes produces the same number of calls if targeted per site or called without targets"

cat >targets.bed <<EOF
q	180	191
q	1002	1013
q	1811	1825
q	1911	1922
q	2344	2355
q	3257	3268
q	4443	4454
q	5003	5014
q	5074	5085
q	5089	5100
q	5632	5646
q	6412	6423
q	8840	8851
q	9245	9265
q	9785	9796
q	10526	10537
q	11255	11266
q	11530	11541
q	12119	12130
EOF

is $(freebayes -f tiny/q.fa tiny/NA12878.chr22.tiny.bam -t targets.bed | grep -v "^#" | wc -l) $by_region "a targets bed file can be used with the same effect as running by region"
rm targets.bed


is $(samtools view -u tiny/NA12878.chr22.tiny.bam | freebayes -f tiny/q.fa --stdin | grep -v "^#" | wc -l) \
    $(freebayes -f tiny/q.fa tiny/NA12878.chr22.tiny.bam | grep -v "^#" | wc -l) "reading from stdin or not makes no difference"

is $(samtools view tiny/NA12878.chr22.tiny.bam | wc -l) $(freebayes -f tiny/q.fa tiny/NA12878.chr22.tiny.bam -d 2>&1 | grep ^alignment: | wc -l) "freebayes processes all alignments in input"

# ensure targeting works even when there are no reads
is $(freebayes -f tiny/q.fa -@ tiny/q.vcf.gz tiny/NA12878.chr22.tiny.bam | grep -v "^#" | wc -l) 19 "freebayes correctly handles variant input"

# ensure that positions at which no variants exist get put in the out vcf
is $(freebayes -f tiny/q.fa -@ tiny/q_spiked.vcf.gz tiny/NA12878.chr22.tiny.bam | grep -v "^#" | cut -f1,2 | grep -P "(\t500$|\t11000$|\t1000$)" | wc -l) 3 "freebayes puts required variants in output"

is $(freebayes -f tiny/q.fa -@ tiny/q_spiked.vcf.gz tiny/NA12878.chr22.tiny.bam -l | grep -v "^#" | cut -f1,2 | grep -P "(\t500$|\t11000$|\t1000$)" | wc -l) 3 "freebayes limits calls to input variants correctly"

is $(freebayes -f tiny/q.fa -@ tiny/q.vcf.gz -l tiny/1read.bam | grep -v "^#" | wc -l) 20 "freebayes reports all input variants even when there is no input data"

# check variant input with region specified
is $(freebayes -f tiny/q.fa -@ tiny/q_spiked.vcf.gz -r q:1-10000 tiny/NA12878.chr22.tiny.bam | grep -v "^#" | cut -f1,2 | grep -P "(\t500$|\t11000$|\t1000$)" | wc -l) 2 "freebayes handles region and variant input"

is $(freebayes -f tiny/q.fa -@ tiny/q_spiked.vcf.gz -r q:1-10000 tiny/NA12878.chr22.tiny.bam -l | grep -v "^#" | cut -f1,2 | grep -P "(\t500$|\t1000$)" | wc -l) 2 "freebayes limits to variant input correctly when region is given"

# check variant input when reading from stdin
is $(freebayes -f tiny/q.fa -@ tiny/q_spiked.vcf.gz --stdin < tiny/NA12878.chr22.tiny.bam | grep -v "^#" | cut -f1,2 | grep -P "(\t500$|\t11000$|\t1000$)" | wc -l) 3 "freebayes handles variant input and reading from stdin"

is $(freebayes -f tiny/q.fa -@ tiny/q_spiked.vcf.gz -l --stdin < tiny/NA12878.chr22.tiny.bam | grep -v "^#" | cut -f1,2 | grep -P "(\t500$|\t11000$|\t1000$)" | wc -l) 3 "freebayes limits to variant input when reading from stdin"

is $(freebayes -f tiny/q.fa -@ tiny/q_spiked.vcf.gz -r q:1-10000 -l --stdin < tiny/NA12878.chr22.tiny.bam | grep -v "^#" | cut -f1,2 | grep -P "(\t500$|\t1000$)" | wc -l) 2 "freebayes handles region, stdin, and variant input"

gzip -c tiny/q.fa >tiny/q.fa.gz
cp tiny/q.fa.fai tiny/q.fa.gz.fai
freebayes -f tiny/q.fa.gz -@ tiny/q_spiked.vcf.gz -r q:1-10000 -l - < tiny/NA12878.chr22.tiny.bam >/dev/null 2>/dev/null
is $? 1 "freebayes bails out when given a gzipped or corrupted reference"
rm tiny/q.fa.gz*

is $(freebayes -f tiny/q.fa tiny/NA12878.chr22.tiny.bam | grep -v "^#" | wc -l) $(freebayes-parallel tiny/q.regions 2 -f tiny/q.fa tiny/NA12878.chr22.tiny.bam | grep -v "^#" | wc -l) "running in parallel makes no difference"

#is $(freebayes -f 'tiny/q with spaces.fa' tiny/NA12878.chr22.tiny.bam | grep -v "^#" | wc -l) $(freebayes-parallel 'tiny/q with spaces.regions' 2 -f 'tiny/q with spaces.fa' tiny/NA12878.chr22.tiny.bam | grep -v "^#" | wc -l) "freebayes handles spaces in file names"


is $(freebayes -f splice/1:883884-887618.fa splice/1:883884-887618.bam | grep ^1 | wc -l) 1 "freebayes can handle spliced reads"

is $(freebayes -f tiny/q.fa tiny/NA12878.chr22.tiny.bam --gvcf | grep '<\*>' | wc -l) 20 "freebayes produces the expected number of lines of gVCF output"

is $(freebayes -f tiny/q.fa tiny/NA12878.chr22.tiny.bam --gvcf --gvcf-chunk 50 | grep '<\*>' | wc -l) 245 "freebayes produces the expected number of lines of gVCF output"

#!/usr/bin/env bash

BASH_TAP_ROOT=bash-tap
. ./bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for freebayes

plan tests 5

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
q 180 191
q 1002 1013
q 1811 1825
q 1911 1922
q 2344 2355
q 3257 3268
q 4443 4454
q 5003 5014
q 5074 5085
q 5089 5100
q 5632 5646
q 6412 6423
q 8840 8851
q 9245 9265
q 9785 9796
q 10526 10537
q 11255 11266
q 11530 11541
q 12119 12130
EOF

is $(freebayes -f tiny/q.fa tiny/NA12878.chr22.tiny.bam -t targets.bed | grep -v "^#" | wc -l) $by_region "a targets bed file can be used with the same effect as running by region"
#rm targets.bed


is $(samtools view -u tiny/NA12878.chr22.tiny.bam | freebayes -f tiny/q.fa --stdin | grep -v "^#" | wc -l) \
    $(freebayes -f tiny/q.fa tiny/NA12878.chr22.tiny.bam | grep -v "^#" | wc -l) "reading from stdin or not makes no difference"

is $(samtools view tiny/NA12878.chr22.tiny.bam | wc -l) $(freebayes -f tiny/q.fa tiny/NA12878.chr22.tiny.bam -d 2>&1 | grep ^alignment: | wc -l) "freebayes processes all alignments in input"

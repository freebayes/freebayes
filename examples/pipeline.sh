#!/bin/bash

outdir=OUTPUT
mkdir -p $outdir
reference=REFERENCE.fasta
bamlist=BAMFILENAMES.txt
cnvmap=CNVMAP.bed
region=REGION

time /share/home/erik/bin/bamtools merge -region $region \
    $(for file in $(cat $bamlist); do echo " -in "$file; done) \
 | time /share/home/erik/bin/ogap -z -R 25 -C 20 -Q 20 -S 0 -f $reference \
 | time /share/home/erik/bin/bamleftalign -f $reference \
 | time /share/home/erik/bin/samtools calmd -EAru - $reference 2>/dev/null \
 | time /share/home/erik/bin/freebayes \
    --min-alternate-count 2 \
    --min-alternate-qsum 40 \
    --pvar 0.0001 \
    --use-mapping-quality \
    --posterior-integration-limits 1,3 \
    --genotype-variant-threshold 4 \
    --site-selection-max-iterations 3 \
    --genotyping-max-iterations 25 \
    --max-complex-gap 3 \
    --cnv-map $cnvmap \
    --stdin \
    --region $region \
    -f $reference \
 | gzip >$outdir/$region.vcf.gz

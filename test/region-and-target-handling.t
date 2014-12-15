#!/bin/bash
# vi: set ft=sh :
test=$(dirname $0)
root=$(dirname $0)/..
source $test/test-simple-bash/lib/test-simple.bash \
    tests 8

PATH=$root/bin:$PATH

ref=$test/`basename $0`.ref
alt=$test/`basename $0`.alt
bam=$test/`basename $0`.bam
bed=$test/`basename $0`.bed
vcf=$test/`basename $0`.vcf

trap 'rm -f $ref* $alt $bam* $bed $vcf' EXIT

# 01234567	start
# ATCGGCTA
# 12345678	end

cat >$ref <<REF
>ref
ATCGGCTA
REF
samtools faidx $ref

cat >$alt <<ALT
>alt
GTTAGGTT
ALT

samtools view -b - >$bam <<SAM
@HD	VN:1.5	SO:coordinate
@SQ	SN:ref	LN:8
alt	0	ref	1	30	1X1=2X1=1X1=1X	*	0	0	GTTAGGTT	*
SAM
samtools index $bam

cat >$bed <<BED
ref	0	1	first_base
ref	2	4	third_and_fourth_base
ref	5	6	sixth_base
ref	7	8	eigth_base
BED

cat >$vcf <<VCF
##fileformat=VCFv4.1
##INFO=<ID=NAME,Number=0,Type=String,Description="Test name">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
ref	1	.	A	G	1234	PASS	NAME=first base
ref	2	.	T	.	1234	PASS	NAME=second base
ref	3	.	C	T	1234	PASS	NAME=third base
ref	4	.	G	A	1234	PASS	NAME=fourth base
ref	5	.	G	.	1234	PASS	NAME=fifth base
ref	6	.	C	G	1234	PASS	NAME=sixth base
ref	7	.	T	.	1234	PASS	NAME=seventh base
ref	8	.	A	T	1234	PASS	NAME=eigth base
VCF

PS4='\n+ '

function run_freebayes() {
    ($root/bin/freebayes "$@" \
        --haplotype-length 0 --min-alternate-count 1 \
        --min-alternate-fraction 0 --pooled-continuous --report-monomorphic \
        --ploidy 1 \
        -f $ref $bam \
        | grep -vE '^#' | cut -f1-5)
}

if [[ -n $TEST_DEBUG ]]; then
    cat $ref >&2
    cat $bed >&2
    cat $vcf >&2
    vcfannotate --bed $bed --key MATCH $vcf >&2
    vcfintersect --bed $bed $vcf >&2
    bedtools intersect -a $vcf -b $bed >&2
fi

[[ -z $(run_freebayes --region ref:4-5 --region ref:6-7) ]]; ok $? 'ref:4-5 ref:6-7 are empty'
[[ -z $(run_freebayes --region ref:4   --region ref:6)   ]]; ok $? 'ref:4   ref:6   are empty'

expected=`cat <<END
ref	6	.	C	G
ref	8	.	A	T
END
`
[[ $(run_freebayes --region ref:5-6 --region ref:7-8) == $expected ]]; ok $? 'ref:5-6 ref:7-8'
[[ $(run_freebayes --region ref:5   --region ref:7)   == $expected ]]; ok $? 'ref:5   ref:7'
[[ $(run_freebayes --region ref:5-)                   == $expected ]]; ok $? 'ref:5-'

expected=`cat <<END
ref	1	.	A	G
ref	3	.	CG	TA
ref	6	.	C	G
ref	8	.	A	T
END
`
[[ $(run_freebayes --targets $bed) == $expected ]]; ok $? "--targets $bed"
[[ $(run_freebayes --region ref)   == $expected ]]; ok $? "--region ref"

[[ $(run_freebayes --region ref:1-20 2>&1) =~ "Target region coordinates (ref 1 19) outside of reference sequence bounds (ref 8)" ]]
    ok $? 'region outside of bounds error'

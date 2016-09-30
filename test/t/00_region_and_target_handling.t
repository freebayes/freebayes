#!/usr/bin/env bash
# vi: set ft=sh :

test=$(dirname $0)/..
root=$(dirname $0)/../..

source $test/test-simple-bash/lib/test-simple.bash \
    tests 11

PATH=$root/bin:$PATH

ref=$test/$(basename $0).ref
alt=$test/$(basename $0).alt
bam=$test/$(basename $0).bam
bed=$test/$(basename $0).bed
vcf=$test/$(basename $0).vcf

trap 'rm -f $ref* $alt $bam* $bed $vcf' EXIT

# 01234567890	start
# ATCGGCTAAAA	ref
# GTTAGGTTAAC	alt
# 12345678901	end

cat >$ref <<REF
>ref
ATCGGCTAAAA
REF
samtools faidx $ref

cat >$alt <<ALT
>alt
GTTAGGTTAAC
ALT

samtools view -S -b - >$bam <<SAM
@HD	VN:1.5	SO:coordinate
@SQ	SN:ref	LN:11
alt	0	ref	1	30	1X1=2X1=1X1=1X2=1X	*	0	0	GTTAGGTTAAC	*
SAM
samtools index $bam

cat >$bed <<BED
ref	0	1	first_base
ref	2	4	third_and_fourth_base
ref	5	6	sixth_base
ref	7	8	eighth_base
ref	10	11	eleventh base
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
ref	8	.	A	T	1234	PASS	NAME=eighth base
ref	9	.	A	.	1234	PASS	NAME=ninth base
ref	10	.	A	.	1234	PASS	NAME=tenth base
ref	11	.	A	C	1234	PASS	NAME=eleventh base
VCF

PS4='\n+ '

function run_freebayes() {
    ($root/bin/freebayes "$@" \
        --haplotype-length 0 --min-alternate-count 1 \
        --min-alternate-fraction 0 --pooled-continuous --report-monomorphic \
        --ploidy 1 \
        -f $ref $bam \
        | grep -vE "^#" | cut -f1-5)
}

if [[ -n $TEST_DEBUG ]]; then
    cat $ref >&2
    cat $bed >&2
    cat $vcf >&2
    vcfannotate --bed $bed --key MATCH $vcf >&2
    vcfintersect --bed $bed $vcf >&2
    bedtools intersect -a $vcf -b $bed >&2
fi


output=$(run_freebayes --region ref:4-5 --region ref:6-7)
ok [ -z "$output" ] "ref:4-5 ref:6-7 are empty" || echo "$output"

output=$(run_freebayes --region ref:4 --region ref:6)
ok [ -z "$output" ] "ref:4 ref:6 are empty" || echo "$output"

expected=$(cat <<END
ref	6	.	C	G
ref	8	.	A	T
END
)
output=$(run_freebayes --region ref:5-6 --region ref:7-8)
ok [ "$output" == "$expected" ] "ref:5-6 ref:7-8" || echo "$output"

output=$(run_freebayes --region ref:5 --region ref:7)
ok [ "$output" == "$expected" ] "ref:5 ref:7" || echo "$output"

expected=$(cat <<END
ref	6	.	C	G
ref	8	.	A	T
ref	11	.	A	C
END
)
output=$(run_freebayes --region ref:5-)
ok [ "$output" == "$expected" ] "ref:5-" || echo "$output"

expected=$(cat <<END
ref	1	.	A	G
ref	3	.	CG	TA
ref	6	.	C	G
ref	8	.	A	T
ref	11	.	A	C
END
)
output=$(run_freebayes)
ok [ "$output" == "$expected" ] "full output" || echo "$output"

output=$(run_freebayes --targets $bed)
ok [ "$output" == "$expected" ] "--targets $bed" || echo "$output"

output=$(run_freebayes --region ref)
ok [ "$output" == "$expected" ] "--region ref" || echo "$output"

# commas as separators: should really be 1,000
expected=$(cat <<END
ref	11	.	A	C
END
)
output=$(run_freebayes --region ref:1,0)
ok [ "$output" == "$expected" ] "ref:1,0-" || echo "$output"

output=$(run_freebayes --region ref:1,0-1,1)
ok [ "$output" == "$expected" ] "ref:1,0-1,1" || echo "$output"

expected="ERROR(freebayes): Target region coordinates (ref 1 20) outside of reference sequence bounds (ref 11) terminating."
output=$(run_freebayes --region ref:1-20 2>&1)
ok [ "$output" == "$expected" ] "region outside of bounds error" || echo "$output"

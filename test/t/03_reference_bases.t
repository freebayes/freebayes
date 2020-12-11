#!/usr/bin/env bash

cd ../test

BASH_TAP_ROOT=bash-tap
source ./bash-tap/bash-tap-bootstrap

root=$(dirname $0)/../..

PATH=../build:$root/build:$root/../build:$root/bin:$PATH

echo "******************"
echo `which freebayes`
ldd `which freebayes`

plan tests 3

ref=$(basename $0).ref
bam=$(basename $0).bam

trap 'rm -f ${ref}* ${bam}*' EXIT

function make_ref() {
    local bases=$1 && shift
    cat >${ref} <<REF
>ref
${bases}
REF
    samtools faidx ${ref}
}

samtools view -S -b - >${bam} <<SAM
@HD	VN:1.5	SO:coordinate
@SQ	SN:ref	LN:31
alt	0	ref	1	30	13=1X1=12X1=1X1=1X	*	0	0	CCCCCCCCCCCCAGTTAAAAAAAAAAAGGTT	AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
SAM
samtools index ${bam}

function run_freebayes() {
    freebayes --haplotype-length 0 --min-alternate-count 1 \
              --min-alternate-fraction 0 --pooled-continuous \
              --ploidy 1 \
              -f $ref $bam \
              2>&1 \
        | grep -vE "^#" | cut -f1-5
}


make_ref "AATCGGCTAZ"
expected='ERROR\(freebayes\): Found non-DNA character Z at position 9'
like "$(run_freebayes)" "${expected}" "freebayes rejects invalid reference"

expected=$(cat <<END
ref	14	.	A	G
ref	16	.	CNNNNNNNNNNN	TAAAAAAAAAAA
ref	29	.	C	G
ref	31	.	A	T
END
)
make_ref "CCCCCCCCCCCCAATCURYKMSWBDHVGCTA"
is "$(run_freebayes)" "${expected}" "freebayes does not put uppercase IUPAC bases in VCFs"

make_ref "CCCCCCCCCCCCaatcurykmswbdhvgcta"
is "$(run_freebayes)" "${expected}" "freebayes does not put lowercase IUPAC bases in VCFs"

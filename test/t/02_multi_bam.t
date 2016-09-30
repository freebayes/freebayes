#!/usr/bin/env bash

BASH_TAP_ROOT=bash-tap
source ./bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for freebayes

plan tests 7

ref=$(basename $0).ref

trap 'rm -f ${ref}* $(basename $0)*.bam*' EXIT

cat >${ref} <<REF
>ref
AATCGGCTA
REF
samtools faidx ${ref}

function run_freebayes() {
    freebayes "$@" \
              --haplotype-length 0 --min-alternate-count 1 \
              --min-alternate-fraction 0 --pooled-continuous --report-monomorphic \
              --ploidy 1 \
              -f $ref $bam \
              2>&1 \
        | grep -vE "^#" | cut -f1-5
}

function make_bam() {
    local id=$1 && shift
    local sm=$1 && shift
    local pl=$1 && shift
    local suffix=${1:-} && shift
    local first_snp=${1:-G} && shift

    local bam="$(basename $0).${id}.${sm}.${pl}${suffix}.bam"
    samtools view -S -b - >${bam} <<SAM
@HD	VN:1.5	SO:coordinate
@SQ	SN:ref	LN:9
@RG	ID:${id}	SM:${sm}	PL:${pl}
alt	0	ref	1	30	1=1X1=2X1=1X1=1X	*	0	0	A${first_snp}TTAGGTT	*	RG:Z:${id}
SAM
    samtools index ${bam}
    echo ${bam}
}

expected=$(cat <<END
ref	2	.	A	G
ref	4	.	CG	TA
ref	7	.	C	G
ref	9	.	A	T
END
)
bam1=$(make_bam "id1" "sample1" "platform1" "a")
bam2=$(make_bam "id1" "sample1" "platform1" "b")
is "$(run_freebayes -f ${ref} ${bam1} ${bam2})" "${expected}" "freebayes calls from two identical BAMs"

bam1=$(make_bam "id1" "sample1" "platform1")
bam2=$(make_bam "id2" "sample2" "platform2")
is "$(run_freebayes -f ${ref} ${bam1} ${bam2})" "${expected}" "freebayes calls from two BAMs with different samples for different read groups"

bam1=$(make_bam "id1" "sample1" "platform1")
bam2=$(make_bam "id2" "sample1" "platform2")
is "$(run_freebayes -f ${ref} ${bam1} ${bam2})" "${expected}" "freebayes calls from two BAMs with different technologies for different read groups"

bam1=$(make_bam "id1" "sample1" "platform1")
bam2=$(make_bam "id1" "sample2" "platform1")
expected='ERROR\(freebayes\): multiple samples \(SM\) map to the same read group \(RG\)'
like "$(run_freebayes -f ${ref} ${bam1} ${bam2})" "${expected}" "freebayes rejects two BAMs with different samples for same read groups"

bam1=$(make_bam "id1" "sample1" "platform1")
bam2=$(make_bam "id1" "sample1" "platform2")
expected='ERROR\(freebayes\): multiple technologies \(PL\) map to the same read group \(RG\)'
like "$(run_freebayes -f ${ref} ${bam1} ${bam2})" "${expected}" "freebayes rejects two BAMs with different technologies for same read groups"

bam1=$(make_bam "id1" "sample1" "platform1")
bam2=$(make_bam "id2" "sample2" "platform2" "" "C")
expected=$(cat <<END
ref	2	.	A	C,G
ref	4	.	CG	TA
ref	7	.	C	G
ref	9	.	A	T
END
)
is "$(run_freebayes -f ${ref} ${bam1} ${bam2})" "${expected}" "freebayes calls multiple alts from two different BAMs"

is "$(run_freebayes -f ${ref} ${bam1} ${bam1})" "Error: Duplicate bam file '02_multi_bam.t.id1.sample1.platform1.bam'" "freebayes rejects two BAMs with the same name"

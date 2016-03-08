#!/usr/bin/env bash
#
# Test that verifies a bug where a high-frequency variant is being dropped
# when the --only-use-input-alleles argument is given.

BASH_TAP_ROOT=bash-tap
. ./bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH

plan tests 5

result_line=$(freebayes -f tiny/chr12_fragment.fa --variant-input tiny/02_only-use-input-alleles_missed_call_target.vcf --only-use-input-alleles tiny/02_only-use-input-alleles_missed_call.bam |grep '^chr12_fragment')

is $(echo "$result_line" |awk -F $'\t' 'BEGIN {OFS = FS} {print $2}') 186 "starting coord" 
is $(echo "$result_line" |awk -F $'\t' 'BEGIN {OFS = FS} {print $4}') CA "reference"
is $(echo "$result_line" |awk -F $'\t' 'BEGIN {OFS = FS} {print $5}') "CAGCTCCA,AA,A" "alts"

# Make an assertion that the floating-point QUAL column is greater than a
# certain threshold... this is a little awkward in bash since it doesn't
# natively support floating-point comparisons.
qual=$(echo "$result_line" |awk -F $'\t' 'BEGIN {OFS = FS} {print $6}')
is $(bc <<< "$qual > 10") 1 "quality sufficiently high"

is $(echo "$result_line" |awk -F $'\t' 'BEGIN {OFS = FS} {print $8}' |grep -P -o '(^|;)AF=(\d|,)+' |sed 's/;//') "AF=0,1,0" "allele frequency"

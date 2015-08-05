#!/bin/bash

if [ $# -lt 3 ];
then
    echo "usage: $0 [regions file] [ncpus] [freebayes arguments]"
    echo
    echo "Run freebayes in parallel over regions listed in regions file, using ncpus processors."
    echo "Will merge and sort output, producing a uniform VCF stream on stdout.  Flags to freebayes"
    echo "which would write to e.g. a particular file will obviously cause problms, so caution is"
    echo "encouraged when using this script."
    echo
    echo "examples:"
    echo
    echo "Run freebayes in parallel on 100000bp chunks of the ref (fasta_generate_regions.py is also"
    echo "located in the scripts/ directory in the freebayes distribution).  Use 36 threads."
    echo
    echo "    freebayes-parallel <(fasta_generate_regions.py ref.fa.fai 100000) 36 -f ref.fa aln.bam >out.vcf"
    echo
    echo "Generate regions that are equal in terms of data content, and thus have lower variance"
    echo "in runtime.  This will yield better resource utilization."
    echo
    echo "    bamtools coverage -in aln.bam | coverage_to_regions.py ref.fa 500 >ref.fa.500.regions"
    echo "    freebayes-parallel ref.fa.500.regions 36 -f ref.fa aln.bam >out.vcf"
    echo
    exit
fi

regionsfile=$1
shift
ncpus=$1
shift

command=("freebayes" "$@")

(
#$command | head -100 | grep "^#" # generate header
# iterate over regions using gnu parallel to dispatch jobs
cat "$regionsfile" | parallel -k -j "$ncpus" "${command[@]}" --region {}
) | vcffirstheader \
    | vcfstreamsort -w 1000 | vcfuniq # remove duplicates at region edges

#!/usr/bin/python

import sys


if len(sys.argv) == 1:
    print "usage: ", sys.argv[0], " <fasta index file> <region size>"
    print "generates a list of freebayes/bamtools region specifiers on stdout"
    print "intended for use in creating cluster jobs"
    exit(1)

fasta_index_file = open(sys.argv[1])
region_size = int(sys.argv[2])

for line in fasta_index_file:

    fields = line.strip().split("\t")

    chrom_name = fields[0]
    chrom_length = int(fields[1])
    region_start = 0

    while region_start < chrom_length:
        start = region_start
        end = region_start + region_size
        if end > chrom_length:
            end = chrom_length
        print chrom_name + ":" + str(region_start) + "-" + str(end)
        region_start = end



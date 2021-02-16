#!/usr/bin/env python3

import argparse
from math import ceil


def generate_regions(fasta_index_file, size, chunks=False, chromosomes=None, bed_files=None):

    if not fasta_index_file.endswith(".fai"):
        fasta_index_file = fasta_index_file + ".fai"

    with open(fasta_index_file, "r") as fasta_index:
        for line in fasta_index:
            fields = line.strip().split("\t")
            chrom_name = fields[0]
            chrom_length = int(fields[1])
            if chromosomes is not None and chrom_name not in chromosomes:
                continue
            region_start = 0
            if chunks is True:
                region_size = ceil(chrom_length / size)  # have to make sure this works
            else:
                region_size = size
            while region_start < chrom_length:
                region_end = region_start + region_size
                if region_end > chrom_length:
                    region_end = chrom_length
                start = str(region_start)
                end = str(region_end)
                if bed_files is not None:
                    region = str(ceil(region_end / region_size))
                    file_path = f"{bed_files}.{chrom_name}.region.{region}.bed"
                    with open(file_path, "w") as f:
                        f.write("\t".join([chrom_name, start, end]))
                else:
                    print(f"{chrom_name}:{start}-{end}")
                region_start = region_end


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Generates a list of freebayes/bamtools region specifiers. Intended "
                                                 "for parallelization or creating cluster jobs.")

    parser.add_argument("--chunks", action="store_true",
                        help="Split the fasta into N chunks rather than into N length pieces")
    parser.add_argument("--chromosomes", nargs="+", default=None,
                        help="List of chromosomes to create chunks for")
    parser.add_argument("--bed", metavar="base name", type=str,
                        help="Write chunks to individual bed files (for Snakemake) instead of stdout.")
    parser.add_argument("fai", metavar="<fasta or fai file>",
                        help="The fasta file to split. Must be indexed.")
    parser.add_argument("region_size", metavar="<N>", type=int,
                        help="Region size")

    args = parser.parse_args()
    generate_regions(args.fai, args.region_size, chunks=args.chunks, chromosomes=args.chromosomes, bed_files=args.bed)

#!/usr/bin/env python3

import argparse
from math import ceil
import sys

def parse_chromosomes_arg(chrom_args):
    """
    Parses region arguments (e.g., '1', '2:10000-', '3:500-2000').
    Supports comma or space separation.
    Returns: dictionary {chrom: [(start, end), ...]}
    """
    if not chrom_args:
        return None

    targets = {}
    
    # Flatten list, handling comma-separated input (e.g., '1,2:100-')
    raw_list = []
    for item in chrom_args:
        raw_list.extend(item.split(','))

    for item in raw_list:
        if not item: continue
        
        chrom = item
        start = 0
        end = None # None means end of chromosome

        # Parse 'chrom:coords' format
        if ':' in item:
            chrom, coords = item.split(':')
            if '-' in coords:
                start_s, end_s = coords.split('-')
                if start_s:
                    start = int(start_s)
                if end_s:
                    end = int(end_s)
            else:
                # If only start position is given (e.g., 2:1000)
                start = int(coords)
        
        if chrom not in targets:
            targets[chrom] = []
        targets[chrom].append((start, end))
    
    return targets

def generate_regions(fasta_index_file, size, chunks=False, target_regions=None, bed_files=None):

    if not fasta_index_file.endswith(".fai"):
        fasta_index_file = fasta_index_file + ".fai"

    with open(fasta_index_file, "r") as fasta_index:
        for line in fasta_index:
            fields = line.strip().split("\t")
            chrom_name = fields[0]
            chrom_length = int(fields[1])

            # Filter: If specific regions are given, skip chromosomes not requested
            if target_regions is not None and chrom_name not in target_regions:
                continue

            # Set regions to process: custom regions or the entire chromosome (0, None)
            regions_to_process = target_regions[chrom_name] if target_regions else [(0, None)]

            for (req_start, req_end) in regions_to_process:
                
                # Determine actual bounds based on request and chromosome length
                actual_start = req_start
                actual_end = req_end if req_end is not None else chrom_length

                # Safety cap: ensure end doesn't exceed chromosome length
                if actual_end > chrom_length:
                    actual_end = chrom_length
                if actual_start >= actual_end:
                    continue # Skip invalid or empty ranges

                target_length = actual_end - actual_start

                # Recalculate region size if chunk mode is active
                if chunks is True:
                    # Split the target sub-region into N chunks
                    region_size = ceil(target_length / size) 
                else:
                    region_size = size

                # Start generating chunks within the defined bounds
                current_pos = actual_start
                while current_pos < actual_end:
                    step_end = current_pos + region_size
                    if step_end > actual_end:
                        step_end = actual_end

                    start_str = str(current_pos)
                    end_str = str(step_end)

                    if bed_files is not None:
                        # Write to BED file logic
                        chunk_idx = str(ceil((step_end - actual_start) / region_size))
                        file_path = f"{bed_files}.{chrom_name}.region.{chunk_idx}.bed"
                        with open(file_path, "w") as f:
                            f.write("\t".join([chrom_name, start_str, end_str]) + "\n")
                    else:
                        print(f"{chrom_name}:{start_str}-{end_str}")

                    current_pos = step_end


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Generates FreeBayes/Bamtools region specifiers "
                                                 "from a FASTA index for parallelization.")

    parser.add_argument("--chunks", action="store_true",
                        help="Split the region into N chunks rather than into fixed length pieces.")
    
    parser.add_argument("--chromosomes", nargs="+", default=None,
                        help="List of chromosomes or specific regions to process. "
                             "Format: '1' or '2:1000-' or '3:100-200'. Can be space or comma separated.")
    
    parser.add_argument("--bed", metavar="base name", type=str,
                        help="Write chunks to individual bed files instead of stdout.")
    parser.add_argument("fai", metavar="<fasta or fai file>",
                        help="The FASTA file to split. Must be indexed.")
    parser.add_argument("region_size", metavar="<N>", type=int,
                        help="Region size (bp) or Number of chunks (if --chunks is used).")

    args = parser.parse_args()

    # Parse arguments before calling generator
    parsed_targets = parse_chromosomes_arg(args.chromosomes)

    generate_regions(args.fai, args.region_size, chunks=args.chunks, target_regions=parsed_targets, bed_files=args.bed)

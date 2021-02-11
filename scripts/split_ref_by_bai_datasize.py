#!/usr/bin/env python3
"""
Create a list of regions by splitting a reference based on the amount of data in bam files.
Uses the `bai` index of the bam files. Useful for submitting jobs of equal size to a cluster.
"""

import sys
import os
import argparse
import time
import logging

import struct
import numpy as np
from scipy import interpolate


DEFAULT_LOGGING_LEVEL = logging.INFO
MAX_LOGGING_LEVEL = logging.CRITICAL

def setup_logger(verbose_level):
    fmt=('%(levelname)s %(asctime)s [%(module)s:%(lineno)s %(funcName)s] :: '
            '%(message)s')
    logging.basicConfig(format=fmt, level=max((0, min((MAX_LOGGING_LEVEL,
                        DEFAULT_LOGGING_LEVEL-(verbose_level*10))))))


def Main(argv):
    tic_total = time.time()

    # parse arguments
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('bamfiles', metavar='BAMFILE', nargs='*')
    parser.add_argument('-L', '--bam-list', nargs='*')
    parser.add_argument('-r', '--reference-fai', help="reference fasta index file", required=True)
    parser.add_argument('-s', '--target-data-size', default='100e6', help="target combined data size of bam files in each region (MB)")
    parser.add_argument('--bai-interval-size', default=16384, type=int, help="Size in baseparis of each interval in the bam index (bai).")
    parser.add_argument('-v', '--verbose', action='count', default=0,
                        help="increase logging verbosity")
    parser.add_argument('-q', '--quiet', action='count', default=0,
                        help="decrease logging verbosity")
    args = parser.parse_args(argv)

    # setup logger
    setup_logger(verbose_level=args.verbose-args.quiet)
    if argv is not None:
        logging.warning('Using passed arguments: '+str(argv))
    logging.info('args: '+str(args))

    # additional argument parsing and datatype handling
    if not args.bamfiles and not args.bam_list:
        logging.error("Must provide an BAMFILE and/or --bam-list argument")
        sys.exit(2)
    args.target_data_size = int(float(args.target_data_size))*1000000
    logging.info('target-data-size: '+str(args.target_data_size)+' bytes')

    # read bam-lists if provided
    if args.bam_list:
        for bamlistfile in args.bam_list:
            with open(bamlistfile,'r') as fh:
                for x in fh:
                    x = x.split('#')[0].strip()
                    if x:
                        args.bamfiles.append(x)
    #logging.info('bam files: '+", ".join(args.bamfiles)) # output complete list of bam files being used

    # read the reference fasta index
    fai_chrom = []
    fai_len = []
    with open(args.reference_fai,'r') as fh:
        for x in fh:
            x = x.strip().split(sep='\t')
            fai_chrom.append(str(x[0]))
            fai_len.append(int(x[1]))


    ## read bai indexes, skipping bin info
    # list by chrom of number of intervals
    n_intvs = np.array([int(np.ceil(x/args.bai_interval_size)) for x in fai_len])
    # list by chrom of lists of interval offsets
    icumsz = [] # cumulative size of data by interval
    for i,n in enumerate(n_intvs):
        icumsz.append(np.zeros((n,), dtype=np.int64))

    for bamfn in args.bamfiles:
        baifn = bamfn+'.bai'
        with open(baifn,'rb') as fh:
            logging.info("processing: "+baifn)
            # filetype magic check
            assert struct.unpack('4s', fh.read(4))[0] == b'BAI\x01'

            # number of reference sequences (chroms)
            n_ref = struct.unpack('i', fh.read(4))[0]
            assert n_ref == len(fai_len), "fasta index and bam index have must have same number of chroms"

            for ci in range(n_ref):
                # skip over the binning index
                n_bin = struct.unpack('i', fh.read(4))[0]
                for bini in range(n_bin):
                    bin_id = struct.unpack('I', fh.read(4))[0]
                    n_chunk = struct.unpack('i', fh.read(4))[0]
                    fh.seek(n_chunk*16, os.SEEK_CUR)
                # read interval index
                n_intv = struct.unpack('i', fh.read(4))[0]
                if n_intv > 0:
                    ioff = np.array(struct.unpack(str(n_intv)+'Q', fh.read(n_intv*8)), dtype=np.int64)
                    while( len(ioff) < len(icumsz[ci]) ):
                        ioff = np.append(ioff, ioff[-1]+1)
                    icumsz[ci] += ioff-ioff[0]



    ## make the list of regions
    regions = []

    for ci,chrom in enumerate(fai_chrom):

        # sanity check last point if there are more than one
        if len(icumsz[ci]) > 1:
            assert icumsz[ci][-1] >= icumsz[ci][-2]

        # tiny chroms just get 1 region
        if len(icumsz[ci]) < 2:
            regions.extend([ (fai_chrom[ci], 0, fai_len[ci]) ])
            continue
        ds = icumsz[ci]
        pos = np.arange(0, ds.shape[0])*args.bai_interval_size

        # estimate total data size for the chrom
        f = interpolate.interp1d(pos, ds, fill_value='extrapolate', kind='linear')
        ds_total = f([fai_len[ci]])[0]

        num_regions = int(np.ceil(ds_total/args.target_data_size))

        # approx equal LENGTH regions
        # tmp = np.linspace(0, fai_len[ci], num=num_regions+1, endpoint=True, dtype=int)

        # approx equal DATA SIZE regions
        f = interpolate.interp1d(ds, pos, fill_value='extrapolate', kind='linear')
        dsx = np.linspace(0, ds_total, num=num_regions+1, endpoint=True, dtype=int)
        tmp = f(dsx).astype(int)
        # ensure we exactly hit the endpoints
        tmp[0] = 0
        tmp[-1] = fai_len[ci]

        regions.extend([ (fai_chrom[ci], tmp[i], tmp[i+1]) for i in range(len(tmp)-1) ])

    ## Output regions file
    for r in regions:
        print(*r, sep='\t')

    logging.info("Number of chroms: {}".format(len(fai_len)))
    logging.info("Number of splits: {}".format(len(regions)-len(fai_len)))
    logging.info("Number of regions: {}".format(len(regions)))

    logging.info("Done: {:.2f} sec elapsed".format(time.time()-tic_total))
    return 0



#########################################################################
# Main loop hook... if run as script run main, else this is just a module
if __name__ == '__main__':
    if 'TESTING_ARGS' in globals():
        sys.exit(Main(argv=TESTING_ARGS))
    else:
        sys.exit(Main(argv=None))


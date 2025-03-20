#!/usr/bin/env python
import argparse
import random
from tqdm import tqdm as tqdm
import numpy as np

def subsample_kraken2(input_path, bs_id_str, mhg, cs, subsampling_fractions, replicate_nr, output_path):
    with open('{}BS18-{}_mhg{}_cs{}.kraken2'.format(input_path, bs_id_str, mhg, str(cs).replace('.', '')), 'r') as inf:
        lines = inf.readlines()

    for sf in tqdm(subsampling_fractions):
        line_nr = int(sf * len(lines))
        for r in range(1, replicate_nr + 1):
            with open('{}BS18-{}_mhg{}_cs{}_sf{}_rep{}.kraken2'.format(output_path, bs_id_str, mhg, str(cs).replace('.', ''), str(sf).replace('.', ''), r), 'w') as outf:
                for l in random.sample(lines, line_nr):
                    outf.write(l)

def main():
    #Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--kraken_path')
    parser.add_argument('--bs_id')
    parser.add_argument('--subsample_path')
    args = parser.parse_args()

    # padding 
    bs_id_str = '%04d' % args.bs_id

    # hard coded parameters
    # mhg_cs = [(2, '00'), (4, '005')]
    mhg_cs = [(2, '00')]
    subsampling_fractions = [round(x, 2) for x in np.arange(0.1, 1, 0.1)]
    replicate_nr = 10
    
    for (mhg, cs) in mhg_cs:
        subsample_kraken2(args.kraken_path, bs_id_str, mhg, cs, subsampling_fractions, replicate_nr, args.subsample_path)

if __name__ == "__main__":
    main()

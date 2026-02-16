#!/usr/bin/env python
import argparse
import random
import subprocess
import numpy as np
from itertools import product

import sys
sys.path.append('/cluster/raid/home/f80878961/scripts/plan-b-omics/notebooks')

from bee_metagenomics_lib import get_kraken2_read_nr

#def subsample_kraken2(kraken2_fn, subsample_kraken2, sf, r):
#    with open(kraken2_fn, 'r') as inf:
#        lines = inf.readlines()
#    line_nr = int(sf * len(lines))
#    with open(subsample_kraken2, 'w') as outf:
#        for l in random.sample(lines, line_nr):
#            outf.write(l)

def subsample_kraken2(kraken2_fn, subsample_kraken2, sf, r, max_read_nr):
    with open(kraken2_fn, 'r') as inf:
        lines = inf.readlines()
    
    # subsample a fraction of the largest library
    sampled_read_nr = int(sf * max_read_nr)
    
    with open(subsample_kraken2, 'w') as outf:
        for l in random.sample(lines, min(sampled_read_nr, len(lines))):
            outf.write(l)
            
def write_k2report(make_kreport_script, subsample_kraken2, ktaxonomy, subsample_k2report, kraken2_path):
    command = ['python', make_kreport_script, '-i', subsample_kraken2, '-t', ktaxonomy, '-o', subsample_k2report]
    result = subprocess.run(command, cwd=kraken2_path, capture_output=False)

def run_bracken(kdb_path, krakdb, subsample_k2report, level, subsample_bracken, subsample_breport, bracken_path):
    command = ['bracken', '-d', '{}{}'.format(kdb_path, krakdb), '-i', subsample_k2report, '-r', '150', '-l', level, '-t', '10', '-o', subsample_bracken, '-w', subsample_breport]
    result = subprocess.run(command, cwd=bracken_path, capture_output=False)

def rarefy_colony(cs_thresholds, subsampling_fractions, replicate_nr, 
    kraken2_path, bs_id, krakdb, readpool, mhg, 
    kdb_path, max_read_nr, make_kreport_script, levels, bracken_path):
    
    for cs, sf, r in product(cs_thresholds, subsampling_fractions, range(1, replicate_nr + 1)):
        print(cs, sf, r)
    
        kraken2_fn = '{}{}_{}_{}_mhg{}_cs{}_sf1_rep1.kraken2'.format(kraken2_path, bs_id, krakdb, readpool, mhg, cs)
        subsample_name = '{}_{}_{}_mhg{}_cs{}_sf{}_rep{}'.format(bs_id, krakdb, readpool, mhg, cs, str(sf).replace('.', ''), r)
        subsample_kraken2_fn = '{}{}.kraken2'.format(kraken2_path, subsample_name)
        subsample_k2report_fn = '{}{}.k2report'.format(kraken2_path, subsample_name)
        ktaxonomy = '{}{}/mydb_taxonomy.txt'.format(kdb_path, krakdb)

        # break when fraction * max_read_nr > library size
        if int(sf * max_read_nr) > get_kraken2_read_nr(kraken2_path, bs_id, krakdb, readpool, mhg, cs, 1, 1):
            break
        
        # subsample kraken2
        subsample_kraken2(kraken2_fn, subsample_kraken2_fn, sf, r, max_read_nr)
        
        # write kreport
        write_k2report(make_kreport_script, subsample_kraken2_fn, ktaxonomy, subsample_k2report_fn, kraken2_path)
        
        for l in levels:
            subsample_bracken_fn = '{}{}_{}.bracken'.format(bracken_path, subsample_name, l)
            subsample_breport_fn = '{}{}_{}.breport'.format(bracken_path, subsample_name, l)
            run_bracken(kdb_path, krakdb, subsample_k2report_fn, l, subsample_bracken_fn, subsample_breport_fn, bracken_path)

#def subsample_kraken2(input_path, bs_id_str, krakdb, readpool, mhg, cs, subsampling_fractions, replicate_nr, output_path):
#    with open('{}BS18-{}_{}_{}_mhg{}_cs{}.kraken2'.format(input_path, bs_id_str, krakdb, readpool, mhg, cs), 'r') as inf:
#        lines = inf.readlines()
#
#    for sf in tqdm(subsampling_fractions):
#        line_nr = int(sf * len(lines))
#        for r in range(1, replicate_nr + 1):
#            with open('{}BS18-{}_{}_{}_mhg{}_cs{}_sf{}_rep{}.kraken2'.format(output_path, bs_id_str, krakdb, readpool, mhg, cs, str(sf).replace('.', ''), r), 'w') as outf:
#                for l in random.sample(lines, line_nr):
#                    outf.write(l)

def main():
    # hard coded parameters
    krakdb, readpool, mhg = ('corent', 'nonbee', 2)
    levels = ['S', 'G', 'F']
    # cs_thresholds = ['005', '025', '050']
    cs_thresholds = ['050']
    subsampling_fractions = [round(x, 2) for x in np.arange(0.1, 1, 0.1)]
    replicate_nr = 1
    
    #Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--kraken2_path')
    parser.add_argument('--bs_id', type=str)
    parser.add_argument('--krakdb', type=str)
    parser.add_argument('--readpool', type=str)
    parser.add_argument('--mhg', type=int)
    parser.add_argument('--kdb_path', type=str)
    parser.add_argument('--max_read_nr', type=int)
    parser.add_argument('--make_kreport_script', type=str)
    parser.add_argument('--bracken_path', type=str)
    args = parser.parse_args()
    
    rarefy_colony(cs_thresholds, subsampling_fractions, replicate_nr,
        args.kraken2_path, args.bs_id, args.krakdb, args.readpool, args.mhg, 
        args.kdb_path, args.max_read_nr, args.make_kreport_script, levels, args.bracken_path)

if __name__ == "__main__":
    main()

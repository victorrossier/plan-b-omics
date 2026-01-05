#!/usr/bin/env python
import numpy as np
from tqdm import tqdm as tqdm
import pandas as pd
import bz2
import argparse

def write_allele_table(allele_count_mat, snp_pos, chr_ids, out_path, allele_name, batch_id):
    allele_count_df = pd.DataFrame(allele_count_mat)
    allele_count_df.insert(0, 'SNP_position', snp_pos)
    allele_count_df.insert(0, 'Chromosome', chr_ids)
    allele_count_df.to_csv('{}count_{}_{}.txt'.format(out_path, allele_name, batch_id), header=False, index=False, sep='\t')

def merge_sync_batch(
    batch_id,
    snp_nr,
    bs_ids,
    count_path,
    out_path):
    '''
    merge sync files from popoolation (1 per colony) to A, T, C, G counts --> then processed together for allele filtering
    '''
    A_count = np.zeros((snp_nr, len(bs_ids)), dtype=np.uint16)
    T_count = np.zeros((snp_nr, len(bs_ids)), dtype=np.uint16)
    C_count = np.zeros((snp_nr, len(bs_ids)), dtype=np.uint16)
    G_count = np.zeros((snp_nr, len(bs_ids)), dtype=np.uint16)
    N_count = np.zeros((snp_nr, len(bs_ids)), dtype=np.uint16)
    del_count = np.zeros((snp_nr, len(bs_ids)), dtype=np.uint16)
    max_int = np.iinfo(np.uint16).max
    chr_ids = [""] * snp_nr
    snp_pos = [""] * snp_nr

    j = 0
    for bsid in tqdm(bs_ids):
        with bz2.open('{}{}.count.bz2'.format(count_path, bsid), 'rt') as inf:
            inf.readline()
            i = 0
            for l in inf:
                ch, pos, ref, alt, counts = l.split()
                atcgnd = np.array(list(map(int, counts.split(':'))))
    
                # fill positional info first
                if j == 0:
                    chr_ids[i] = ch
                    snp_pos[i] = pos
                # check consistent chromosome and position
                else:
                    assert ch == chr_ids[i], 'inconsistent chromosome id'
                    assert pos == snp_pos[i], 'inconsistent position'
            
                if max(atcgnd) > max_int:
                    print(bs_id, ch, pos, atcgnd)
                    
                A_count[i, j] = atcgnd[0]
                T_count[i, j] = atcgnd[1]
                C_count[i, j] = atcgnd[2]
                G_count[i, j] = atcgnd[3]
                N_count[i, j] = atcgnd[4]
                del_count[i, j] = atcgnd[5]     
                i += 1
                
            assert i == len(A_count) ## check that snp_nr was correctly inputed
        j += 1

    write_allele_table(A_count, snp_pos, chr_ids, out_path, 'A', batch_id)
    write_allele_table(T_count, snp_pos, chr_ids, out_path, 'T', batch_id)
    write_allele_table(C_count, snp_pos, chr_ids, out_path, 'C', batch_id)
    write_allele_table(G_count, snp_pos, chr_ids, out_path, 'G', batch_id)
    write_allele_table(N_count, snp_pos, chr_ids, out_path, 'N', batch_id)
    write_allele_table(del_count, snp_pos, chr_ids, out_path, 'del', batch_id)

def main():
    #Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--batch_id')
    parser.add_argument('--snp_nr', type=int)
    parser.add_argument('--bs_ids', type=str)
    parser.add_argument('--count_path', type=str)
    parser.add_argument('--out_path', type=str) 
    args = parser.parse_args()

    bs_ids = args.bs_ids.split(',')
    
    merge_sync_batch(args.batch_id, args.snp_nr, bs_ids, args.count_path, args.out_path) 

if __name__ == "__main__":
    main()

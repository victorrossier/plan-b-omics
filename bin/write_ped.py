#!/usr/bin/env python
import argparse

import numpy as np
from tqdm import tqdm as tqdm
import pandas as pd

# def write_ped(bim_file, chunk_bs_ids, snp_nr, bgs_file, batch_j, chunk_size, plink_path, prefix):
# 
#     # contains REF and ALT alleles
#     bim_df = pd.read_csv(bim_file, sep='\t')
#     
#     ped_mat = np.full((len(chunk_bs_ids), 2 * snp_nr), '', dtype='U1')
#     
#     with open(bgs_file, 'r') as inf:
#         
#         # index of the first colony of the batch/chunk
#         start_j = batch_j * chunk_size
#     
#         # check we have the same bs ids in header than in batch
#         header_bsids = inf.readline()
#         assert header_bsids.rstrip().split(',')[3:][start_j: start_j + len(chunk_bs_ids)] == chunk_bs_ids, 'bs_id mismatch'
#         
#         for i in tqdm(range(snp_nr)):
#             # parse number of ref allele per SNP (0, 1, 2 as queen is diploid)
#             ref_numbers = inf.readline().rstrip().split(',')[3:]
#     
#             for j in range(start_j, start_j + len(chunk_bs_ids)):
#                 ref_nr = ref_numbers[j]
#                 if ref_nr == '2':
#                     ped_mat[j - start_j, 2 * i] = bim_df.iloc[i, 5]
#                     ped_mat[j - start_j, 2 * i + 1] = bim_df.iloc[i, 5]
#                 elif ref_nr == '1': 
#                     ped_mat[j - start_j, 2 * i] = bim_df.iloc[i, 5]
#                     ped_mat[j - start_j, 2 * i + 1] = bim_df.iloc[i, 4]
#                 elif ref_nr == '0': 
#                     ped_mat[j - start_j, 2 * i] = bim_df.iloc[i, 4]
#                     ped_mat[j - start_j, 2 * i + 1] = bim_df.iloc[i, 4]
#                 else: 
#                     print('ERROR')
# 
#     ped_df = pd.DataFrame({
#         'FID': chunk_bs_ids, 
#         'IID': chunk_bs_ids,
#         'PID': ['0'] * len(chunk_bs_ids),
#         'MID': ['0'] * len(chunk_bs_ids),
#         'SEX': ['0'] * len(chunk_bs_ids), 
#         'PHENOTYPE': ['-9'] * len(chunk_bs_ids)
#     })
#     ped_df = ped_df.join(pd.DataFrame(ped_mat))
#     ped_df.to_csv('{}{}_{}.ped'.format(plink_path, prefix, f'{batch_j:03}'), header=False, index=False, sep='\t')

def write_ped_pop(
    plink_path, prefix, pop, chunk_bs_ids, batch_j, chunk_size, prob_threshold):

    bim_like = pd.read_csv('{}{}.bim_like'.format(plink_path, prefix), sep='\t', header=None)
    bgs_file = '{}{}_{}.bgs'.format(plink_path, prefix, pop)
    aa_prob_file = '{}{}_{}_AA.prob'.format(plink_path, prefix, pop)
    ar_prob_file = '{}{}_{}_AR.prob'.format(plink_path, prefix, pop)
    rr_prob_file = '{}{}_{}_RR.prob'.format(plink_path, prefix, pop)
    snp_nr = len(bim_like)
    
    start_j = batch_j * chunk_size
    ped_mat = np.full((len(chunk_bs_ids), 2 * snp_nr), '', dtype='U1')
    
    inf = open(bgs_file, 'r')
    inf_aa = open(aa_prob_file, 'r')
    inf_ar = open(ar_prob_file, 'r')
    inf_rr = open(rr_prob_file, 'r')

    # check we have the same bs ids in header than in batch
    header_bsids = inf.readline()
    assert header_bsids.rstrip().split(',')[3:][start_j: start_j + len(chunk_bs_ids)] == chunk_bs_ids, 'bs_id mismatch'
    inf_aa.readline()
    inf_ar.readline()
    inf_rr.readline()

    # iterate snps and colonies
    for i in tqdm(range(snp_nr)):
        ref_numbers = inf.readline().rstrip().split(',')[3:]
        aa_probs = inf_aa.readline().rstrip().split(',')[3:]
        ar_probs = inf_ar.readline().rstrip().split(',')[3:]
        rr_probs = inf_rr.readline().rstrip().split(',')[3:]

        for j in range(start_j, start_j + len(chunk_bs_ids)):
            ref_nr = ref_numbers[j]

            # homozygote reference
            if (ref_nr == '2') and (float(rr_probs[j]) >= prob_threshold):
                ped_mat[j - start_j, 2 * i] = bim_like.iloc[i, 5]
                ped_mat[j - start_j, 2 * i + 1] = bim_like.iloc[i, 5]
            # heterozygote
            elif (ref_nr == '1') and (float(ar_probs[j]) >= prob_threshold):
                ped_mat[j - start_j, 2 * i] = bim_like.iloc[i, 5]
                ped_mat[j - start_j, 2 * i + 1] = bim_like.iloc[i, 4]
            # homozygote alternative
            elif (ref_nr == '0') and (float(aa_probs[j]) >= prob_threshold):
                ped_mat[j - start_j, 2 * i] = bim_like.iloc[i, 4]
                ped_mat[j - start_j, 2 * i + 1] = bim_like.iloc[i, 4]
            # below prob_threshold
            elif ref_nr in {'0', '1', '2'}:
                ped_mat[j - start_j, 2 * i] = '0'
                ped_mat[j - start_j, 2 * i + 1] = '0'
            else: 
                print('ERROR')
            
    inf.close()
    inf_aa.close()
    inf_ar.close()
    inf_rr.close()

    ped_df = pd.DataFrame({
        'FID': ['0'] * len(chunk_bs_ids), 
        'IID': chunk_bs_ids,
        'PID': ['0'] * len(chunk_bs_ids),
        'MID': ['0'] * len(chunk_bs_ids),
        'SEX': ['0'] * len(chunk_bs_ids), 
        'PHENOTYPE': ['-9'] * len(chunk_bs_ids)
    })
    ped_df = ped_df.join(pd.DataFrame(ped_mat))
    ped_df.to_csv('{}{}_{}_{}.ped'.format(plink_path, prefix, pop, batch_j), header=False, index=False, sep='\t')
    
def main():
    #Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--chunk_bs_ids', type=str)
    parser.add_argument('--batch_id', type=int)
    parser.add_argument('--chunk_size', type=int)
    parser.add_argument('--plink_path', type=str) 
    parser.add_argument('--prefix', type=str) 
    parser.add_argument('--pop', type=str) 
    parser.add_argument('--prob_threshold', type=float) 
    args = parser.parse_args()

    chunk_bs_ids = args.chunk_bs_ids.rstrip().split(',')
    batch_j = args.batch_id - 1
    
    write_ped_pop(args.plink_path, args.prefix, args.pop, chunk_bs_ids, batch_j, args.chunk_size, args.prob_threshold)
    # write_ped(args.bim_file, chunk_bs_ids, args.snp_nr, args.bgs_file, batch_j, args.chunk_size, args.plink_path, args.prefix)

if __name__ == "__main__":
    main()
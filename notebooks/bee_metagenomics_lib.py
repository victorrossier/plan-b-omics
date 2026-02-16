#!/usr/bin/env python
# coding: utf-8

import os
import glob
import gzip
import subprocess
import pandas as pd

################################################################################################################################################################
# 1. ete3 taxonomy helper function (if needed have stuff from PhD too...)
################################################################################################################################################################

def get_desc_taxa(node, cand_taxa):
    desc_taxa = []
    if node.sci_name in cand_taxa:
        desc_taxa.append(node.sci_name)
        return desc_taxa
    for ch in node.children:
        desc_taxa.extend(get_desc_taxa(ch, cand_taxa))
    return desc_taxa


################################################################################################################################################################
# 2. Slurm monitoring
################################################################################################################################################################

def write_bs_ids_file(bs_ids_file, bs_ids):
    with open(bs_ids_file, 'w') as outf:
        outf.write('\n'.join(bs_ids))

def read_bs_ids_file(bs_ids_file):
    with open(bs_ids_file, 'r') as inf:
        lines = inf.read().split('\n')
        array_id2bs_id = dict(zip(range(1, len(lines) + 1), lines))
    return array_id2bs_id, {v:k for k, v in array_id2bs_id.items()}

def write_array_str(array_ids_to_run, parallel_job_nr=100):
    '''
    convoluted way to build array string (as BSUB argument cannot be too long)
    '''
    if len(array_ids_to_run) == 1:
        return str(array_ids_to_run[0])
    array_str_list = []
    start = array_ids_to_run[0]
    for i in range(1, len(array_ids_to_run)):
        curr = array_ids_to_run[i]
        prev = array_ids_to_run[i - 1]
        if curr != prev + 1:
            if start == prev:
                array_str_list.append(str(prev))
            else:
                array_str_list.append('{}-{}'.format(start, prev))
            start = curr
    # end 
    if start == prev:
        array_str_list.append(str(prev))
    else:
        array_str_list.append('{}-{}'.format(start, prev + 1)) ## this + 1 is a hack
    start = curr
    
    return '{}%{}'.format(','.join(array_str_list), parallel_job_nr)

def get_status_job_array_ids(status='RUNNING', job_name='pp'):
    command = ['squeue', '-t', status, '-n', job_name, '-o', '%.18i']
    result = subprocess.run(command, capture_output=True, text=True)
    job_array_ids = []
    # to handle pending jobs...
    for job_str in result.stdout.split()[1:]:
        job_str = job_str.split('_')[1]
        if job_str.startswith('['):
            split_job_str = job_str[1:].split('%')[0].split(',')
            for js in split_job_str:
                if len(js.split('-')) > 1:
                    x = js.split('-')
                    job_array_ids.extend(list(range(int(x[0]), int(x[1]) + 1)))
                else:
                    job_array_ids.append(int(js))
        else:
            job_array_ids.append(int(job_str))
    return job_array_ids

def get_jobs_to_restart(step_name, completed_jobs_step_before, array_id2bs_id, step_path, expected_files, job_id="*"):

    # capture completed, running, pending, and failed/to_run jobs
    running_jobs = get_status_job_array_ids(status='RUNNING', job_name=step_name)
    pending_jobs = get_status_job_array_ids(status='PENDING', job_name=step_name)
    completed_jobs = []
    
    for array_id in completed_jobs_step_before:
        bs_id = array_id2bs_id[array_id]
        
        # ensure latest log file is there and finishes by DONE
        x = glob.glob('{}/{}_{}_{}.out'.format(step_path, step_name, job_id, array_id))
        if len(x) != 1:
            if len(x) > 1: 
                print(x)
            continue
    
        out_file = x[0]
        with open(out_file, 'r') as inf:
            lines = inf.readlines()
            no_error = True
            for l in lines:
                if l.startswith('(ERR)'):
                    no_error = False
            done = lines and lines[-1].startswith('DONE')

        # expected files (glob format) 
        expected_files_correct = True
        for exp_file in expected_files:
            fn = exp_file.replace('BSID', bs_id)
            if not os.path.exists(fn) or not os.path.getsize(fn) > 100:
                expected_files_correct = False
                break
        
        if done and no_error and expected_files_correct:
            completed_jobs.append(array_id)
            
    jobs_to_run = list(set(completed_jobs_step_before).difference(completed_jobs + running_jobs + pending_jobs))
    return jobs_to_run, completed_jobs

def clean_failed_jobs(failed_array_ids, array_id2bs_id, exp_files, step_path, step_name):
    for array_id in failed_array_ids:
        bs_id = array_id2bs_id[array_id]
        # remove expected output files
        for exp_file in exp_files:
            fn = exp_file.replace('BSID', bs_id)
            if os.path.exists(fn):
                os.remove(fn)
        # and slurm .out and .err files
        slurm_files = glob.glob('{}/{}_*_{}.*'.format(step_path, step_name, array_id))
        for f in slurm_files:
            if os.path.exists(f):
                os.remove(f)

################################################################################################################################################################
# 3. aligner parsers
################################################################################################################################################################

# should merge these two...
def parse_flagstat(flagstat_fn):
    """
    Parse a samtools flagstat output file.
    Returns number of reads for:
        - total
        - with_itself_and_mate_mapped
        - singletons
    """
    results = {
        "total": None,
        "with_itself_and_mate_mapped": None,
        "singletons": None
    }
    with open(flagstat_fn, "r") as inf:
        for line in inf:
            if "in total" in line:
                results["total"] = int(line.split()[0])
            
            elif "with itself and mate mapped" in line:
                results["with_itself_and_mate_mapped"] = int(line.split()[0])

            elif "singletons" in line:
                results["singletons"] = int(line.split()[0])
                
    return results
    
def parse_flagstats_ont(flagstats_fn):
    results = {
        'in total': None, 
        'primary': None, 
        'secondary': None, 
        'supplementary': None, 
        'duplicates': None, 
        'primary duplicates': None, 
        'mapped': None, 
        'primary mapped': None
    }
    with open(flagstats_fn, "r") as inf:
        results['in total'] = float(inf.readline().split()[0])
        results['primary'] = float(inf.readline().split()[0])
        results['secondary'] = float(inf.readline().split()[0])
        results['supplementary'] = float(inf.readline().split()[0])
        results['duplicates'] = float(inf.readline().split()[0])
        results['primary duplicates'] = float(inf.readline().split()[0])
        results['mapped'] = float(inf.readline().split()[0])
        results['primary mapped'] = float(inf.readline().split()[0])
    return results
    
def parse_nanostats(nanostats_fn):
    results = {
        'Mean read length': None, 
        'Mean read quality': None, 
        'Median read length': None, 
        'Median read quality': None, 
        'Number of reads': None, 
        'Read length N50': None, 
        'STDEV read length': None,
        'Total bases': None
    }
    with open(nanostats_fn, "r") as inf:
        inf.readline()
        results['Mean read length'] = float(inf.readline().split(':')[1].replace(',', '').strip())
        results['Mean read quality'] = float(inf.readline().split(':')[1].replace(',', '').strip())
        results['Median read length'] = float(inf.readline().split(':')[1].replace(',', '').strip())
        results['Median read quality'] = float(inf.readline().split(':')[1].replace(',', '').strip())
        results['Number of reads'] = float(inf.readline().split(':')[1].replace(',', '').strip())
        results['Read length N50'] = float(inf.readline().split(':')[1].replace(',', '').strip())
        results['STDEV read length'] = float(inf.readline().split(':')[1].replace(',', '').strip())
        results['Total bases'] = float(inf.readline().split(':')[1].replace(',', '').strip())

    return results

def parse_stats(stats_fn):
    
    results = {
        'sequences': None, # number of reads
        'reads mapped': None,
        'reads unmapped': None, 
        'total length': None, # number of bases total
        'bases mapped': None # number of bases from read that mapped
    }
    with open(stats_fn, "r") as inf:
        for l in inf:
            stat, value = l.split(':')
            if stat in results:
                results[stat] = float(value.split()[0].strip())

    return results

def count_lines_gz(gzip_fn):
    ''' faster'''
    with gzip.open(gzip_fn, "rb") as inf:
        return sum(buf.count(b"\n") for buf in iter(lambda: inf.read(1024*1024), b""))

def get_total_read_nr_flagstat(bowtie2_path, bs_id, idx_name):
    '''
    total number of paired reads
    '''
    return parse_flagstat('{}{}_{}_mapped.flagstat'.format(bowtie2_path, bs_id, idx_name))['total'] / 2

def get_unmapped_read_nr_flagstat(bowtie2_path, bs_id, idx_name):
    '''
    number of paired reads not mapping (i.e., to bee reference genome)
    '''
    flagstat_res = parse_flagstat('{}{}_{}_mapped.flagstat'.format(bowtie2_path, bs_id, idx_name))
    return (flagstat_res['total'] - flagstat_res['with_itself_and_mate_mapped'] - 2 * flagstat_res['singletons']) / 2

################################################################################################################################################################
# 3. Kraken2 and Bracken parsers and helprs
################################################################################################################################################################

def parse_kreport(file_path, bs_id, krakdb, readpool, mhg, cs, sf, r):
    kreport_fn = '{}{}_{}_{}_mhg{}_cs{}_sf{}_rep{}.k2report'.format(file_path, bs_id, krakdb, readpool, mhg, cs, str(sf).replace('.', ''), r)
    if not os.path.exists(kreport_fn): 
        print('missing {}'.format(kreport_fn))
    with open(kreport_fn, 'r') as inf:
        ucseqs_nr = int(inf.readline().split()[1])
        cseqs_nr = int(inf.readline().split()[1])
    return ucseqs_nr, cseqs_nr

def parse_breport(bracken_path, bs_id, krakdb, readpool, mhg, cs, sf, r, level, taxa):
    '''get read numbers for taxa'''
    breport_fn = '{}{}_{}_{}_mhg{}_cs{}_sf{}_rep{}_{}.breport'.format(bracken_path, bs_id, krakdb, readpool, mhg, cs, str(sf).replace('.', ''), r, level)
    if not os.path.exists(breport_fn): 
        print('missing {}'.format(breport_fn))
    breport_df = pd.read_csv(breport_fn, header=None, sep='\t')
    breport_df[5] = [x.lstrip() for x in breport_df[5].to_list()]
    return breport_df[breport_df[5].isin(taxa)].set_index(5)[1].to_dict()
    
def get_kraken2_read_nr(kraken2_path, bs_id, krakdb, readpool, mhg, cs, sf, r):
    '''
    number of sampled read nr = number of nonbee reads = number of reads that did not map to the bee reference genomes
    using kraken2 report instead of bowtie2 flagstats to consider subsampling (was done at kraken2 step)
    '''
    uc_k2, c_k2 = parse_kreport(kraken2_path, bs_id, krakdb, readpool, mhg, cs, sf, r)
    return uc_k2 + c_k2

def get_bracken_classified_read_nr(bracken_path, bs_id, krakdb, readpool, mhg, cs, sf, r, level):
    '''
    equals reads classified by Bracken at root and below
    '''
    return parse_breport(bracken_path, bs_id, krakdb, readpool, mhg, cs, sf, r, level, ['root'])['root']

def get_classified_read_nr(bowtie2_path, kraken2_path, bracken_path, bs_id, krakdb, readpool, mhg, cs, sf, r, level, idx_name):
    '''
    includes reads mapping to the bee genome normalized by the sampling proportion

    if sampling fraction = 1, then no subsampling was perform and thus (k2_sampled_read_nr / k2_total_read_nr) = 1
    '''
    c_brack = get_bracken_classified_read_nr(bracken_path, bs_id, krakdb, readpool, mhg, cs, sf, r, level)
    k2_sampled_read_nr = get_kraken2_read_nr(kraken2_path, bs_id, krakdb, readpool, mhg, cs, sf, r)
    k2_total_read_nr = get_kraken2_read_nr(kraken2_path, bs_id, krakdb, readpool, mhg, cs, 1, 1)
    mapped_read_nr = get_total_read_nr_flagstat(bowtie2_path, bs_id, idx_name) - k2_total_read_nr ## 
    return c_brack + (k2_sampled_read_nr / k2_total_read_nr) * mapped_read_nr 

# def get_uc_read_nr_kraken2(kraken_path, bs_id, krakdb, readpool, mhg, cs, sf, r, bracken_path, level):
#     kraken_ucread_nr, kraken_cread_nr = parse_kreport(kraken_path, bs_id, krakdb, readpool, mhg, cs, sf, r)
#     return kraken_ucread_nr + kraken_cread_nr - parse_breport(bracken_path, bs_id, krakdb, readpool, mhg, cs, sf, r, level, ['root'])['root']
# 
# def get_nonbee_read_nr_kraken2(kraken_path, bs_id, krakdb, readpool, mhg, cs, sf, r):
#     "simply sum of kraken2 classified and unclassified here"
#     kraken_ucread_nr, kraken_cread_nr = parse_kreport(kraken_path, bs_id, krakdb, readpool, mhg, cs, sf, r)
#     return kraken_ucread_nr + kraken_cread_nr

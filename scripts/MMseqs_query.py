#!/home/andrew/miniconda3/envs/pygemme/bin/python

import pyGEMME
from pathlib import Path
import os
import shutil
import argparse

def main():
    parser = argparse.ArgumentParser(description='Run mmseqs2 query')
    parser.add_argument('-j','--jobname', type=str, default = None, help='jobname')
    parser.add_argument('-o','--output_dir', type=str, help='output directory', default='examples')
    args = parser.parse_args()
    
    if args.jobname is None:
        raise ValueError('Please provide a jobname')
    if args.output_dir is None:
        raise ValueError('Please provide an output directory')

    jobname = args.jobname
    output_dir = args.output_dir

    results_dir = Path(f'{output_dir}/{jobname}')
    if not results_dir.exists():
        results_dir.mkdir(parents=True)

    query_sequence = pyGEMME.get_query_sequence_uniprot(jobname, results_dir)

    msa_lines = pyGEMME.run_mmseqs2(query_sequence, jobname, use_env=True, use_filter=True)
    if os.path.exists(results_dir.joinpath(f'{jobname}_env')):
        shutil.rmtree(results_dir.joinpath(f'{jobname}_env'))
    shutil.move(f'{jobname}_env', results_dir)

    raw_msa_file = results_dir.joinpath(f'{jobname}_env/uniref.a3m')
    clean_msa = pyGEMME.remove_inserts_a3m(jobname,raw_msa_file)
    clean_msa_file = results_dir.joinpath(f'{jobname}_env/uniref.fasta')
    with open(clean_msa_file,'w') as fid:
        fid.writelines(clean_msa)

if __name__ == '__main__':
    main()
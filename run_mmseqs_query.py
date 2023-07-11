from pyGEMME.mmseqs_query import get_query_sequence_uniprot, run_mmseqs2, remove_inserts_a3m
from pathlib import Path
import shutil
import argparse

def main():
    parser = argparse.ArgumentParser(description='Run mmseqs2 query')
    parser.add_argument('-j','--jobname', type=str, help='jobname')
    parser.add_argument('-o','--output_dir', type=str, help='output directory', default='examples')
    args = parser.parse_args()

    jobname = args.jobname
    output_dir = args.output_dir

    results_dir = Path(f'{output_dir}/{jobname}')
    if not results_dir.exists():
        results_dir.mkdir(parents=True)

    query_sequence = get_query_sequence_uniprot(jobname, results_dir)

    msa_lines = run_mmseqs2(query_sequence, jobname, use_env=True, use_filter=True)
    shutil.move(f'{jobname}_env', results_dir)

    raw_msa_file = results_dir.extend(f'{jobname}_env/uniref.a3m')
    clean_msa = remove_inserts_a3m(jobname,raw_msa_file)
    clean_msa_file = results_dir.extend(f'{jobname}_env/uniref.fasta')
    with open(clean_msa_file,'w') as fid:
        fid.writelines(clean_msa)

if __name__ == '__main__':
    main()
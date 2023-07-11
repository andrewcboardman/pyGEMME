#!/usr/bin/env python

from multiprocessing import Pool
import pyGEMME
from Bio import AlignIO
import pandas as pd
import argparse
import os


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-t','--test',default=False,action='store_true')
    parser.add_argument('-i','--input',type=str,default=None,action='store',help='MSA input file')
    parser.add_argument('-o','--output',type=str,default=None,action='store')
    parser.add_argument('-j','--jet_trace_file',type=str,default=None,action='store',help='JET trace file')
    parser.add_argument('-p','--plot_output',default=False,action='store_true',help='Plot output')
    parser.add_argument('--use_jet_trace',default=False,action='store_true',help='Use JET trace')
    args = parser.parse_args()

    if args.test:
        msa_file = 'examples/blat_ecolx/BLAT_ECOLX.fasta'
        jet_res_file = 'examples/blat_ecolx/results/trace.csv'
        output_dir = 'examples/blat_ecolx/results'

    else:
        if args.input is None:
            raise ValueError('Please specify an MSA input file!')
        msa_file = args.input # Set MSA filepath

        if args.use_jet_trace & (args.jet_trace_file is None):
            raise ValueError('Please specify an JET input file or disable use_jet_trace!')
        elif args.use_jet_trace:
            jet_res_file = args.jet_trace_file

        if args.output is None:
            raise ValueError('Please specify an output path!')
        output_dir = args.output

    ## Load MSA
    msa = next(AlignIO.parse(msa_file,format='fasta'))

    ## Check if output directory exists
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    if args.use_jet_trace:
        ## Load pre-computed JET trace
        trace = pd.read_csv(jet_res_file)
    else:
        ## Estimate evolutionary trace 
        trace = pyGEMME.estimate_trace(msa)
        trace.to_csv(f'{output_dir}/trace.csv',index=False)

    ## Estimate evolutionary distances using trace
    msa_array, msa_array_binary, d_evol = pyGEMME.estimate_distances(msa, trace.trace.values)
    d_evol.to_csv(f'{output_dir}/d_evol.csv',index=False)

    ## Predict GEMME fitness scores
    fitness = pyGEMME.predict_fitness(msa_array, trace.trace.values, d_evol.d_evol.values)
    fitness.to_csv(f'{output_dir}/fitness.csv',index=False)

    ## Visualise GEMME predictions (optional)
    if args.plot_output or args.test:
        pyGEMME.make_plots(trace, d_evol, fitness, output_dir)

if __name__ == '__main__':
    main()
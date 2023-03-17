"""make diagnostic plots to check JET output"""
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
from Bio import Align, AlignIO, Seq, SeqRecord
from Bio.Phylo import draw
from pyGEMME import trace, model
from scipy.stats import spearmanr

def test_trace():
    msa_file = 'examples/blat_ecolx/aliBLAT.fasta'
    trace_output_file = 'tests/blat_ecolx_trace.csv'
    msa, msa_subsamples, trees, traces, mean_trace_sig = trace.calculate_mean_trace(msa_file, trace_output_file)

    pc_id = trace.pc_id_to_query(msa)
    _ = plt.hist(pc_id, bins=50, alpha = 0.5, label='Full Alignment')
    msa_subsample_w_query = Align.MultipleSeqAlignment([msa[0]] + [sr for sr in msa_subsamples[0]])

    _ = plt.hist(trace.pc_id_to_query(msa_subsample_w_query),bins=50, alpha = 0.5, label='Subsample')

    trees[0].ladderize()
    plt.rc('font', size=0) 
    draw(trees[0])

    # Check consensus sequences of tree
    depth_index = trace.get_depth_index(trees[0])
    consensus_sequences, consensus_node_names = trace.get_consensus_tree(msa_subsamples[0], trees[0], depth_index)

    consensus_aln = []
    for clade in trees[0].find_clades():
        sr = SeqRecord.SeqRecord(
            seq= Seq.Seq(''.join(consensus_sequences[consensus_node_names.index(clade.name)])),
            id = clade.name)
        consensus_aln.append(sr)

    consensus_aln = Align.MultipleSeqAlignment(consensus_aln)
    with open('tests/test_subsample_msa_consensus.fa','w') as fid:
        AlignIO.write([consensus_aln],fid,format='fasta')

    # Check trace calculation at a chosen depth 
    trace = np.zeros(consensus_sequences.shape[1])
    max_depth = depth_index.depth.max()
    for depth in range(max_depth):
        ids_by_depth = depth_index[depth_index.depth==depth].name.values      
        ids_by_depth = [id for id in ids_by_depth if id.startswith('Inner')]
        consensus_by_depth = consensus_sequences[[consensus_node_names.index(id) for id in ids_by_depth]]
        trace[(trace<1) & ((consensus_by_depth != '-').sum(axis=0) > 1)] = depth

    trace[trace<1] = max_depth

    # visualise a chosen depth
    depth = 3
    ids_by_depth = depth_index[depth_index.depth==depth].name.values      
    ids_by_depth = [id for id in ids_by_depth if id.startswith('Inner')]
    consensus_by_depth = consensus_sequences[[consensus_node_names.index(id) for id in ids_by_depth]]
    for id in ids_by_depth:
        print(id)
        print(''.join(consensus_sequences[consensus_node_names.index(id)]))
    print(''.join(np.where(trace == depth, 'X','-')))

    # Check traces against JET results
    jet_res_file = 'tests/BLAT_model/BLAT_jet.res'
    df = pd.read_csv(jet_res_file,delim_whitespace=True)
    plt.subplots()

    plt.rcdefaults()
    plt.plot(df.trace)
    for trace_sig in traces:
        plt.plot(trace_sig, alpha = 0.1)
    
    plt.plot(mean_trace_sig)
    print(spearmanr(df.trace, mean_trace_sig))

    # check for convergence to mean value
    corr = []
    for i in range(2,len(traces)):
        mean_trace_sig_partial	= np.stack(traces[:i]).mean(axis=0)
        corr.append(spearmanr(mean_trace_sig, mean_trace_sig_partial)[0])

    plt.subplots()
    plt.plot(corr)

    if max(corr) > 0.99:
        iter_thresh = min([i+1 for i in range(len(corr)) if corr[i]>0.99])
        print(f'99% correlation to mean value after {iter_thresh} iterations')
    else:
        print(f'Has not converged after {len(traces)} iterations')


def test_model_concordance():

    msa_file = 'tests/BLAT_model/BLAT_A.fasta'
    msa = next(AlignIO.parse(msa_file,format='fasta'))
    jet_res_file = 'tests/BLAT_model/BLAT_jet.res'
    jet_trace = pd.read_csv(jet_res_file,delim_whitespace=True)
    predict_tracejet, d_evol = model.predict_gemme(msa, jet_trace.trace.values)


    gemme_predictions_tidy = pd.read_csv('tests/BLAT_model/model/BLAT_predictions.csv',index_col=0)
    compare = gemme_predictions_tidy.merge(predict_tracejet, on = ['wt','pos','mut'])
    print(compare[['combi_x','combi_y']].corr().iloc[0,1])
    print(compare[['norm_epi_x','norm_epi_y']].corr().iloc[0,1])
    print(compare[['norm_ind_x','norm_ind_y']].corr().iloc[0,1])
    print(compare[['epi_x','epi_y']].corr().iloc[0,1])
    print(compare[['ind_x','ind_y']].corr().iloc[0,1])
    d_evol_gemme = pd.read_csv('tests/BLAT_model/model/BLAT_evolDist.mat',delim_whitespace=True).V1.values
    print(np.max(d_evol[1:]-d_evol_gemme))
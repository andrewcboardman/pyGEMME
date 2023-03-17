# -*- coding: utf-8 -*-

# Copyright (c) 2018: Elodie Laine
# This code is part of the gemme package and governed by its license.
# Please see the LICENSE.txt file included as part of this package.

import os
import re
import subprocess
import shutil
from Bio import AlignIO
import numpy as np
import pandas as pd
from scipy import stats
from scipy.spatial.distance import cdist
import matplotlib.pyplot as plt
from multiprocessing import Pool
import re
from functools import partial
from tqdm import tqdm
from Bio import Align, Seq, SeqRecord
from Bio.Phylo import draw
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from typing import Literal, Any

alphabet: list[str] = [c.upper() for c in ("a","c","d","e","f","g","h","i","k","l","m","n","p","q","r","s","t","v","w","y")]

from .find_resources import find_path, find_template_conf, find_jet_matrices, find_tool

def edit_jet_config(input_file, output_file, param_dict) -> None:
	with open(input_file,'r') as fid:
		input_config: list[str] = [line.rstrip() for line in fid]
	output_config = []
	
	for line in input_config:
		line_: list[str] = line.split('\t')
		change_line: bool = False
		# Check if we need to change this line; if so prepare a new line
		for key, value in param_dict.items():
			if line_[0] == key:
				change_line = True
				new_line = '\t\t'.join([line_[0], str(value)] + line_[2:])
		# add a new line if we need to
		if change_line:
			output_config.append(new_line)
		else:
			output_config.append(line)

	with open(output_file,'w') as fid:
		for line in output_config:
			fid.write(line + '\n')	


def extract_fasta_alignment_query(input_file:str) -> tuple[str]:
	'''Extract the reference sequence from the input alignment'''
	with open(input_file,"r") as fIN:
		header = next(fIN)
		if header[0]!=">":
			raise Exception('bad FASTA format')
		else:
			reference_name	= re.compile("[^A-Z0-9a-z]").split(header[1:])[0]

		reference_seq = ''
		line = next(fIN)
		n_lines = 0
		while line[0] != '>':
			reference_seq += line.strip().strip(".").strip("-")
			n_lines += 1
			line = next(fIN)
	
	return reference_name, reference_seq


def create_fasta_file(seq, name, output_file):
	with open(output_file,'w') as fid:
		fid.write(f'>{name}\n')
		fid.write(seq+'\n')


def create_pdb_file(seq, output_file):
	"""Create a PDB file with dummy CA atoms based on the reference sequence"""
	d = {'C': 'CYS', 'D': 'ASP', 'S': 'SER',  'Q': 'GLN', 'K': 'LYS',
	'I': 'ILE', 'P': 'PRO', 'T': 'THR', 'F': 'PHE', 'N': 'ASN', 
	'G': 'GLY',  'H': 'HIS', 'L': 'LEU', 'R': 'ARG', 'W': 'TRP', 
        'A': 'ALA', 'V': 'VAL', 'E': 'GLU', 'Y': 'TYR', 'M': 'MET'}

	fOUT = open(output_file,'w')
	i = 1
	for let in seq:
		fOUT.write('ATOM%7d  CA  %s A%4d      43.524  70.381  46.465  1.00   0.0\n'%(i,d[let],i))
		i += 1
	fOUT.close()
	return None


def sanitise_input_fasta_alignment(fFile, output_file):
	"""Remove tabs and spaces from names of aligned sequences"""
	with open(fFile,'r') as fin:
		with open(output_file,'w') as fout:
			names = []
			for line in fin:
				if line[0] == '>':
					name = line[1:]
					name = name.split('\t')[0]
					name = name.split(' ')[0]
					while name in names:
						name = name + '_1'
					fout.write('>' + name + '\n')
				else:
					fout.write(line)


def JET_realign(
        input_file,
		output_dir,
		params):
	""""Run JET to get a tree representation of alignment"""
	
	n_iter = params['n_iter']
	N_seqs_max = params['N_seqs_max']
	retrieval_method=params['retrieval_method']	
	if retrieval_method == 'input':
		input_format = params['input_format']
	else:
		raise Exception('Only local alignment input currently supported')

	# Get configuration file
	# Print and adjust maximum sequence number
	template_jet_config_file = find_template_conf()
	jet_config_file = f"{output_dir}/jet.conf"
	jet_config_params = {
		'results':N_seqs_max,
		'muscle':find_tool('muscle'),
		'substMatrix':find_jet_matrices()
	}
	edit_jet_config(template_jet_config_file, jet_config_file, jet_config_params)
	
	# Get the query sequence 
	query_name, query_seq = extract_fasta_alignment_query(input_file)
	query_sequence_file = f'{output_dir}/{query_name}_query.fasta'
	create_fasta_file(query_seq, query_name, query_sequence_file)

	# Make a dummy PDB file for JET
	dummy_pdb_file = f'{output_dir}/{query_name}.pdb'
	create_pdb_file(query_seq, dummy_pdb_file)	
	print(f'Made dummy PDB file at {dummy_pdb_file}...')


	# Get input alignment
	if input_format == 'fasta':
		# sanitise FASTA inputs
		clean_alignment_file = f'{output_dir}/{query_name}_A.fasta'
		sanitise_input_fasta_alignment(input_file, clean_alignment_file)
	else:
		raise Exception('Only FASTA input currently supported')
	# # BLAST input
	# 	clean_alignment_file = f'{output_dir}/{query_name}_A.psiblast'
	# 	subprocess.call(f"cp {input_file} {clean_alignment_file}",shell=True)
	# 	jet_input = f" -r input -b {clean_alignment_file}"

	# Logging
	jet_output_file = f'{output_dir}/{query_name}.out'

	# Build JET command
	jet_locations = f'{find_path()}:{find_path()}/jet/extLibs/vecmath.jar'
	jetcmd = ["java",
		"-Xmx1000m",
		"-cp", jet_locations,
		"jet.JET",
		"-p", "J",
		"-d", "chain",
		'-n', str(n_iter),
		'-r', retrieval_method,
		'-c', jet_config_file,
		'-i', dummy_pdb_file,
		'-f', clean_alignment_file,
		'-o', output_dir, 
		'>', jet_output_file
	]
	
	reCode=subprocess.call(' '.join(jetcmd),shell=True)

	# Clean up results
	JET_output_dir = f"{output_dir}/{query_name}"
	JET_results_file = f"{JET_output_dir}/{query_name}_jet.res"
	os.remove(f"{output_dir}/{query_name}.pdb")
	if os.path.isfile(JET_results_file):
		os.rename(JET_results_file,f"{output_dir}/{query_name}_jet.res")
	if os.path.exists(JET_output_dir):
		shutil.rmtree(JET_output_dir)


	# make output files
	output = {
		'query_name':query_name,
		'query_sequence_file':query_sequence_file,
		'jet_alignment_file': f"{output_dir}/{query_name}_A.fasta",
		'jet_results_file':f"{output_dir}/{query_name}_jet.res",		
		'input_alignment':input_file,
		'jet_cmd': jetcmd
	}

	return output

###### Functions for computation of trace in Python
def pc_id_to_query(msa,alphabet= alphabet, gaps_adjust = False, position_weights=None):
    """Percentage ID of MSA sequences to first sequence in MSA using scipy.spatial.cdist"""
    msa_array = np.array(msa)

    N = len(msa)
    L = len(msa[0])
    nchar = len(alphabet)

    msa_array_binary = np.zeros((*msa_array.shape,nchar),dtype='int8')
    for i, a in enumerate(alphabet):
        msa_array_binary[:,:,i] = msa_array == a    
    msa_array_binary_reshape = np.reshape(msa_array_binary, (N, L * nchar)) 

    if position_weights is not None:
        # raise Exception('Not implemented yet!')
        # if len(position_weights) != L:
        #     raise Exception('Weights have different length')
        position_weights_char = np.repeat(position_weights, repeats=nchar)

        # Number of mismatches per residue
        # This counts 2 for a mismatch and 1 for a gap
        # Correct this later to just count mismatches / L 
        dist_to_query = cdist(
            msa_array_binary_reshape[0].reshape(1, L * nchar), 
            msa_array_binary_reshape, 
            metric='hamming',
            w = position_weights_char
            )
        
        n_gaps = (~msa_array_binary.any(axis=2) * position_weights).sum(axis=1)
        

    else:
        dist_to_query = cdist(
            msa_array_binary_reshape[0].reshape(1, L * nchar), 
            msa_array_binary_reshape, 
            metric='hamming'
            )
        n_gaps = (~msa_array_binary.any(axis=2)).sum(axis=1)
    dist_to_query = dist_to_query[0] * L * nchar
    
    # Adjust for number of gaps in aligned sequences
    if gaps_adjust:
        dist_to_query = dist_to_query - n_gaps
    
    # Adjust for doublecount of mismatches 
    dist_to_query = dist_to_query / 2

    # Final pc id
    pc_id = 1 - dist_to_query / L
    return pc_id


def sample_by_bin(bin_labels, bins, n_sample = 100):
	"""Stratified sampling of MSA by identity to query

	Args:
		bin_labels (_type_): _description_
		bins (_type_): _description_
		n_sample (int, optional): _description_. Defaults to 100.

	Yields:
		_type_: _description_
	"""
	for bin in bins:
		bin_index = np.where(bin_labels == bin)[0]
		n_sample = min(n_sample, len(bin_index))
		bin_index_sample = list(np.random.choice(bin_index, n_sample, replace=False))
		yield bin_index_sample


def strat_subsample_msa(msa, N_sample, pc_id,  bins = [0.2,0.4,0.6,0.8,0.9], random_seed = 1):
    """Stratified subsample from MSA, uniformly between each of the bins defined by bins"""
    np.random.seed(random_seed)
    N_sample_bin = N_sample // (len(bins) - 1)
    pc_id_binned = pd.cut(pc_id,bins = bins)
    sample_index = []
    for cat_sample_index in sample_by_bin(pc_id_binned, pc_id_binned.categories, N_sample_bin):
        sample_index.extend(cat_sample_index)        
    msa_sample = Align.MultipleSeqAlignment([msa[int(i)] for i in sample_index])
    return msa_sample


def random_subsample_msa(msa, N_sample, bins = [0.2,0.4,0.6,0.8,0.9], random_seed = 1):
    """random subsample from MSA, not biased by percent ID to query"""
    np.random.seed(random_seed)
    sample_index = np.random.choice(np.arange(len(msa)), N_sample, replace=False)
    msa_sample = Align.MultipleSeqAlignment([msa[int(i)] for i in sample_index])
    return msa_sample


def get_tree(msa_sample, method = 'nj', distance_model = 'blosum62'):
    """construct a phylogenetic tree from a subsampled MSA"""
    calculator = DistanceCalculator(distance_model)
    constructor = DistanceTreeConstructor(calculator, method)
    tree = constructor.build_tree(msa_sample)
    return tree


def get_consensus(msa_array):
    """From a numpy array representing an MSA, find a consensus sequence"""
    # one-hot encode residues 
    msa_array_binary = np.zeros((*msa_array.shape,20),dtype='bool')
    for i, a in enumerate(alphabet):
        msa_array_binary[:,:,i] = msa_array == a

    # find consensus columns
    consensus = (msa_array_binary.mean(axis=0)==1)

    # index into array of amino acids
    # represent lack of consensus with a gap 
    aa_array = np.repeat(np.array(alphabet + ['-']).reshape(1, 21), 263, axis=0)
    consensus_plus = np.concatenate((
        consensus, 
        ~consensus.any(axis=1).reshape((263, 1))
    ), axis=1) # L, 21
    return aa_array[consensus_plus]


def get_depth_index(tree):
    """get a dataframe containing depth labels for each node of a tree"""
    depths = tree.depths(unit_branch_lengths=True)
    df = pd.DataFrame(dict(
        name = [branch.name for branch in depths.keys()],
        depth = depths.values()
    ))
    df = df.sort_values('depth')
    return df


def lookup_by_names(tree):
    """create a lookup dictionary from names to clades"""
    names = {}
    for clade in tree.find_clades():
        if clade.name:
            if clade.name in names:
                raise ValueError("Duplicate key: %s" % clade.name)
            names[clade.name] = clade
    return names


def get_consensus_tree(msa, tree, depth_index):
    """Using an MSA and an inferred tree, walk up from the deepest branches and generate consensus sequences for inner nodes"""  
    # lookup table for subsampled MSA
    lookup = lookup_by_names(tree)

    # Make an array for inner node consensus sequences
    L = len(msa[0])
    inner_nodes = [clade for clade in tree.find_clades() if not clade.is_terminal()]
    n_inner_nodes = len(inner_nodes)
    inner_node_names = [clade.name for clade in inner_nodes]
    consensus_sequences = np.concatenate((
        np.array(msa),
        np.zeros((n_inner_nodes,L), dtype='<U1')
    ))
    msa_ids = [sr.id for sr in msa]
    consensus_node_names = msa_ids + inner_node_names

    # iterate over possible depths, starting from longest branch
    # compute consensus sequences at each inner node
    # represent consensus gap or lack of consensus with a gap character
    for node in depth_index.sort_values('depth',ascending=False).name.values:
        if node.startswith('Inner'):
            child_ids = [clade.name for clade in lookup[node]]
            child_sequences = consensus_sequences[[consensus_node_names.index(child_id) for child_id in child_ids]]
            consensus_sequences[consensus_node_names.index(node)] = get_consensus(child_sequences)

    return consensus_sequences, consensus_node_names


def get_trace(msa, tree):
    """Use the consensus sequences for the tree to assign a trace to each residue:
	this is the level of the tree at which this residue can be used to distinguish clades"""
    # Index tree by depth
    depth_index = get_depth_index(tree)

    # Get consensus sequences for inner nodes
    consensus_sequences, consensus_node_names = get_consensus_tree(msa, tree, depth_index)

    # iterate over possible depths, starting from root
    # get all ids at this depth and filter for inner nodes
    # compare their consensus sequences
    # if more than 1 tree has a consensus character at position i and depth d
    # and position i has not had a trace level assigned yet
    # position i has trace level d
    trace = np.zeros(consensus_sequences.shape[1])
    max_depth = depth_index.depth.max()
    for depth in range(max_depth):
        ids_by_depth = depth_index[depth_index.depth==depth].name.values      
        ids_by_depth = [id for id in ids_by_depth if id.startswith('Inner')]
        consensus_by_depth = consensus_sequences[[consensus_node_names.index(id) for id in ids_by_depth]]
        trace[(trace<1) & ((consensus_by_depth != '-').sum(axis=0) > 1)] = depth

    trace[trace<1] = max_depth
    # Correct for error in JET implementation
    trace[trace< max_depth] -= 1
    trace_sig = (max_depth - trace) / max_depth
    return trace_sig


def get_subsample_trace(msa, pc_id, N_sample, random_seed, sample_method = 'random', tree_method = 'nj', tree_metric = 'blosum62'):
	"""Subsample an MSA, given percentage ID, and """

	# Get subsample from MSA
	print('Subsample', random_seed)
	if sample_method == 'random':
		msa_subsample = random_subsample_msa(msa, N_sample, random_seed = random_seed)
	elif sample_method == 'stratified':
		msa_subsample = strat_subsample_msa(msa, N_sample, pc_id, random_seed= random_seed)

	# Create tree
	tree = get_tree(msa_subsample, method = tree_method, metric = tree_metric)

	# Get trace from subsampled MSA
	trace_sig = get_trace(msa_subsample, tree)

	return msa_subsample, tree, trace_sig


def calculate_mean_trace(msa_file, output_file, parallel = True, n_workers=10, n_tasks=20):
	"""Calculate mean trace by taking n_tasks (default = 20) subsamples from MSA"""
	msa = next(AlignIO.parse(msa_file,format='fasta'))
	for sr in msa:
		sr.seq = re.sub('\.','-',str(sr.seq))
	N = len(msa)
	N_sample = round(np.sqrt(N))

	# Get id to query sequence
	pc_id = pc_id_to_query(msa)
	print('Got PC id. Subsampling...')

	get_subsample_trace_p = partial(get_subsample_trace, msa, pc_id, N_sample)
	if parallel:
		with Pool(n_workers) as p:
			subsample_traces = p.map(get_subsample_trace_p, range(n_tasks))

	msa_subsamples = []
	trees = []
	traces = []
	for msa_subsample, tree, trace_sig in subsample_traces:
		msa_subsamples.append(msa_subsample)
		trees.append(tree)
		traces.append(trace_sig)


	mean_trace_sig = np.stack(traces).mean(axis=0)
	std_trace_sig = np.stack(traces).std(axis=0)

	df_trace = pd.DataFrame(dict(
		wt = msa[0], 
		pos = np.arange(1, len(msa[0])+1), 
		trace_mean = mean_trace_sig,
		trace_std= std_trace_sig
	))

	df_trace.to_csv(output_file)
	return msa, msa_subsamples, trees, traces, mean_trace_sig


def weighted_id(msa_file, jet_res_file, max_seqs = None, mode='query', alphabet=alphabet):
	""" Another function to compute weighted IDs, this time pairwise"""
	msa = list(AlignIO.parse(msa_file,format='fasta'))
	if max_seqs is None:
		max_seqs = len(msa)

	trace = pd.read_csv(jet_res_file,delim_whitespace=True)

	msa_ids = [sr.id for sr in list(msa[0])]
	msa_seqs = [sr.seq for sr in list(msa[0])]

	N = len(msa_seqs)
	L = len(msa_seqs[0])
	nchar = len(alphabet)

	msa_array = np.zeros((N,L),dtype='<U1')
	for i, seq in enumerate(msa_seqs):
		msa_array[i, :] = list(seq)

	msa_array[~np.isin(msa_array,alphabet)] = '-'
	msa_array_binary = np.zeros((*msa_array.shape,20),dtype='int8')

	for i, a in enumerate(alphabet):
		msa_array_binary[:,:,i] = msa_array == a

	position_weights = trace.trace.values ** 2
	assert position_weights.shape[0] == L


	msa_array_binary_reshape = np.reshape(msa_array_binary, (N, L * nchar))
	msa_array_weighted = msa_array_binary_reshape * position_weights.repeat(nchar)

	if mode == 'query':
		query = msa_array_weighted[0].reshape(1, L * nchar)
		distmat = cdist(query, msa_array_weighted[1:max_seqs], metric='cityblock')
		df = pd.DataFrame(dict(
			id1 = msa_ids[0],
			id2 = msa_ids[1:max_seqs],
			dist=distmat.squeeze()
		))
	elif mode == 'pairwise':
		distmat = cdist(msa_array_weighted[:max_seqs], msa_array_weighted[:max_seqs], metric='cityblock')
		r, c = np.triu_indices(max_seqs, 1, max_seqs)
		df = pd.DataFrame(dict(
			id1 = [msa_ids[r_] for r_ in r],
			id2 = [msa_ids[c_] for c_ in c],
			dist=distmat[r,c]
		))
	return df
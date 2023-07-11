import time
import numpy as np
from typing import List
from scipy.spatial.distance import cdist
import pandas as pd
import re
from functools import partial
from multiprocessing import Pool
from Bio import AlignIO, Align
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

###### Functions for computation of trace in Python

ALPHABET: List[str] = [c.upper() for c in ("a","c","d","e","f","g","h","i","k","l","m","n","p","q","r","s","t","v","w","y")]

def pc_id_to_query(msa, alphabet= ALPHABET, gaps_adjust = False, position_weights=None):
    """Fast percentage ID of MSA sequences to first sequence in MSA using scipy.spatial.cdist"""
    msa_array = np.array(msa)

    N = len(msa)
    L = len(msa[0])
    nchar = len(alphabet)

    # Convert to one-hot array of integers, looping over alphabet
    msa_array_binary = np.zeros((*msa_array.shape,nchar),dtype='int8')
    for i, a in enumerate(alphabet):
        msa_array_binary[:,:,i] = msa_array == a    
    
    # reshape to 2D array
    msa_array_binary_reshape = np.reshape(msa_array_binary, (N, L * nchar)) 

    if position_weights is None:
        position_weights = np.ones(L)
    if len(position_weights) != L:
            raise Exception('Weights have different length to MSA!')
    
    # Repeat weights for each character
    position_weights_char = np.repeat(position_weights, repeats=nchar)
    
    # Calculate Hamming (edit) distance to query sequence
    # This counts 2 for a mismatch and 1 for a gap, as we have one-hot encoded the sequences
    dist_to_query = cdist(
        msa_array_binary_reshape[0].reshape(1, L * nchar), 
        msa_array_binary_reshape, 
        metric='hamming',
        w = position_weights_char
        )
    
    # Un-normalise for number of characters
    # Adjust for doublecount of mismatches
    dist_to_query = dist_to_query[0] * nchar / 2
        
    if gaps_adjust:
        # Number of gaps per residue; if enabled, use this to correct for gaps counted above
        n_gaps = (~msa_array_binary.any(axis=2) * position_weights).sum(axis=1)
        dist_to_query = dist_to_query - n_gaps / 2

    # Convert dissimilarity to percent identity
    pc_id = 1 - dist_to_query
    return pc_id


def sample_by_bin(bin_labels, bins, n_sample = 100):
    """Stratified sampling of MSA by identity to query
    """
    for bin in bins:
        bin_index = np.where(bin_labels == bin)[0]
        n_sample = min(n_sample, len(bin_index))
        bin_index_sample = list(np.random.choice(bin_index, n_sample, replace=False))
        yield bin_index_sample


def subsample_msa(msa, N_sample, random_seed = 1, bins = [0.2,0.4,0.6,0.8,0.9], sample_method = 'random',pc_id = None):
    """subsample an MSA, either randomly or stratified by identity to query"""
    np.random.seed(random_seed)

    if sample_method == 'random':
        sample_index = np.random.choice(np.arange(len(msa)), N_sample, replace=False)
        msa_sample = Align.MultipleSeqAlignment([msa[int(i)] for i in sample_index])
    elif sample_method == 'stratified':
        if pc_id is None:
            pc_id = pc_id_to_query(msa)
        N_sample_bin = N_sample // (len(bins) - 1)
        pc_id_binned = pd.cut(pc_id,bins = bins)
        sample_index = []
        for cat_sample_index in sample_by_bin(pc_id_binned, pc_id_binned.categories, N_sample_bin):
            sample_index.extend(cat_sample_index)        
        msa_sample = Align.MultipleSeqAlignment([msa[int(i)] for i in sample_index])
    return msa_sample


def get_tree(msa_sample, method = 'nj', distance_model = 'blosum62'):
    """construct a phylogenetic tree from a subsampled MSA"""
    calculator = DistanceCalculator(distance_model)
    constructor = DistanceTreeConstructor(calculator, method)
    tree = constructor.build_tree(msa_sample)
    return tree


def get_consensus_sequence(msa_array, alphabet = ALPHABET):
    """From a numpy array representing an MSA, find a consensus sequence"""
    N, L = msa_array.shape
    nchar = len(alphabet)

    # one-hot encode residues 
    msa_array_binary = np.zeros((*msa_array.shape, nchar),dtype='bool')
    for i, a in enumerate(alphabet):
        msa_array_binary[:,:,i] = msa_array == a

    # find consensus columns
    consensus = (msa_array_binary.mean(axis=0)==1)

    # index into array of amino acids
    # represent lack of consensus with a gap 
    aa_array = np.repeat(np.array(alphabet + ['-']).reshape(1, nchar+1), L, axis=0)
    consensus_plus = np.concatenate((
        consensus, 
        ~consensus.any(axis=1).reshape((L, 1))
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
            consensus_sequences[consensus_node_names.index(node)] = get_consensus_sequence(child_sequences)

    return consensus_sequences, consensus_node_names


def get_trace(msa, tree):
    """The evolutionary trace of a position in a phylogenetic tree is defined by the minimum depth of the tree at which the amino acid at that position is fixed. Positions with a trace of 0 (and a trace level of 1) are perfectly conserved across the whole tree, while positions with a trace of max_depth (and a trace level of 0) are not conserved across even the smallest subtree.
    To calculate evolutionary trace for a tree, we therefore need to calculate consensus sequences for each inner node of the tree, and then work down from the root to the deepest to find the point at which each position becomes variable.
    """

    # Index tree by depth
    depth_index = get_depth_index(tree)

    # Get consensus sequences for inner nodes
    consensus_sequences, consensus_node_names = get_consensus_tree(msa, tree, depth_index)
    
    # compare their consensus sequences
    # if more than 1 tree has a consensus character at position i and depth d
    # and position i has not had a trace level assigned yet
    # position i has trace level d
    trace = np.zeros(consensus_sequences.shape[1])

    # iterate from the root to the maximum depth of the tree
    max_depth = depth_index.depth.max()
    for depth in range(max_depth):
        # Select inner nodes of the tree at this depth
        ids_by_depth = depth_index[depth_index.depth==depth].name.values      
        ids_by_depth = [id for id in ids_by_depth if id.startswith('Inner')]
        # get consensus sequences for these nodes
        consensus_by_depth = consensus_sequences[[consensus_node_names.index(id) for id in ids_by_depth]]
        # get all residues that are consensus characters at this depth
        consensus_residues = (consensus_by_depth != '-').sum(axis=0) > 1
        # Set trace level for residues that are consensus characters at this depth and have not been already assigned a trace level
        unassigned_trace_residues = trace < 1
        trace[unassigned_trace_residues & consensus_residues] = depth

    # Set unassigned trace levels to maximum depth
    trace[trace<1] = max_depth
    trace[trace < max_depth] -= 1
    trace_sig = (max_depth - trace) / max_depth

    return trace_sig


def get_subsample_trace(msa, pc_id, N_sample, random_seed, sample_method = 'random', tree_method = 'nj', tree_metric = 'blosum62'):
    """Subsample an MSA, """

    # Get subsample from MSA
    print('Subsample', random_seed)
    msa_subsample = subsample_msa(msa, N_sample, random_seed= random_seed, sample_method = sample_method, pc_id=pc_id)

    # Create tree
    tree = get_tree(msa_subsample, method = tree_method, distance_model = tree_metric)

    # Get trace from subsampled MSA
    trace_sig = get_trace(msa_subsample, tree)

    return msa_subsample, tree, trace_sig


def estimate_trace(msa, parallel = True, n_workers=10, n_tasks=20):
    """Calculate mean trace by taking n_tasks (default = 20) subsamples from MSA"""
    
    # Clean up MSA by converting '.' characters to '-'
    for sr in msa:
        sr.seq = re.sub('\.','-',str(sr.seq))
    
    # Get number of sequences in MSA
    N = len(msa)
    # Subsample size is sqrt(N)
    N_subsample = round(np.sqrt(N))

    # Get id to query sequence
    t0 = time.time()
    pc_id = pc_id_to_query(msa)
    t1 = time.time()
    print(f'Calculated all identities to query in {t1-t0:.2f}s.')

    # Subsample MSA, create tree, and get trace for each subsample
    print(f'Generating {n_tasks} subsamples of size {N_subsample}...')
    get_subsample_trace_p = partial(get_subsample_trace, msa, pc_id, N_subsample)
    if parallel:
        with Pool(n_workers) as p:
            subsample_traces = p.map(get_subsample_trace_p, range(n_tasks))
    else:
        subsample_traces = map(get_subsample_trace_p, range(n_tasks))

    # Get mean trace from subsamples
    msa_subsamples = []
    trees = []
    traces = []
    for msa_subsample, tree, trace_sig in subsample_traces:
        msa_subsamples.append(msa_subsample)
        trees.append(tree)
        traces.append(trace_sig)
    mean_trace_sig = np.stack(traces).mean(axis=0)
    std_trace_sig = np.stack(traces).std(axis=0)

    # Format as dataframe
    df_trace = pd.DataFrame(dict(
        wt = msa[0], 
        pos = np.arange(1, len(msa[0])+1), 
        trace = mean_trace_sig,
        trace_std= std_trace_sig
    ))

    return df_trace

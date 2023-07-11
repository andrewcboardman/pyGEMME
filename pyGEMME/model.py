
import numpy as np
import pandas as pd


alphabet = [c.upper() for c in ("a","c","d","e","f","g","h","i","k","l","m","n","p","q","r","s","t","v","w","y")]
coarsegrained_alphabet = [
    'IMVL',
    'FWY',
    'G',
    'P',
    'CAST',
    'NHQED',
    'RK'
]


def estimate_distances(msa, trace):
    """Estimate evolutionary distances weighted using trace"""
    msa_array = np.array(msa)
    N, L = msa_array.shape
    nchar = len(alphabet)

    if trace.shape[0] != L:
        raise ValueError('Trace length does not match MSA length!')

    
    msa_array_binary = np.zeros((N, L, nchar),dtype='bool')
    for i, aa in enumerate(alphabet):
        msa_array_binary[:,:,i] = msa_array == aa

    msa_array_binary_reshape = np.reshape(msa_array_binary, (N, L * nchar)) 
    position_weights = trace**2
    position_weights_char = np.repeat(position_weights, repeats=nchar)

    msa_array_weighted_positions = msa_array_binary_reshape * position_weights_char
    sim  = msa_array_weighted_positions @ msa_array_binary_reshape[0] 
    max_sim = (np.ones(L) * position_weights).sum()  
    d_evol = max_sim - sim

    df_d_evol = pd.DataFrame(dict(
        seq = [seq.id for seq in msa],
        d_evol = d_evol
    ))
    return msa_array, msa_array_binary, df_d_evol


def independent_model(msa_array, trace, alphabet = alphabet, reduced_alphabet = coarsegrained_alphabet, reduced = False):
    """Independent model"""
    N, L = msa_array.shape

    if not reduced:
        nchar = len(alphabet)
        msa_array_binary = np.zeros((N, L, nchar),dtype='bool')
        for i, aa in enumerate(alphabet):
            msa_array_binary[:,:,i] = msa_array == aa

        aa_counts = msa_array_binary.sum(axis=0)
        aa_counts_pseudo = np.maximum(aa_counts, 1)
        wt_counts = aa_counts[msa_array_binary[0]]
        
        ind_pred = np.log(aa_counts_pseudo / wt_counts.reshape(L, 1))
        norm_ind = ind_pred * trace.reshape((L, 1))

    if reduced:
        nchar_reduced = len(coarsegrained_alphabet)
        msa_array_binary_reduced = np.zeros((N, L, nchar_reduced),dtype=bool)
        for i, aa_set in enumerate(coarsegrained_alphabet):
            msa_array_binary_reduced[:,:,i] = np.isin(msa_array, list(aa_set))

        aa_counts_reduced = msa_array_binary_reduced.sum(axis=0)
        aa_counts_pseudo_reduced = np.maximum(aa_counts_reduced, 1)
        wt_counts_reduced = aa_counts_reduced[msa_array_binary_reduced[0]]

        # Convert to predictions for full alphabet
        nchar = len(alphabet)
        ind_pred_reduced = np.log(aa_counts_pseudo_reduced / wt_counts_reduced.reshape(L, 1))
        converter = np.zeros((nchar_reduced, nchar))
        for i, aa_set in enumerate(coarsegrained_alphabet):
            converter[i, :] = np.isin(np.array(alphabet), list(aa_set))
        ind_pred_reduced = ind_pred_reduced @ converter
        norm_ind = ind_pred_reduced * trace.reshape((L, 1))

    return norm_ind


def epistatic_model(msa_array, d_evol, trace,  thresh = 5):
    N, L = msa_array.shape
    nchar = len(alphabet)
    msa_array_binary = np.zeros((N, L, nchar),dtype='bool')
    for i, aa in enumerate(alphabet):
        msa_array_binary[:,:,i] = msa_array == aa

    # find all distances to query
    seq, pos, aa = np.where(msa_array_binary)   
    max_dist = d_evol.max()
    epi = np.ones((N, L, nchar)) * (max_dist)
    epi[seq, pos, aa] = d_evol[seq]

    # Minimum distance to query
    epi_min = epi[1:].min(axis=0) 
    epi_secondmin = np.partition(epi, 1, axis=0)[1]
    epi_min_combi = np.where(
        epi_secondmin - epi_min > thresh,
        epi_secondmin, epi_min
    )
    
    # Epistatic model
    aa_counts = msa_array_binary.sum(axis=0)
    epi_pred = np.select(
        condlist = [
            msa_array_binary[0],
            aa_counts<=1,
            aa_counts>1
        ],
        choicelist=[
            0,
            -100,
            epi_min_combi            
        ]
    )
      
    
    # the actual normalisation 
    norm_epi = epi_pred / epi_pred.max()
    norm_epi[norm_epi < 0] = 1
    norm_epi = -norm_epi * trace.reshape((L, 1))

    return norm_epi


def calc_norm_factor(msa_array, trace):
    """Extra normalisation factor for epistatic model, origin unknown"""
    N, L = msa_array.shape
    nchar = len(alphabet)
    msa_array_binary = np.zeros((N, L, nchar),dtype='bool')
    for i, aa in enumerate(alphabet):
        msa_array_binary[:,:,i] = msa_array == aa
    aa_counts = msa_array_binary.sum(axis=0)
    norm_factor = np.log(aa_counts.sum(axis=1).max())
    norm_factor *= L/np.sum(trace**2)
    return norm_factor


def predict_fitness(
        msa_array,
        trace,
        d_evol,
        thresh = 5, 
        alpha = 0.6,
        alphabet = alphabet,
        coarsegrained_alphabet = coarsegrained_alphabet
        ):
    """Python reimplementation of GEMME modelling"""
    N, L = msa_array.shape
    nchar = len(alphabet)

    # Independent model
    norm_ind = independent_model(msa_array, trace, alphabet, coarsegrained_alphabet, reduced = False)

    # Epistatic model
    norm_epi = epistatic_model(msa_array, d_evol, trace, thresh = thresh)

    # combine independent and epistatic models
    norm_ind_reduced = independent_model(msa_array, trace, alphabet, coarsegrained_alphabet, reduced = True)    
    norm_factor = calc_norm_factor(msa_array, trace)
    # Some bizarre normalisation 
    combi = alpha * norm_epi * norm_factor + (1 - alpha) * norm_ind_reduced  

    df = pd.DataFrame(dict(
        wt = msa_array[0].repeat(nchar),
        mut = np.tile(alphabet, L),
        pos = (np.arange(L) + 1).repeat(nchar),
        combi = combi.reshape(L * nchar),
        norm_epi = norm_epi.reshape(L * nchar),
        norm_ind = norm_ind.reshape(L * nchar),
    ))
    return df
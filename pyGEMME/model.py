
import numpy as np
import pandas as pd
from tqdm import tqdm
from itertools import product


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
    sim = msa_array_weighted_positions @ msa_array_binary_reshape[0] 
    max_sim = (np.ones(L) * position_weights).sum()  
    d_evol = max_sim - sim

    df_d_evol = pd.DataFrame(dict(
        seq = [seq.id for seq in msa],
        d_evol = d_evol
    ))
    return df_d_evol


def predict_fitness(
        msa,
        trace,
        d_evol,
        thresh = 5, 
        alpha = 0.6,
        alphabet = alphabet,
        coarsegrained_alphabet = coarsegrained_alphabet
        ):
    """Estimate conservation of MSA using tree traces"""
    N = len(msa)
    L = len(msa[0])
    nchar = len(alphabet)
    
    coarsegrain_dict = {}
    for aa in alphabet:
        for i, group in enumerate(coarsegrained_alphabet):
            if aa in list(group):
                coarsegrain_dict[aa] = i
    
    # Convert to numpy array
    msa_array = np.array(msa)
    # Count matches to wild type AA
    wt_counts = (msa_array == msa_array[0,:]).sum(axis=0)

    # Coarse grain MSA
    msa_array_coarse = np.zeros((N,L), dtype=int)
    for aa in alphabet:
        msa_array_coarse[np.where(msa_array==aa)] = coarsegrain_dict[aa]
    # Count matches to wild type AA group   
    wt_counts_coarse = (msa_array_coarse == msa_array_coarse[0,:]).sum(axis=0) 
    
    output = []
    
    for i, j in tqdm(product(range(L),range(nchar)), total = L * nchar):
        pos = i + 1
        aa_mut = alphabet[j]
        
        if aa_mut == msa_array[0,i]:
            output.append({'pos': pos,
            'aa_mut': aa_mut        
        })
            
        else:
                
            # count number of matches in full alphabet
            matches = msa_array[:,i]==aa_mut
            # add pseudocount
            mut_count_pseudo = np.maximum(matches.sum(), 1)
            # take log ratio of mut count to wt count
            ind_pred = np.log(mut_count_pseudo / wt_counts[i])
            
            # count number of matches in coarsegrained alphabet
            matches_reduced = msa_array_coarse[:,i] == coarsegrain_dict[aa_mut]
            # add pseudocount
            mut_count_reduced = np.maximum(matches_reduced.sum(), 1)
            # take log ratio of mut count to wt count
            ind_pred_reduced = np.log(mut_count_reduced / wt_counts_coarse[i])
            
            # distance for matches
            d_evol_match = d_evol[matches]
            
            if len(d_evol_match) == 0:
                # If no matches, set to max
                d_evol_min = np.max(d_evol)
            elif len(d_evol_match) == 1:
                # If one match, take that distance
                d_evol_min = np.min(d_evol_match)
            else:
                # If more than one match, check the two best
                # If second smallest distance is more than thresh larger than smallest,
                # Take second smallest as true smallest (???)
                d_evol_first, d_evol_second = list(np.partition(d_evol_match, 1)[:2])
                if d_evol_second - d_evol_min > thresh:
                    d_evol_min = d_evol_second
                else: 
                    d_evol_min = d_evol_first
                    
            output.append({
                'pos': pos,
                'aa_mut': aa_mut,
                'd_evol_min': d_evol_min,
                'ind_pred': ind_pred,
                'ind_pred_reduced': ind_pred_reduced         
            })
       
    df = pd.DataFrame(output)
    # Normalise epistatic prediction to be between 0 and 1
    df['d_evol_norm'] = -df.d_evol_min / np.max(df.d_evol_min)  
    
    # scale by max possible value for ind_pred
    gap_counts = (msa_array == '-').sum(axis=0)
    epi_norm_factor = np.log((N - gap_counts - wt_counts).max())
    df['epi_pred'] = df.d_evol_norm  * epi_norm_factor
    
    # scale all predictions by trace
    df['trace'] = df.pos.map(lambda x: trace[x-1])
    df['ind_pred'] = df.ind_pred * df.trace
    df['ind_pred_reduced'] = df.ind_pred_reduced * df.trace
    df['epi_pred'] = df.epi_pred * df.trace    
    
    # combine epistatic and reduced independent with weighting factor alpha
    df['combi'] = alpha * df.epi_pred + (1 - alpha) * df.ind_pred_reduced
    return df
import subprocess
import numpy as np
import pandas as pd
from Bio import SeqIO, AlignIO
import os 
import glob
from scipy import stats
from .trace import JET_realign
from .find_resources import (
    get_default_params,
    get_rscript_filename, 
    get_subs_matrix_filename, 
    get_alphabet_filename
)

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


# Run Rscript to compute predictions
def GEMME_models_fit_predict(
        jet_output,
        output_dir,
        params,
        mutations_file=None):
    
    query_name = jet_output['query_name']
    query_sequence_file = jet_output['query_sequence_file']
    input_alignment_file = jet_output['jet_alignment_file']
    jet_results_file = jet_output['jet_results_file']

    subs_matrix=params["subs_matrix"]
    alphabet=params["alphabet"]
    alpha = params["alpha"]

    model_path = get_rscript_filename('computePred')
    model_fn_path = get_rscript_filename('pred')
    subs_matrix_file = get_subs_matrix_filename(subs_matrix)
    alphabet_file = get_alphabet_filename(alphabet)

    rcmd = [
        "Rscript --save",
        model_path,
        query_name,
        query_sequence_file,
        input_alignment_file,
        jet_results_file,
        subs_matrix_file,
        alphabet_file,
        output_dir,
        model_fn_path
    ]

    if mutations_file is not None:
        rcmd.append("FALSE")
        rcmd.append(mutations_file)
    else:
        rcmd.append("TRUE")
        rcmd.append("none")
    rcmd = ' '.join(rcmd)
    reCode=subprocess.call(rcmd,shell=True)
    print('R command executed succesfully. Tidying output...')

    tidy_predictions = tidy_gemme_predictions(
        output_dir, 
        query_name, 
        query_sequence_file, 
        offset=0
    )
    tidy_predictions.to_csv(f'{output_dir}/{query_name}_predictions.csv')
    for file in glob.glob(f'{output_dir}/*.txt'):
        os.remove(file)

    return None


def tidy_prediction_matrix(prediction_matrix):
    df = pd.melt(prediction_matrix.reset_index(), id_vars='index')
    df.columns = ['mut','query_pos','pred']
    df['mut'] = df.mut.str.upper()
    df['query_pos'] = df.query_pos.str.slice(start=1).astype(int)
    df = df.set_index(['mut','query_pos'])
    return df.pred


def tidy_gemme_predictions(output_dir, query_name, query_seq_file, offset):
    read_r_matrix = lambda file: pd.read_csv(file, delim_whitespace=True)

    norm_pred_combi = read_r_matrix(f'{output_dir}/{query_name}_normPred_evolCombi.txt')
    norm_pred_epi = read_r_matrix(f'{output_dir}/{query_name}_normPred_evolEpi.txt')
    norm_pred_ind = read_r_matrix(f'{output_dir}/{query_name}_normPred_evolInd.txt')
    pred_epi = read_r_matrix(f'{output_dir}/{query_name}_pred_evolEpi.txt')
    pred_ind = read_r_matrix(f'{output_dir}/{query_name}_pred_evolInd.txt')
    pssm = read_r_matrix(f'{output_dir}/{query_name}_pssm.txt')
    pssm80 = read_r_matrix(f'{output_dir}/{query_name}_pssm80.txt')
    pssm60 = read_r_matrix(f'{output_dir}/{query_name}_pssm60.txt')

    gemme_predictions_tidy = pd.DataFrame(dict(
        combi=tidy_prediction_matrix(norm_pred_combi), 
        norm_epi=tidy_prediction_matrix(norm_pred_epi),
        norm_ind=tidy_prediction_matrix(norm_pred_ind),
        epi=tidy_prediction_matrix(pred_epi),
        ind=tidy_prediction_matrix(pred_ind),
        pssm=tidy_prediction_matrix(pssm),
        pssm80=tidy_prediction_matrix(pssm80),
        pssm60=tidy_prediction_matrix(pssm60),
        ))
    gemme_predictions_tidy = gemme_predictions_tidy.reset_index()

    # add conservation
    conservation = read_r_matrix(f'{output_dir}/{query_name}_conservation.txt')
    conservation = conservation.T
    conservation.index = conservation.index.astype(int)
    gemme_predictions_tidy.merge(conservation, left_on='query_pos',right_index=True)

    # add query residues
    sr = SeqIO.read(query_seq_file,format='fasta')
    wt_seq = list(sr.seq)
    assert len(wt_seq) == gemme_predictions_tidy.query_pos.max()
    gemme_predictions_tidy['wt'] = gemme_predictions_tidy.query_pos.map(lambda x: wt_seq[x-1])

    # Offset before naming mutants
    gemme_predictions_tidy['pos'] = gemme_predictions_tidy.query_pos + offset
    gemme_predictions_tidy['mutant'] = gemme_predictions_tidy.wt + gemme_predictions_tidy.pos.astype(str) + gemme_predictions_tidy.mut
    return gemme_predictions_tidy


def transform_fit_predict(alignment_file, output_dir, params=None, mutations_file=None):
    
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    
    if params is None:
        params = get_default_params()

    jet_output = JET_realign(
        alignment_file,
        output_dir=output_dir,
        params=params['jet']
        )

    GEMME_models_fit_predict(
        jet_output,
        output_dir=output_dir,
        params = params['model'],
        mutations_file=mutations_file
        )


def predict_gemme(
        msa, 
        trace,
        thresh = 5, 
        alpha = 0.6,
        alphabet = alphabet,
        coarsegrained_alphabet = coarsegrained_alphabet
        ):
    msa_array = np.array(msa)
    N, L = msa_array.shape
    nchar = len(alphabet)

    msa_array_binary = np.zeros((N, L, nchar),dtype='bool')
    for i, aa in enumerate(alphabet):
        msa_array_binary[:,:,i] = msa_array == aa

    aa_counts = msa_array_binary.sum(axis=0)
    aa_counts_pseudo = np.maximum(aa_counts, 1)
    wt_counts = aa_counts[msa_array_binary[0]]
    ind_pred = np.log(aa_counts_pseudo / wt_counts.reshape(L, 1))


    msa_array_binary_reshape = np.reshape(msa_array_binary, (N, L * nchar)) 
    position_weights = trace**2
    position_weights_char = np.repeat(position_weights, repeats=nchar)

    msa_array_weighted_positions = msa_array_binary_reshape * position_weights_char
    sim  = msa_array_weighted_positions @ msa_array_binary_reshape[0] 
    max_sim = (np.ones(L) * position_weights).sum()  
    d_evol = max_sim - sim

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
    
    epi_pred = np.select(
        condlist = [
            msa_array_binary[0],
            aa_counts<=1,
            aa_counts>1
        ],
        choicelist=[
            0,
            -100,# * 1.15, # don't like this but it works - could use NA
            epi_min_combi            
        ]
    )
    
    # Some bizarre normalisation 
    epi_pred_norm = epi_pred / np.sum(trace**2) 
    epi_pred_norm[epi_pred_norm < 0] = 1
    #epi_pred = np.nan_to_num(epi_pred,nan=1)

    
    # the actual normalisation 
    norm_epi = epi_pred / epi_pred.max()
    norm_epi[norm_epi < 0] = 1
    norm_epi = -norm_epi * trace.reshape((L, 1))

    # This is sensible thank Jesus
    norm_ind = ind_pred * trace.reshape((L, 1))

    # Weird ass reduced alphabet calculation
    nchar_reduced = len(coarsegrained_alphabet)
    msa_array_binary_reduced = np.zeros((N, L, nchar_reduced),dtype=bool)
    for i, aa_set in enumerate(coarsegrained_alphabet):
        msa_array_binary_reduced[:,:,i] = np.isin(msa_array, list(aa_set))
    aa_counts_reduced = msa_array_binary_reduced.sum(axis=0)
    aa_counts_pseudo_reduced = np.maximum(aa_counts_reduced, 1)
    wt_counts_reduced = aa_counts_reduced[msa_array_binary_reduced[0]]
    ind_pred_reduced = np.log(aa_counts_pseudo_reduced / wt_counts_reduced.reshape(L, 1))
    converter = np.zeros((nchar_reduced, nchar))
    for i, aa_set in enumerate(coarsegrained_alphabet):
        converter[i, :] = np.isin(np.array(alphabet), list(aa_set))
    ind_pred_reduced = ind_pred_reduced @ converter

    # combine 
    norm_epi_combine = epi_pred / epi_pred.max()
    norm_epi_combine[norm_epi_combine < 0] = 1
    norm_factor = -np.log(aa_counts.sum(axis=1).max())
    alpha = 0.6
    combi = alpha * norm_epi_combine * norm_factor + (1 - alpha) * ind_pred_reduced    
    norm_combi = combi * trace.reshape((L, 1))

    df = pd.DataFrame(dict(
        wt = msa_array[0].repeat(nchar),
        mut = np.tile(alphabet, L),
        pos = (np.arange(L) + 1).repeat(nchar),
        combi = norm_combi.reshape(L * nchar),
        norm_epi = norm_epi.reshape(L * nchar),
        norm_ind = norm_ind.reshape(L * nchar),
        epi = epi_pred_norm.reshape(L * nchar),
        ind = ind_pred.reshape(L * nchar),
    ))
    return df, d_evol


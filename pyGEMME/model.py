import subprocess
import pandas as pd
from Bio import SeqIO
import os 
import glob
from .run_jet import JET_realign
from .find_resources import (
    get_default_params,
    get_rscript_filename, 
    get_subs_matrix_filename, 
    get_alphabet_filename
)


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

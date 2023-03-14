import os
import shutil
import subprocess
from pyGEMME import model, find_resources, run_jet


def test_find_jet():
    path = find_resources.find_path()
    assert 'jet' in os.listdir(path)

def test_run_jet():
    path = find_resources.find_path()
    jet_locations = f'{path}:{path}/jet/extLibs/vecmath.jar'
    jetcmd = ['java', '-Xmx1000m', '-cp', jet_locations, 'jet.JET']
    proc = subprocess.Popen(jetcmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    proc.wait()
    assert proc.stderr.read().decode().split('\n')[0] == 'missing option(s)'

def test_JET_example():
    input_file = 'example/aliBLAT.fasta'
    output_dir = 'tests/BLAT_JET'
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    os.mkdir(output_dir)

    params = dict(
        n_iter=1,
        N_seqs_max=1000,
        retrieval_method='input',
        input_format='fasta'
    )

    jet_output = run_jet.JET_realign(
        input_file,
		output_dir,
        params
        )
    print(jet_output)
    assert jet_output['query_name'] == 'BLAT'

    
def test_model_example():
    query_name='BLAT'
    query_sequence_file = 'tests/BLAT_model/BLAT_query.fasta'
    jet_alignment_file = 'tests/BLAT_model/BLAT_A.fasta'
    jet_results_file = 'tests/BLAT_model/BLAT_jet.res'
    input_files = {
		'query_name':query_name,
		'query_sequence_file':query_sequence_file,
		'jet_alignment_file': jet_alignment_file,
		'jet_results_file':jet_results_file
	}

    output_dir = 'tests/BLAT_model/model/'
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    os.mkdir(output_dir)

    params = dict(
        subs_matrix="blosum62p", 
		alphabet="lz-bl.7",
        alpha=0.6
    )

    
    model.GEMME_models_fit_predict(
		input_files,
		output_dir,
		params = params,
		mutations_file=None)
    
    assert os.path.exists(output_dir + query_name + '_predictions.csv')

def test_pipeline():
    input_file = 'example/aliBLAT.fasta'
    output_dir = 'tests/BLAT_pipeline/'
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    os.mkdir(output_dir)


    model.transform_fit_predict(
        alignment_file=input_file, 
        output_dir=output_dir
    )
    assert os.path.exists(output_dir + 'BLAT_predictions.csv')


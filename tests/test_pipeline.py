import os
import shutil
import subprocess
from pyGEMME import model, jet, find_resources


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

    jet_output = jet.JET_realign(
        input_file,
		output_dir,
        n_iter=1,
        N_seqs_max=1000,
        mode='fasta',
        retrieval_method='input')
    print(jet_output)
    assert jet_output['query_name'] == 'BLAT'

    
def test_model_example():
    prot='BLAT'
    wt_sequence_file = 'tests/BLAT_model/BLAT_query.fasta'
    input_alignment_file = 'tests/BLAT_model/BLAT_A.fasta'
    jet_results_file = 'tests/BLAT_model/BLAT_jet.res'
    output_dir = 'tests/BLAT_model/model/'
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    os.mkdir(output_dir)

    model.GEMME_models_fit_predict(
		prot,
		wt_sequence_file,
		input_alignment_file,
		jet_results_file, 
		output_dir,
		subs_matrix="blosum62p", 
		alphabet="lz-bl.7",
		mutations_file=None)
    
    assert os.path.exists(output_dir + prot + '_conservation.txt')

import subprocess
import pkg_resources
from .jet import JET_realign
from .find_resources import find_path, get_subs_matrix_filename, get_alphabet_filename

# Run Rscript to compute predictions
def GEMME_models_fit_predict(
		prot,
		wt_sequence_file,
		input_alignment_file,
		jet_results_file, 
		output_dir,
		subs_matrix="blosum62p", 
		alphabet="lz-bl.7",
		mutations_file=None):
	
	model_path = str(find_path()) + '/computePred.R'
	model_fn_path = str(find_path()) + '/pred.R'
	subs_matrix_file = str(get_subs_matrix_filename(subs_matrix))
	alphabet_file = str(get_alphabet_filename(alphabet))

	rcmd = [
		"Rscript --save",
		model_path,
	  	prot,
		wt_sequence_file,
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
	return(reCode)

def transform_fit_predict(alignment_file, jobname, n_iter, N_seqs_max, mode, retrieval_method, jet_conf_file, model_subs_matrix, model_alphabet_file, mutations_file, keep_tmp_files):

    jet_output = JET_realign(
        alignment_file,
        output_dir=jobname,
        n_iter=n_iter,
        N_seqs_max=N_seqs_max,
        mode = mode,
        retrieval_method=retrieval_method,
        jet_conf_file=jet_conf_file,
        keep_tmp_files=keep_tmp_files
        )

    GEMME_models_fit_predict(
        jet_output['query_name'],
        jet_output['query_sequence_file'],
        jet_output['jet_alignment_file'],
        jet_output['jet_results_file'],
        output_dir=jobname,
        subs_matrix=model_subs_matrix, 
        alphabet_file=model_alphabet_file,
        mutations_file=mutations_file
        )
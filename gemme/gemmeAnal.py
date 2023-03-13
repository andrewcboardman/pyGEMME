# -*- coding: utf-8 -*-

# Copyright (c) 2018: Elodie Laine
# This code is part of the gemme package and governed by its license.
# Please see the LICENSE.txt file included as part of this package.

import os
import re
import subprocess
import shutil


def extract_query(input_file:str):
	'''Extract the reference sequence from the input alignment'''
	with open(input_file,"r") as fIN:
		header = next(fIN)
		if header[0]!=">":
			raise Exception('bad FASTA format')
		else:
			reference_name = re.compile("[^A-Z0-9a-z]").split(header[1:])[0]

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
		fid.write(seq)

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
        n_iter,
        N_seqs_max,
        mode,
        retrieval_method='input',
        jet_conf_file='custom.conf',
        keep_tmp_files=False):
	""""Run JET to get a tree representation of alignment"""

	# Get configuration file
	# Print and adjust maximum sequence number
	shutil.copy(jet_conf_file,f"{output_dir}/jet.conf")
	subprocess.call(f"sed -i 's/results\t\t5000/results\t\t{N_seqs_max}/' {output_dir}/jet.conf",shell=True) 

	# Get the sequence 
	query_name, query_seq = extract_query(input_file)
	query_sequence_file = f'{output_dir}/{query_name}_query.fasta'
	create_fasta_file(query_seq, query_name, query_sequence_file)

	# Make a dummy PDB file for JET
	dummy_pdb_file = f'{output_dir}/{query_name}.pdb'
	create_pdb_file(query_seq, dummy_pdb_file)	
	print(f'Made dummy PDB file at {dummy_pdb_file}...')
	
	# Build & run JET command
	if retrieval_method=="input":
		if mode == 'fasta':
			# sanitise FASTA inputs
			clean_alignment_file = f'{output_dir}/{query_name}_A.fasta'
			sanitise_input_fasta_alignment(input_file, clean_alignment_file)
		# BLAST input
		else:
			clean_alignment_file = f'{output_dir}/{query_name}_A.psiblast'
			subprocess.call(f"cp {input_file} {clean_alignment_file}",shell=True)
		# prepare jet command
		jet_input = f" -r input -f {clean_alignment_file}"
	# Other retrieval method (specified in JET)			
	else:
		jet_input = f"-r {retrieval_method}"

	jetcmd = ["java",
		"-Xmx1000m",
		"-cp .:./jet/extLibs/vecmath.jar",
		"jet.JET",
		f"-c {output_dir}/jet.conf",
		f"-i {dummy_pdb_file}",
		f"-o {output_dir}",
		"-p J",
		"-d chain",
		f"-n {n_iter}",
		jet_input,
		f"> {output_dir}/{query_name}.out"
	]
	jetcmd = ' '.join(jetcmd)
	print(jetcmd)
	reCode=subprocess.call(jetcmd,shell=True)

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
		'jet_alignment_file':query_name+'_A.fasta',
		'jet_results_file':query_name+"_jet.res",		
		'input_alignment':input_file,
	}

	return output

# Run Rscript to compute predictions
def GEMME_models_fit_predict(
		prot,
		wt_sequence_file,
		input_alignment_file,
		jet_results_file, 
		output_dir,
		blosum62_file="blosum62p.txt", 
		alphabet_file="./alphabets/lz-bl.7.txt",
		mutFile=None):

	rcmd = [
		"Rscript --save ./computePred.R",
	  	prot,
		wt_sequence_file,
		input_alignment_file,
		jet_results_file,
		blosum62_file,
		alphabet_file,
		output_dir
	]

	if mutFile is not None:
		rcmd.append("FALSE")
		rcmd.append(mutFile)
	else:
		rcmd.append("TRUE")
		rcmd.append("none")
	rcmd = ' '.join(rcmd)
	reCode=subprocess.call(rcmd,shell=True)
	return(reCode)

# Remove temporary files
def remove_tmp_files(prot,bFile,fFile):

	if bFile!='':
		if bFile!=prot+"_A.psiblast":
			os.remove(prot+"_A.psiblast")
	else:
		if os.path.isfile(prot+"/"+prot+"_A.psiblast"):
			os.rename(prot+"/"+prot+"_A.psiblast",prot+"_A.psiblast")
	if fFile!='':
		if fFile!=prot+"_A.fasta":
			if os.path.isfile(prot+"_A.fasta"):
				os.remove(prot+"_A.fasta")		
#	if os.path.isfile(prot+"/"+prot+"_jet.res"):
#		os.rename(prot+"/"+prot+"_jet.res",prot+"_jet.res")
	os.remove(prot+".pdb")
	dir_name = prot+"/"
	if os.path.isdir(dir_name):
		for f in os.listdir(dir_name): 
			f_path = os.path.join(dir_name, f)
			if os.path.isfile(f_path):
				os.remove(f_path)
		os.rmdir(dir_name)


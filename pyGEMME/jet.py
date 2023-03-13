# -*- coding: utf-8 -*-

# Copyright (c) 2018: Elodie Laine
# This code is part of the gemme package and governed by its license.
# Please see the LICENSE.txt file included as part of this package.

import os
import re
import subprocess
import shutil
from .find_resources import find_path, find_template_conf, find_jet_matrices, find_tool

def edit_jet_config(input_file, output_file, param_dict):
	with open(input_file,'r') as fid:
		input_config = [line.rstrip() for line in fid]
	output_config = []
	
	for line in input_config:
		line_ = line.split('\t')
		change_line = False
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


def extract_fasta_alignment_query(input_file:str):
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
        retrieval_method='input'):
	""""Run JET to get a tree representation of alignment"""

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
	if mode == 'fasta':
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

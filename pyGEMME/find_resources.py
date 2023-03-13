import pyGEMME
from importlib.resources import files


def find_path():
	return files(pyGEMME)

def find_template_conf():
	return files(pyGEMME).joinpath('jet.conf')

def find_tool(name):
	return files(pyGEMME.tools).joinpath(name)

def find_jet_matrices():
	return files(pyGEMME.jet_matrices)

def get_subs_matrix_filename(name='blosum62p'):
	return files(pyGEMME.matrix).joinpath(f'{name}.txt')

def get_alphabet_filename(name="lz-bl.7"):
	return files(pyGEMME.alphabets).joinpath(f'{name}.txt')
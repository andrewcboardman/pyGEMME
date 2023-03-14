import pyGEMME
import json
from importlib.resources import files


def find_path():
	return files(pyGEMME)

def get_default_params():
    with open(files(pyGEMME.data).joinpath("default_params.json"),'r') as fid:
        params = json.load(fid)
    return params

def find_template_conf():
	return files(pyGEMME.data).joinpath('jet.conf')

def find_tool(name):
	return files(pyGEMME.tools).joinpath(name)

def find_jet_matrices():
	return files(pyGEMME.data.jet_matrices)

def get_subs_matrix_filename(name='blosum62p'):
	return files(pyGEMME.data.matrix).joinpath(f'{name}.txt')

def get_alphabet_filename(name="lz-bl.7"):
	return files(pyGEMME.data.alphabets).joinpath(f'{name}.txt')

def get_rscript_filename(name='pred'):
	return str(files(pyGEMME.pred).joinpath(f'{name}.R'))

def get_subs_matrix_filename(name='blosum62p'):
	return str(files(pyGEMME.data.matrix).joinpath(f'{name}.txt'))

def get_alphabet_filename(name="lz-bl.7"):
	return str(files(pyGEMME.data.alphabets).joinpath(f'{name}.txt'))
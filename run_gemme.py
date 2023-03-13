# -*- coding: utf-8 -*-

# Copyright (c) 2018: Elodie Laine
# This code is part of the gemme package and governed by its license.
# Please see the LICENSE.txt file included as part of this package.

import argparse
from pathlib import Path
from datetime import datetime
from pyGEMME.model import transform_fit_predict


def check_argument_groups(parser, arg_dict, group, argument):
    """
    Check for use of arguments.
    Raise an error if the parser uses an argument that belongs to an argument
    group (i.e. mode) but the group flag is not used or if the argument is
    required in the argument group and it's not used.
    Notes:
    - argument and group should be strings and should start with - or --
    - group should be a flag with action 'store_true'
    - argument shouldn't be a flag, and should be None if it isn't used.
    >>> import argparse
    >>> parser = argparse.ArgumentParser()
    >>> parser.add_argument('--phylo', action='store_true')
    _StoreTrueAction(...)
    >>> parser.add_argument('--inseq')
    _StoreAction(...)
    >>> args = parser.parse_args(['--phylo', '--inseq', 'not_none'])
    >>> check_argument_groups(parser, vars(args), '--phylo', '--inseq', True)
    """
    c = 0
    for myArg in argument:
        arg_name = myArg.replace("-", "")
        if arg_dict[arg_name]!='':
            c = c+1
            group_name = group.replace("-", "")
            if arg_dict[group_name]=="input":
                if c!=1:
                    parser.error("gemme requires " + str(argument) +
                        " if " + group + " is set to input.")
            else:
                if c>0:
                    parser.error("gemme requires " + group +
                " to be set to input if " + str(argument) + " is used.")
    return None

def parse_command_line():
    """
    Parse command line.
    It uses argparse to parse phylosofs' command line arguments and returns the
    argparse parser.
    """
    parser = argparse.ArgumentParser(
        prog="gemme",
        description="""
        GEMME (Global Epistasis Model for predicting Mutational Effects) 
        is a tool to predict mutational outcomes based on sequence analysis
        """,
        epilog="""
        It has been developed at LCQB (Laboratory of Computational and
        Quantitative Biology), UMR 7238 CNRS, Sorbonne Université.
        If you use it, please cite:
        Laine E, Karami Y, Carbone A. Predicting Protein Mutational Landscapes
        using Evolutionary Conservation and Global Epistasis
        """,
    )

    parser.add_argument(
    	'input', metavar='aliFile', type=str,
        help='input alignment file in FASTA format (first sequence = ungapped query sequence)',
    )

    parser.add_argument(
        '-m', '--mutations',
        help='text file containing the mutations of interest. Each line of the file should contain a mutation (e.g. D136R) or combination of mutations separated by commas and ordered according to their positions in the sequence (e.g. D136R,V271A)',
        default=None
    )

    parser.add_argument(
        '-k','--keep_tmp_files',default=False, action='store_true'
    )

    retMet_args = parser.add_argument_group('Conservation levels calculation', """
        Arguments used for computing the conservation levels.
        """)
    
    retMet_args.add_argument(
        '-n', '--nIter',
        help='number of iterations to compute the conservation levels',
        default='1'
    )

    retMet_args.add_argument(
        '-N', '--NSeqs_max',
        help='maximum number of sequences to compute the conservation levels',
        default='20000'
    )
        
    retMet_args.add_argument(
        '-r', '--retrievingMethod',
        help='mode to retrieve related sequences for conservation levels calculation (input, local or server)',
        default='local'
    )

    retMet_args.add_argument(
        '-b', '--blastFile',
        help='psiblast file containing related sequences',
        default=None
    )

    retMet_args.add_argument(
        '-f', '--fastaFile',
        help='fasta file containing related sequences',
        default=None
    )

    args = parser.parse_args()

    arg_dict = vars(args)

    #check_argument_groups(parser, arg_dict, '--retrievingMethod', ['--blastFile','--fastaFile'])
 
    # Check flag arguments
    if args.input is None:
        parser.error("gemme requires an input alignment.")

    return args


def run(input_file,mutations_file,retrieval_method, blast_file=None,fasta_file=None,n_iter=1,N_seqs_max=20000, keep_tmp_files=False):
    """
    Run pipeline
    """

    # make jobname
    timestamp = datetime.now().strftime("%Y_%m_%d-%I_%M_%S_%p")
    jobname = Path(input_file).stem + timestamp
    Path(jobname).mkdir(exist_ok=True)

    # check input alignment mode
    if blast_file is None and fasta_file is None:
        alignment_file = input_file
        mode = 'fasta'
    elif blast_file is None:
        alignment_file = fasta_file
        mode = 'fasta'
    elif fasta_file is None:
        alignment_file = blast_file
        mode = 'blast'
    else:
        raise Exception('BLAST and FASTA input given simultaneously!')

    transform_fit_predict(
        alignment_file,
        output_dir=jobname,
        n_iter=n_iter,
        N_seqs_max=N_seqs_max,
        mode = mode,
        retrieval_method=retrieval_method,
        jet_conf_file='jet.conf',
        keep_tmp_files=keep_tmp_files,
        model_subs_matrix="blosum62p", 
		alphabet_file="pyGEMME/alphabets/lz-bl.7.txt",
        mutations_file=mutations_file
        )


def main():
    args = parse_command_line()
    run(args.input,args.mutations,args.retrievingMethod,args.blastFile,args.fastaFile,args.nIter,args.NSeqs_max,args.keep_tmp_files)

if (__name__ == '__main__'):
    main()


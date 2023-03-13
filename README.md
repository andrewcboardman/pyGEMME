# pyGEMME

pyGEMME is a Python package which bundles together GEMME, JET and muscle3.6 into one easily installable package.

GEMME is implemented in Python and R (https://www.python.org/, https://cran.r-project.org/). See Laine E, Karami Y, Carbone A. GEMME: A Simple and Fast Global Epistatic Model Predicting Mutational Effects,  Molecular Biology and Evolution, Volume 36, Issue 11, November 2019, Pages 2604–2619.

## Installation

To install pyGEMME using conda, download [environment.yml](https://raw.githubusercontent.com/andrewcboardman/pyGEMME/main/environment.yml) and run `conda create -n pyGEMME -f environment.yml` followed by `source setup.sh`. You will also need the Java Runtime Environment; on linux systems, you can install this using `apt install openjdk-9-jre`.

Other dependencies
Joint Evolutionary Trees: http://www.lcqb.upmc.fr/JET2/
seqinr R package: https://cran.r-project.org/web/packages/seqinr/index.html


# Usage notes

The inputAli.fasta is a mandatory argument that corresponds to the input multiple sequence 
alignment file, in FASTA format. The query sequence is taken as the first sequence in the alignment.

By default, GEMME will predict the effect of all possible single mutations at all positions in the 
query sequence. Alternatively, a set of single or multiple mutations can be given with the option -m.
Eachline of the file should contain a mutation (e.g. D136R) or combination of mutations separated 
by commas and ordered according to their positions in the sequence (e.g. D136R,V271A).

GEMME calls JET2 to compute evolutionary conservation levels. By default, JET2 will retrieve a set
of sequences related to the query, independent from the input set, according to specific criteria.
The retrieval method used in JET2 is PSI-BLAST, which can perform the search either locally (by 
default) or remotely (-r server). Alternatively, the user can provide her/his own psiblast file  
(-r input-b pFile) or her/his own multiple sequence alignment in FASTA format (-r input -f fFile).
JET is run in its iterative mode, iJET, 10 times and the final conservation levels are the maxium 
values obtained over the 10 iterations. 
JET2 configuration file is: default.conf.
JET2 output file is: myProt_jet.res.

By default, GEMME will output mutational effects predictions obtained from the global epistatic model,
the independent model, and a combination of those two using a reduced alphabet (alphabets/lw-i.11.txt):
myProt_pred_evolEpi.txt
myProt_normPred_evolEpi.txt
myProt_pred_evolInd.txt
myProt_normPred_evolInd.txt
myProt_normPred_evolCombi.txt
The values of interest are the normalized predictions (normPred). Each file contains a 20 x n matrix, 
where n is the number of positions in the query sequence.
If the user provides her/his own list of mutations, then only the global epistatic model will be run 
and the output file will contain 2 columns, the first one with the mutations, the second one with the 
normalized predicted effects.



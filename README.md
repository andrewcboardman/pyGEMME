# pyGEMME

## About

pyGEMME is a Python reimplementation of the [GEMME](http://www.lcqb.upmc.fr/GEMME/Home.html) package which is lightweight and built in Numpy and Biopython. It is bundled with a notebook for querying the MMSeqs server, using code from the ColabFold package.

## Installation & usage

The relevant packages can be installed using the Conda package manager by running `conda install -n pyGEMME --file environment.yml`.

To test the full pyGEMME on the provided MSA (BLAT_ECOLX), run the command `python run_gemme.py -t -p`.

To run pyGEMME on a new MSA, first generate the MSA (e.g. from the MMSeqs server using `mmseqs_qeury.ipynb`) and then run the command `python run_gemme.py -i MSA_FILE -o OUTPUT_DIR -p`, replacing MSA_FILE and OUTPUT_DIR with the locations of the desired input MSA file and output directories respectively;
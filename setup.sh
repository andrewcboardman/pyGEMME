CRAN_MIRROR="https://cran.ma.imperial.ac.uk/"

conda env create -f environment.yml
Rscript -e "install.packages('seqinr', repos='$CRAN_MIRROR')"
Rscript -e "install.packages('RColorBrewer', repos='$CRAN_MIRROR')"

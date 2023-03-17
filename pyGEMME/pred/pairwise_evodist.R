
library("seqinr")

# amino acid full alphabet
aa = c("a","c","d","e","f","g","h","i","k","l","m","n","p","q","r","s","t","v","w","y")


alignment_file = "tests/BLAT_model/BLAT_A.fasta"
jet_file = "tests/BLAT_model/BLAT_jet.res"
output_file = "tests/BLAT_model/model/BLAT_pairwise_evolDist.mat"


# read alignment and convert to matrix
ali=as.matrix.alignment(read.alignment(alignment_file,format="fasta"))
# convert all non-amino acids to "-"
ali[!ali%in%aa] = "-"
N = dim(ali)[[1]]
npos = dim(ali)[[2]]

# read evolutionary traces computed by JET
jet=read.table(jet_file,head=TRUE)
if(sum(colnames(jet)=="traceMax")==1){trace=jet[,"traceMax"]}else{trace=jet[,"trace"]}

npos = 10
distTrace = matrix(0, N, N)
encoded_cols = array(0, c(20, npos, N))
for (i in 1:npos) {
    for (j in 1:20) {
        # convert to binary matrix with respect to the query (0: identical, 1: different)
        encoded_cols[i, j ,] = ifelse(ali[,i] == aa[j], 0, 1)
    }
# dist_trace_ = encoded_col %*% t(encoded_col) * trace[i]^2
# distTrace = distTrace + dist_trace_
}
dim(encoded_cols) <- c(npos, 20 * N)
dist_trace_ = encoded_cols %*% t(encoded_cols) * trace[i]^2
distTrace = distTrace + dist_trace_

#write.table(distTrace, output_file)
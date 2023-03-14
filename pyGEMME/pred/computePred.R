# Copyright (c) 2018: Elodie Laine
# This code is part of the gemme package and governed by its license.
# Please see the LICENSE.txt file included as part of this package.

# Usage: Rscript --save runPred.R XXXX

library("seqinr")

# amino acid full alphabet
aa = c("a","c","d","e","f","g","h","i","k","l","m","n","p","q","r","s","t","v","w","y")



# get arguments 
args = commandArgs(trailingOnly = TRUE)
prot = args[1]
wt_sequence_file = args[2]
alignment_file = args[3]
jet_file = args[4]
subs_matrix_file = args[5]
reduced_alphabet_file = args[6]
output_dir = args[7]
function_path = args[8]
mutfile_does_not_exist = args[9]
mutfile = args[10]

print(function_path)
source(function_path)

alpha = 0.6


##################Read inputs

# read alignment and convert to matrix
ali=as.matrix.alignment(read.alignment(alignment_file,format="fasta"))
ali[!ali%in%aa] = "-"
npos = dim(ali)[[2]]
# convert to binary matrix with respect to the query (0: identical, 1: different)
binAli = t(apply(ali,1,f<-function(x){return(x!=ali[1,])}))
pId =apply(binAli,1,f<-function(x){return(1-sum(x)/npos)})
# Compute conservation
aliCons = ali[pId>0.6,]
aliVeryCons = ali[pId>0.8,]
N = c(dim(ali)[[1]],dim(aliCons)[[1]],dim(aliVeryCons)[[1]]) # number of sequences
#indicesQ = getIndicesAli(ali$seq,1)


# BLOSUM62  matrix
subs_matrix = as.matrix(read.table(subs_matrix_file),row.names=1)
rownames(subs_matrix) = tolower(rownames(subs_matrix))
colnames(subs_matrix) = tolower(colnames(subs_matrix))
subs_matrix = subs_matrix[aa, aa]
subs_matrix_rowsums = apply(subs_matrix,1,sum)

# Alphabet
alphabet = read.table(reduced_alphabet_file,sep="")$V1

# read evolutionary traces computed by JET
jet=read.table(jet_file,head=TRUE)
if(sum(colnames(jet)=="traceMax")==1){trace=jet[,"traceMax"]}else{trace=jet[,"trace"]}
#traceAli = sweep(binAli, MARGIN=2, trace, `*`)
# compute evolutionary distances of all sequences with respect to the query
distTrace = binAli[2:N[1],] %*% trace^2

# Read wild type sequence
wt=read.fasta(wt_sequence_file)
n = length(wt[[1]])
wt = wt[[1]][1:n]

# Set working directory before writing output
setwd(output_dir)

###################


#resAliCons = computePSSM(aliCons)
res = list(
  computePSSM(ali,N[1], npos, subs_matrix_rowsums, subs_matrix, aa),
  computePSSM(aliCons,N[2], npos, subs_matrix_rowsums, subs_matrix, aa),
  computePSSM(aliVeryCons,N[3], npos, subs_matrix_rowsums, subs_matrix, aa)
  )
write.table(res[[1]][[3]],paste0(prot,"_pssm.txt"))
write.table(res[[2]][[3]],paste0(prot,"_pssm60.txt"))
write.table(res[[3]][[3]],paste0(prot,"_pssm80.txt"))

##### Modelling

# compute sequence counts for each position and each aa
nbSeqs = computeNbSeqs(ali)
# compute gap counts for each position
nbGaps = N[1] - apply(nbSeqs,2,sum)
# compute log-odd ratios between mutated and wt sequence counts
predInd = computePredNbSeqs(wt,nbSeqs,aa)
rownames(predInd)=aa

predEpi=computePredSimple(ali,distTrace,wt,5)
rownames(predEpi)=aa
evolDist = predEpi/sum(trace^2)
evolDist[is.na(evolDist)] = 1



# output the conservation values
conservation = rbind(trace,KL=res[[1]][[2]],SE=res[[1]][[1]],gap=nbGaps/N[1],KL60=res[[2]][[2]],SE60=res[[2]][[1]],KL80=res[[3]][[2]],SE80=res[[3]][[1]])


if(mutfile_does_not_exist) {
  
  print("running independent model...")
  normPredInd = normalizePredWithNbSeqsZero(predInd,trace,wt,aa)
  rownames(normPredInd)=aa

  print("running global epistatic model...")
  normPred=normalizePred(predEpi, trace, wt,aa)
  rownames(normPred)=aa
  
  print("running combined model...")
  normPredCombi = normalizePredWithNbSeqsPC(predEpi,trace,wt,alpha,nbSeqs,alphabet,aa)
  rownames(normPredCombi)=aa

  
  print("done")
} else { 

  print("reading selected mutations...")
  res=readMut(mutfile)
  pos = res[[1]]
  subsaa = res[[2]]
  rawMut = res[[3]]
  print("done")

  normPredInd = normalizePredWithNbSeqsZeroSelMult(predInd, trace, wt, list(pos,subsaa))
  names(normPredInd)=rawMut

  normPred=normalizePredSelMult(predEpi, trace, wt, list(pos,subsaa))
  names(normPred)=rawMut

  normPredCombi = normalizePredWithNbSeqsPCSelMult(predEpi,trace,wt,list(pos,subsaa),alpha,nbSeqs,alphabet,aa)
  names(normPredCombi)=rawMut
}


## Write outputs to file

write.table(conservation,paste0(prot,"_conservation.txt"))
# output the sequence counts log-odd ratios
write.table(predInd,paste0(prot,"_pred_evolInd.txt"))
# output the predicted mutational effects based on sequence counts (conservation at the bottom)
write.table(normPredInd,paste0(prot,"_normPred_evolInd.txt"))

# output the evolutionary distances between the query and the closest variants
write.table(evolDist,paste0(prot,"_pred_evolEpi.txt"))
# output the predicted mutational effects based on evolutionary distance (conservation at the bottom)
write.table(normPred,paste0(prot,"_normPred_evolEpi.txt"))
# output the predicted mutational effects based on sequence counts (conservation at the bottom)
write.table(normPredCombi,paste0(prot,"_normPred_evolCombi.txt"))


# Generate plots
if(mutfile_does_not_exist){
  print("generating the plots...")
  plotMatBlues(paste(prot,"_normPred_evolEpi",sep=""))
  plotMatGreens(paste(prot,"_normPred_evolInd",sep=""))
  plotMatOranges(paste(prot,"_normPred_evolCombi",sep=""))
  print("done")
}


	

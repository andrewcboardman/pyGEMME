# Copyright (c) 2018: Elodie Laine
# This code is part of the gemme package and governed by its license.
# Please see the LICENSE.txt file included as part of this package.

# module containing all R functions to manipulate data and compute pedictions


###########################################################################################
###################### Functions to manipulate alignment data #############################
###########################################################################################

# Given an alignment seqAli, the number of the sequence of interest (i) 
# and a position (pos) in the sequence, get the corresponding 
# index/position in the alignment (used by countSeq, see below)
getIndexMatAli<-function(matAli, i, pos){
        seq = matAli[i,]
        countAll = 0
        count = 0
        while(count<pos){
                countAll = countAll + 1
                if(seq[countAll]!="-"){count = count + 1}
        }
        return(countAll)
}

# Given an alignment seqAli, the number of the sequence of interest (i) 
# and a vector of positions (pos) in the sequence, get the corresponding 
# indices/positions in the alignment (used by countSeq, see below)
getIndexAli<-function(seqAli, i, pos){
  seq = strsplit(seqAli[i][[1]],"")[[1]]
  res = c()
  countAll = 0
  count = 0
  for(p in pos){
    while(count<p){
      countAll = countAll + 1
      if(seq[countAll]!="-"){count = count + 1}
    }
    res = c(res,countAll)
  }
  return(res)
}

# Given an alignment seqAli and the number of the sequence of interest (i) 
# get all the indices/positions in the alignment 
getIndicesAli<-function(seqAli, i){
  seq = strsplit(seqAli[i][[1]],"")[[1]]
  n = sum(seq!="-"&seq!=".")
  res = c()
  countAll = 0
  count = 0
  for(p in 1:n){
    while(count<p){
      countAll = countAll + 1
      if(seq[countAll]!="-"&seq[countAll]!="."){count = count + 1}
    }
    res = c(res,countAll)
  }
  return(res)
}

# Compute the weights for sequences in the alignment
computeWeights<-function(matAli,theta=0.2){
  nbSeq = dim(matAli)[[1]]
  nbPos = dim(matAli)[[2]]
  res=c()
  matD=matrix(0,nr=nbSeq,nc=nbSeq)
  for(i in 1:(nbSeq-1)){
    for(j in (i+1):nbSeq){
      matD[i,j]=sum(matAli[i,]!=matAli[j,])/nbPos
      matD[j,i]=matD[i,j]}
    val=1/(sum(matD[i,]<theta))
    res=c(res,val)
    #print(val)
  }
  return(res)
}

computeNbSeqs<-function(mat,gap=FALSE){
  if(gap){aa=c(aa,"-")}
  occ = apply(mat,2,table)
  n = length(occ)
  res = matrix(0,nc=n,nr=length(aa))
  colnames(res) = 1:n
  rownames(res) = aa
  for(i in 1:n){
    reducedAA = intersect(aa,names(occ[[i]]))
    res[reducedAA,i] = occ[[i]][reducedAA]
  }
  return(res)
}

# compute variability levels as the Shanon entropy
computeConsSE<-function(mat,N){
  mat[mat==0] = 1
  return(apply(mat/N,2,function(x){return(-sum(x*log2(x)))}))
}

# compute conservation levels as the Kullback-Leibler divergence
computeConsKL<-function(mat, blosum62){
  n = apply(mat,2,sum)
  mat = t(t(mat)/n)
  mat = mat + 0.000001
  bg = blosum62[,1][order(rownames(blosum62))]
  return(apply(mat,2,function(x){return(sum(x*log2(x/bg)))}))
}

computeSeqWeights<-function(mat,N,npos){
  # compute the number of different observed amino acids in each column
  r = apply(mat,2,function(x){return(length(unique(x)))})
  # compute occurrence of each letter at each position
  occMat = computeNbSeqs(mat,TRUE)
  #print(occMat)
  # compute weights for each sequence
  weightMat=matrix(nr=N,nc=npos)
  for(k in 1:N){
    midx=cbind(mat[k,],1:npos)
    weightMat[k,] = 1 / (occMat[midx] * r)
  }
  indObs = sum(r)/npos
  return(list(weightMat,indObs))
}

computePseudoCounts<-function(freqmat,npos, aa, bg, bgp){
  PC = matrix(nr=length(aa),nc=npos)
  print(dim(PC))
  rownames(PC)= aa
  
  print(dim(PC))
  print(dim(bgp))
  for(a in aa){
    PC[a,] = apply(freqmat,2,function(x){return(sum(x/bg*bgp[a,]))})
  }
  return(PC)
}

computePSSM<-function(mat,N,npos, bg, bgp, aa){
  if(is.null(dim(mat))){return(list(0,0,0))}
  # compute sequence weights
  res  = computeSeqWeights(mat,N,npos)
  weights = res[[1]]
  indObs = res[[2]]
  # extend amino acid alphabet with gaps
  aa_gap=c(aa,"-")
  # compute  weighted occurrences
  occMat = matrix(0,nr=length(aa_gap),nc=npos)
  rownames(occMat) = aa_gap
  for(i in 1:npos){
    counts = tapply(weights[,i],mat[,i],sum)
    occMat[names(counts),i] = counts
  }
  # distribute gaps
  occMat = occMat[1:20,] + t(occMat["-",]%*%t(bg))
  # compute pseudo-counts
  PC = computePseudoCounts(occMat,npos, aa, bg, bgp)
  # distribute pseudo-counts according to an empirical (?) weight
  occMat = (occMat*(indObs-1) + PC) /indObs
  # normalize to relative frequencies
  n = apply(occMat,2,sum)
  freq = t(t(occMat)/n)
  # compute Shannon entropy
  SE = apply(freq,2,function(x){return(-sum(x*log2(x)))})
  # compute Kullback-Leibler divergence
  KL = apply(freq,2,function(x){return(sum(x*log2(x/bg)))})
  # divide by bg freqs and convertto log-scores
  pssm = 2 * log2(freq/bg)
  return(list(SE,KL,pssm))
}

computePSSM2<-function(mat,blosum62){
  occMat = computeNbSeqs(mat)
  n = apply(occMat,2,sum)
  pssm = t(t(occMat)/n)
  pssm = pssm + 0.000001
  bg = blosum62[,1][order(rownames(blosum62))]
  res = 2 * log2(pssm/bg)
  return(list(pssm,res))
}

computeNbSeqsAlph<-function(nbSeqs,AA, aa){
  facAA = rep(NA,20)
  for(Let in AA){
          Let = toString(Let)
          lets = strsplit(Let,"")[[1]]
          for(let in lets){
                  facAA[aa==tolower(let)]=Let
          }
  }
  newNbSeqs = apply(nbSeqs,2,function(x){return(tapply(x,facAA,sum))})
  return(newNbSeqs)
}

readMut<-function(fname){
	rawMut = read.table(fname,colClass="character")$V1
	n = length(rawMut)
	pos = list()
	subsaa = list()
	for(i in 1:n){
		mut = strsplit(rawMut[i],",")[[1]]
		tmpPos = c()
		tmpSubsaa = c()
		for(m in mut){
			le=nchar(m)
			tmpPos = c(tmpPos,as.numeric(substr(m,2,le-1)))
			tmpSubsaa = c(tmpSubsaa,substr(m,le,le))
		}
		pos[[i]] = tmpPos
                subsaa[[i]] = tmpSubsaa
	}
	return(list(pos,subsaa,rawMut))
}

which.class<-function(a,AA){
  # splits up alphabet
        for(Let in AA){
                Let = toString(Let)
                lets = tolower(strsplit(Let,"")[[1]])
                if(sum(lets==a)>0){return(Let)}
        }
}

###########################################################################################
###################### Functions to compute predictions  ##################################
###########################################################################################

# log-odd ratio between the mutated and wild-type sequence counts
computePredNbSeqs<-function(wt, nbseqs,aa){
  n = length(wt)
  pred = matrix(nc=n,nr=20)
  rownames(pred)=aa
  for(i in 1:n){
    for(a in aa){
       if(a!=wt[i]){
         pred[a,i] = log(max(1,nbseqs[a,i])/nbseqs[wt[i],i])}
      else{
        pred[a,i] = 0}
    }
  }
  return(pred)
}

computePredNbSeqsSelMult<-function(wt, nbseqs, listMut){
  # check mutation list
  pos = listMut[[1]]
  let = listMut[[2]]
  N = length(pos)
  if(N!=length(let)){print("Problem!!!! not the same number of positions and aas !!")}
  # compute predictions
  pred = list()
  for(i in 1:N){
    myPos = pos[[i]]
    myAA = tolower(let[[i]])
    n = length(myPos)
    predTMP = c()
    # treat individual mutations
    for(k in 1:n){
    #  subsStr=paste(subsStr,paste(j[k],myAA[k],sep=""),",")
      predTMP = c(predTMP,log(max(1,nbseqs[myAA[k],myPos[k]])/nbseqs[wt[myPos[k]],myPos[k]]))
    }
    pred[[i]]=predTMP
  }
  return(pred)
}

normalizePredWithNbSeqsZero<-function(pred, trace, wt,aa){
  n = length(trace)
  normPred = matrix(nc=n,nr=20)
  rownames(normPred) = aa
  for(i in 1:n){
    normPred[,i] = pred[,i] * trace[i]
    normPred[aa==wt[i],i] = NA
  }
  return(normPred)
}

normalizePredWithNbSeqsZeroSelMult<-function(pred, trace, wt, listMut){
  pos = listMut[[1]]
  let = listMut[[2]]
  N = length(pos)
  normPred = rep(NA,N)
  for(i in 1:N){
    myPos = pos[[i]]
    myAA = tolower(let[[i]])
    n = length(myPos)
    predTMP = c()
    for(k in 1:n){predTMP=c(predTMP,pred[myAA[k],myPos[k]])}
    normPred[i] = sum(predTMP * trace[myPos])
  }
  return(normPred)
}

normalizePredWithNbSeqs<-function(pred, trace, wt, nbseqs, alpha, aa){
  n = length(trace)
  normPred = matrix(nc=n,nr=20)
  rownames(normPred) = aa
  maxVal = max(apply(nbseqs,2,sum))
  for(i in 1:n){
    normPred[,i] = pred[,i] / max(pred[!is.na(pred)]) * log(1/maxVal) * (-1)
    normPred[is.na(normPred[,i]),i] = - log(1/maxVal)
    for(a in aa){
      normPred[a,i] = alpha * normPred[a,i] - (1-alpha) * log(max(1,nbseqs[a,i])/nbseqs[wt[i],i])
    }
    normPred[,i] = normPred[,i] * trace[i]
    normPred[aa==wt[i],i] = NA
  }
  return(-normPred)
}

findMinDist<-function(vec,delta){
  le = length(vec)
    # if no sequence was found with the substituting aa
    # then no predicted value is computed  
    if(le==0){res = NA}
    else{
    	# if only one sequence was found
    	# then just take this distance
    	if(le==1){res = NA}
    	else{
    		# if the minimum distance is null
            # we are in the presence of the wild type
            if(vec[1]==0){res = vec[1]}
            	# otherwise, it means we have at least two non-zero values
                else{
                	# if there is no significant gap between the first and second values 
                	# take the first value 
                	if((vec[2]-vec[1])<=delta){res = vec[1]}
                    	# otherwise, take the second one
                	else{res = vec[2]}
                }
            }
        }
    return(res)
}

# Compute all evolutionary distances with all homologous 
# sequences in an alignment (mat) with respect to all positions
# and all possible amino acid substitutions
# returns a matrix with aa rows and nbPos columns
# for each aa and each position, we retain only the smallest distance
computePredSimple<-function(mat, distTrace, wt, thresh){
  # number of positions
  leSeq = length(wt)
  # result matrix
  pred = matrix(nc=leSeq,nr=20)
  # count nb seqs
  #N = dim(mat)[[2]]
  # the resulting matrix contains trace values
  # only in mutated sequences and at relevant positions
 # matAli = convertAliToTrace(mat,trace,wt)
  # go over all positions
  for(i in 1:leSeq){
    res = vector()
    for(a in aa){
      if(a!=wt[i]){
        sel = which(mat[,i]==a)
        if(length(sel)>0){sortedDist=sort(distTrace[sel-1])}
        else{sortedDist = c()}
        res = c(res,findMinDist(sortedDist,thresh))
      }
      else{
        res = c(res,0)
      }
    }
    pred[,i] = res
  }
  return(pred)
}

# Normalize predicted values 
# scale between 0 and 2, with 2 replacing NAs, 
# and then multiply by the trace of each position
# if we have not found any sequence with the substitution of interest
# we assume the effect is maximal (2)
# then the effects are scaled according to the trace of the positions
# the more conserved a position, the more important the effect
normalizePred<-function(pred, trace, wt, aa){
	n = length(trace)
	normPred = matrix(nc=n,nr=20)
	for(i in 1:n){
		normPred[,i] = pred[,i] / max(pred[!is.na(pred)])
		normPred[is.na(normPred[,i]),i] = 1
		normPred[,i] = normPred[,i] * trace[i]
		normPred[aa==wt[i],i] = NA
	}
	return(-normPred)
}

normalizePredWithNbSeqsPC<-function(pred, trace, wt, alpha, nbSeqs, alphabet,aa){
	n = length(trace)
	normPred = matrix(nc=n,nr=20)
	rownames(normPred) = aa
	nbseqs = computeNbSeqsAlph(nbSeqs,alphabet,aa)
	maxRefVal = max(apply(nbseqs,2,sum))
	# for each position in the sequence
	for(i in 1:n){
		normPred[,i] = pred[,i] / max(pred[!is.na(pred)]) * log(1/maxRefVal) * (-1)
		normPred[is.na(normPred[,i]),i] = - log(1/maxRefVal)
		# for each possible substitution
		for(a in aa){
			A = which.class(wt[i],alphabet)
			B = which.class(a,alphabet)
			normPred[a,i] = alpha * normPred[a,i] - (1-alpha) * log(max(1,nbseqs[B,i])/nbseqs[A,i])
		}
		normPred[,i] = normPred[,i] * trace[i]
		normPred[aa==wt[i],i] = NA
	}
	return(-normPred)
}

normalizePredWithNbSeqsPCSelMult<-function(pred, trace, wt, listMut, alpha, nbSeqs, alphabet,aa){
  # check mutation list
  pos = listMut[[1]]
  let = listMut[[2]]
  N = length(pos)
  if(N!=length(let)){print("Problem!!!! not the same number of positions and aas !!")}
  normPred = rep(NA,N)
  allVal = unlist(pred)
  maxVal = max(allVal[!is.na(allVal)])
  nbseqs = computeNbSeqsAlph(nbSeqs,alphabet,aa)
  maxRefVal = max(apply(nbseqs,2,sum))
  # for each position in the sequence
  for(i in 1:N){
    predTMP=c()
    myPos = pos[[i]]
    myAA = tolower(let[[i]])
    n = length(myPos)
    for(k in 1:n){
	    A = which.class(wt[myPos[k]],alphabet)
	    B = which.class(myAA[k],alphabet)
	    ratio = log(max(1,nbseqs[B,myPos[k]])/nbseqs[A,myPos[k]])
	    val = pred[myAA[k],myPos[k]] / maxVal * log(1/maxRefVal) * (-1)
	    if(is.na(val)){val=- log(1/maxRefVal)}
	    predTMP=c(predTMP, alpha * val - (1-alpha) * ratio)
    }
    normPred[i] = sum(predTMP * trace[myPos])
  }
  return(-normPred)
}

normalizePredSelMult<-function(pred, trace, wt, listMut){
	pos = listMut[[1]]
 	let = listMut[[2]]
  	N = length(pos)
	allVal = unlist(pred)
	maxVal = max(allVal[!is.na(allVal)])
	normPred=rep(NA,N)
	for(i in 1:N){
		myPos = pos[[i]]
    		myAA = tolower(let[[i]])
    		n = length(myPos)
    		predTMP = c()
    		for(k in 1:n){
			val = pred[myAA[k],myPos[k]] / maxVal
			if(is.na(val)){val = 1}
			predTMP=c(predTMP,val)}
		normPred[i] = sum(predTMP * trace[myPos])
	}
	return(-normPred)
}


library("RColorBrewer")

##############################################
### functions for ploting the results
##############################################
plotMatOranges<-function(pred_name){
  
	pred = as.matrix(read.table(paste(pred_name, ".txt", sep="")))
	sel=seq(1,dim(pred)[[2]])
        Iwidth = round(sqrt(dim(pred)[[2]])) + 3
	prop=99
	pred[is.na(pred)]=0
	minVal = floor(min(pred))
	maxVal = ceiling(max(pred))
	p = (100-prop)/100
	cutVal = quantile(pred,c(0,p))[2]
	jpeg(paste(pred_name, '.jpg', sep=""), width = Iwidth, height = 5, units = 'in', res = 300)
	print(c(-maxVal,-cutVal,-minVal))
	image(sel,1:20,-t(pred[20:1,sel]),
    breaks=c(seq(-maxVal,-cutVal,le=80),-minVal),
    col=c(
      colorRampPalette(brewer.pal(9,"Greys"))(80)[1:10],
      colorRampPalette(brewer.pal(9,"Oranges"))(80)[11:70],
      colorRampPalette(brewer.pal(9,"Greys"))(80)[61:70]
      ),
      yaxt="n",xaxt="n",ylab="",xlab="",bty="n")
	axis(1,c(sel[1]-1,sel[seq(0,length(sel),by=10)]),las=2)
	axis(2,1:20,toupper(rownames(pred)[20:1]),las=2,tick=FALSE)
	dev.off()
}

plotMatGreens<-function(pred_name){

	pred = as.matrix(read.table(paste(pred_name, ".txt", sep="")))
        Iwidth = round(sqrt(dim(pred)[[2]])) + 3
	sel=seq(1,dim(pred)[[2]])
	prop=99
        pred[is.na(pred)]=0
        minVal = floor(min(pred))
        maxVal = ceiling(max(pred))
        p = (100-prop)/100
        cutVal = quantile(pred,c(0,p))[2]
        jpeg(paste(pred_name, '.jpg', sep=""), width = Iwidth, height = 5, units = 'in', res = 300)
        print(c(-maxVal,-cutVal,-minVal))
        image(sel,1:20,-t(pred[20:1,sel]),breaks=c(seq(-maxVal,-cutVal,le=80),-minVal),
          col=c(colorRampPalette(
            brewer.pal(9,"Greys"))(80)[1:10],colorRampPalette(brewer.pal(9,"Greens"))(80)[11:70],colorRampPalette(brewer.pal(9,"Greys"))(80)[61:70]),yaxt="n",xaxt="n",ylab="",xlab="",bty="n")
        axis(1,c(sel[1]-1,sel[seq(0,length(sel),by=10)]),las=2)
        axis(2,1:20,toupper(rownames(pred)[20:1]),las=2,tick=FALSE)
        dev.off()
}

plotMatBlues<-function(pred_name){
  
	pred = as.matrix(read.table(paste(pred_name, ".txt", sep="")))
	sel=seq(1,dim(pred)[[2]])
        Iwidth = round(sqrt(dim(pred)[[2]])) + 3
	prop=99 #95
	pred[is.na(pred)]=0
	minVal = floor(min(pred))
	maxVal = ceiling(max(pred))
	p = (100-prop)/100
	cutVal = quantile(pred,seq(0,p,by=p))[2]
	print(c(-maxVal,-minVal,-cutVal))
	jpeg(paste(pred_name, '.jpg', sep=""), width = Iwidth, height = 5, units = 'in', res = 300)
	#image(sel,1:20,-t(pred[,sel]),breaks=c(seq(-maxVal,-cutVal,le=80),-minVal),col=c(colorRampPalette(brewer.pal(9,"Greys"))(80)[11:20],colorRampPalette(brewer.pal(9,"Blues"))(80)[11:70],colorRampPalette(brewer.pal(9,"Greys"))(80)[71:80]),yaxt="n",xaxt="n",ylab="",xlab="",bty="n")
	image(sel,1:20,-t(pred[20:1,sel]),breaks=c(seq(-maxVal,-cutVal,le=80),-minVal),
    col=c(
      colorRampPalette(
        brewer.pal(9,"Greys"))(80)[1:10],
      colorRampPalette(brewer.pal(9,"Blues"))(80)[11:70],colorRampPalette(brewer.pal(9,"Greys"))(80)[61:70]),yaxt="n",xaxt="n",ylab="",xlab="",bty="n")
	axis(1,c(sel[1]-1,sel[seq(0,length(sel),by=10)]),las=2)
	axis(2,1:20,toupper(rownames(pred)[20:1]),las=2,tick=FALSE)
	dev.off()
}

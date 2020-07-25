rm(list=ls())
library(gasper)
library(ggplot2)
source("utils.R")

if (as.numeric(R.version$minor)>6){
  RNGkind(sample.kind = "Rounding")
}

# set the random seed for reproductibility
set.seed(0)
A <- full(pittsburgh$sA)
n <- nrow(A)
L <- laplacian_mat(A)
val1 <- eigendec(L)
evalues <- val1$evalues
evectors <- val1$evectors
lmax <- max(evalues)
b <- 2
kmax <- floor(log(lmax)/log(b)) + 2
tf <- tight_frame(evalues, evectors)
diagWWt <- colSums(t(tf)^2)

W <- tf
diagWWt <- colSums(t(W)^2)
WWt<-W%*%t(W)

load("pittsburgh/pittsburghGEO.Rdata")
f <- pittsburghGEO$f1

sigma <- 0.004
sumsqf<-sum(f^2)

# Choose the threshold function, soft=1, JS=2
beta<-2

start_time <- Sys.time()
# Monte-Carlo
MC <- 25
cnt<-1
# Experiment on different number of blocks per scales
nBlockVal<-40:55
Rge<-matrix(rep(0,3*length(nBlockVal)*MC),ncol=3)
res<-data.frame(matrix(rep(0,MC*length(nBlockVal)*3*6),ncol=6))
colnames(res)<-c("MCit","nBlock","global-snr","block-snr","thresh_opt","hatsigma")
for(it in 1:MC){
  noise <- rnorm(n, sd = sigma)
  y <- f + noise
  wcn <- analysis(y,tf)
  wcf <- analysis(f,tf)
  hatsigma2 <- 0.5*sum(A*outer(as.vector(y),as.vector(y),'-')^2)/sum(diag(L))
  hatsigma<-sqrt(hatsigma2)

  # Compute coordinatewise with uniform threshold to compare with
  # comparaison coordinatewise global
  nThresh2<-100
  m1 <- max(abs(wcn)) # max tresh global
  Threshval2<-seq(0,m1,length.out=nThresh2)
  sure<-SURE_MSEthresh(wcn, 
                    wcf,
                    Threshval2, 
                    diagWWt,
                    beta, 
                    sigma,
                    sqrt(hatsigma2),
                    policy = "uniform")
  fOracle<-synthesis(sure$wc[,which.min(sure$res$MSE)],tf)
  fSURE<-synthesis(sure$wc[,which.min(sure$res$SURE)],tf)
  fhatSURE<-synthesis(sure$wc[,which.min(sure$res$hatSURE)],tf)

  for(nb in nBlockVal){
    # iteration and nbBlock
    res[cnt,1]<-it
    res[cnt,2]<-nb
    res[cnt+1,1]<-it
    res[cnt+1,2]<-nb
    res[cnt+2,1]<-it
    res[cnt+2,2]<-nb
    # coordinatewise global
    res[cnt,3]<-SNR(f,fOracle)
    res[cnt,6]<-SNR(f,y)#hatsigma
    res[cnt+1,3]<-SNR(f,fSURE)
    res[cnt+1,6]<-SNR(f,y)#hatsigma
    res[cnt+2,3]<-SNR(f,fhatSURE)
    res[cnt+2,6]<-SNR(f,y)#hatsigma

    # Comment/Uncomment for contiguous or geometric blocks
    # construction de nBlock blocs contigus

    block<-cut(1:n,nb,labels=FALSE)
    shift<-nb*rep(0:kmax,each=n)
    blockScales<-rep(block,kmax+1)+shift

    # construction nBlockX*nBlockY*(J+1) blocs pavant RxR
    # nBlockX<-nb
    # nBlockY<-nb
    # nBlock<-nBlockX*nBlockY
    # blockX<-cut(x1,nBlockX,labels=F)
    # blockY<-cut(x2,nBlockY,labels=F)
    # block<-nBlockX*(blockY-1)+blockX
    # shift<-nBlockX*nBlockY*rep(0:kmax,each=n)
    # blockScales<-rep(block,kmax+1)+shift

    # Compute Oracle MSE, SURE and hatSURE with grid search
    nThresh<-100
    m<-0.05 # maximal threshold
    Threshval<-seq(0,m,length.out=nThresh)
    partialMSE<-partialER<-partialDof<-matrix(rep(0,nb*(kmax+1)*nThresh),ncol=nThresh)
    for(i in 1:max(blockScales)){
      if(any(blockScales==i)){
        sumi<-sum(diagWWt[blockScales==i])
        sumij<-t(wcn[blockScales==i,1])%*%WWt[blockScales==i,blockScales==i]%*%wcn[blockScales==i,1]
        normBlock<-sqrt(sum(wcn[blockScales==i,1]^2))
        for(t in 1:nThresh){
          twc<-blockthresh(wcn[blockScales==i,1],Threshval[t],beta)
          partialMSE[i,t]<-sum((wcf[blockScales==i,1]-twc)^2)
          partialER[i,t]<-min(1,Threshval[t]^beta/normBlock^beta)^2*normBlock^2
          partialDof[i,t]<-(normBlock>=Threshval[t])*((1-Threshval[t]^beta/normBlock^beta)*sumi+beta*Threshval[t]^beta/normBlock^(beta+2)*sumij)
        }
      }
    }

    MSE<-colSums(partialMSE)
    ER<-colSums(partialER)
    Dof<-colSums(partialDof)
    SURE<-ER+2*sigma^2*Dof-n*sigma^2
    hatSURE<-ER+2*hatsigma2*Dof-n*hatsigma2

    # minimizers
    minMSE<-Rge[it,1]<-which.min(MSE)
    res[cnt,5]<-minMSE*m/nThresh
    minSURE<-Rge[it,2]<-which.min(SURE)
    res[cnt+1,5]<-minSURE*m/nThresh
    minhatSURE<-Rge[it,3]<-which.min(hatSURE)
    res[cnt+2,5]<-minhatSURE*m/nThresh

    # best denoised for fixed nBlock
    wcdenoised<-matrix(rep(0,3*n*(kmax+1)),ncol=3)
    for(i in 1:max(blockScales)){
      if(any(blockScales==i)) {
        wcdenoised[blockScales==i,1]<-blockthresh(wcn[blockScales==i,1],Threshval[minMSE],beta)
        wcdenoised[blockScales==i,2]<-blockthresh(wcn[blockScales==i,1],Threshval[minSURE],beta)
        wcdenoised[blockScales==i,3]<-blockthresh(wcn[blockScales==i,1],Threshval[minhatSURE],beta)
      }
    }

    res[cnt,4]<-SNR(f,synthesis(wcdenoised[,1],tf))
    res[cnt+1,4]<-SNR(f,synthesis(wcdenoised[,2],tf))
    res[cnt+2,4]<-SNR(f,synthesis(wcdenoised[,3],tf))
    cnt<-cnt+3
  }
  print(it)
}
end_time<-Sys.time()

ind <- seq(1,length(nBlockVal)*3*MC,
           length(nBlockVal)*3)
SNRin <- mean(res$hatsigma[ind])
SNRinsd <- sd(res$hatsigma[ind])
print(paste0("Mean Input SNR=",round(SNRin,2),
             "+-",round(SNRinsd,2),"dB"))
resblock <- res$`block-snr`[seq(1,3*length(nBlockVal)*MC,3)]
dim(resblock) <- c(length(nBlockVal),MC)
resblock <- t(resblock)
meanblockperf <- apply(resblock,2,mean)
sdblockperf <- apply(resblock,2,sd)
meancorperf <- mean(res$`global-snr`[ind])
sdcorperf <- sd(res$`global-snr`[ind])
print(paste0("Mean best Block SNR=",round(max(meanblockperf),2),
             "+-",round(max(sdblockperf),2),"dB",
             " (optimal mean block L=",nBlockVal[which.max(meanblockperf)],")"))
print(paste0("Mean coordinate SNR=",round(meancorperf,2),
             "+-",round(sdcorperf,2),"dB"))

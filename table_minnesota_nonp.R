rm(list=ls())
library(gasper)
library(genlasso)
source("utils.R")

if (as.numeric(R.version$minor)>6){
  RNGkind(sample.kind = "Rounding")
}

x1 <- minnesota$xy[,1]
n <- length(x1)
A <- minnesota$A

t0 <- Sys.time()
L <- laplacian_mat(A)
val1 <- eigendec(L)
evalues <- val1$evalues
evectors <- val1$evectors
lmax <- max(evalues)
b <- 2
kmax <- floor(log(lmax)/log(b)) + 2
tf <- tight_frame(evalues, evectors, b=b)
diagWWt <- colSums(t(tf)^2)

# Reproduction of the table (left and right)
name <- "f2" #"f1"
load("minnesota/fminnesota.Rdata")
if (name == "f1") {
  f <- fminnesota$f1
}
if (name == "f2") {
  f <- fminnesota$f2
}

wcf <- analysis(f, tf)

AA <- as(A, "dgCMatrix")
grA <- graph_from_adjacency_matrix(adjmatrix=AA,
                                   mode = "undirected")

# Sigma values for the case f1 eta=0.01 and A squared
if (name == "f1") {
  sigmal <- c(0.005, 0.01, 0.02)
  #- increase maxiter to ensure GTF convergence
  maxiter <- 3000
}
# Sigma values for the case f2 eta=0.001 and A^4
if (name == "f2") {
  sigmal <- c(0.001, 0.002, 0.004)
  maxiter <- 2000
}

MC <- 10
set.seed(0)
yobs <- array(0, dim=c(n,length(sigmal),MC))
for (i in 1:length(sigmal)){
  for (it in 1:MC){
    noise <- rnorm(n, sd = sigmal[i])
    yobs[,i,it] <- f + noise
  }
}

SNRin <- array(0, dim=c(length(sigmal),MC))
fglob1 <- fglob2 <- fglsb1 <- fglsb2 <- array(0, dim=c(length(sigmal),MC))
flvob1 <- flvob2 <- flvsb1 <- flvsb2 <- array(0, dim=c(length(sigmal),MC))
for (i in 1:length(sigmal)){
  for (it in 1:MC){
    # noisy wc
    wcn <- analysis(yobs[,i,it],tf)
    # global coordinatewise
    thresh <- sort(abs(wcn))
    gl_thresh_b1 <- SURE_MSEthresh(wcn, 
                                   wcf, 
                                   thresh, 
                                   diagWWt, 
                                   b=1, 
                                   sigmal[i], 
                                   NA, 
                                   policy = "dependent")
    gl_thresh_b2 <- SURE_MSEthresh(wcn,
                                   wcf, 
                                   thresh, 
                                   diagWWt, 
                                   b=2, 
                                   sigmal[i], 
                                   NA,
                                   policy = "dependent")
    # level dependent coordinatewise
    lev_thresh_b1 <- lev_thresh_b2 <- list()
    wclevSURE_b1 <- wclevMSE_b1 <- rep(0, length(wcn))
    wclevSURE_b2 <- wclevMSE_b2 <- rep(0, length(wcn))
    for (k in 0:kmax){
      indscale <- seq(k*n+1, (k+1)*n)
      wc <- wcn[indscale]
      thresh_set <- sort(abs(wc))#[seq(1,length(wc),10)]
      lev_thresh_b1[[k+1]] <- SURE_MSEthresh(wc, 
                                         wcf[indscale], 
                                         thresh_set,
                                         diagWWt[indscale], 
                                         b=1, 
                                         sigmal[i], 
                                         NA,
                                         policy = "dependent")
      wclevMSE_b1[indscale] <- lev_thresh_b1[[k+1]]$wc[,lev_thresh_b1[[k+1]]$min[1]]
      wclevSURE_b1[indscale] <- lev_thresh_b1[[k+1]]$wc[,lev_thresh_b1[[k+1]]$min[2]]

      lev_thresh_b2[[k+1]] <- SURE_MSEthresh(wc, 
                                         wcf[indscale], 
                                         thresh_set,
                                         diagWWt[indscale], 
                                         b=2, 
                                         sigmal[i], 
                                         NA,
                                         policy = "dependent")
      wclevMSE_b2[indscale] <- lev_thresh_b2[[k+1]]$wc[,lev_thresh_b2[[k+1]]$min[1]]
      wclevSURE_b2[indscale] <- lev_thresh_b2[[k+1]]$wc[,lev_thresh_b2[[k+1]]$min[2]]
    }
    #global estimators according to each criteria
    hatf_gl_oracle_b1 <- synthesis(gl_thresh_b1$wc[,gl_thresh_b1$min[1]], tf)
    hatf_gl_sigma_b1  <- synthesis(gl_thresh_b1$wc[,gl_thresh_b1$min[2]], tf)
    hatf_gl_oracle_b2 <- synthesis(gl_thresh_b2$wc[,gl_thresh_b2$min[1]], tf)
    hatf_gl_sigma_b2  <- synthesis(gl_thresh_b2$wc[,gl_thresh_b2$min[2]], tf)
    #level dependent estimators according to each criteria
    hatf_lev_oracle_b1 <- synthesis(wclevMSE_b1, tf)
    hatf_lev_sigma_b1  <- synthesis(wclevSURE_b1, tf)
    hatf_lev_oracle_b2 <- synthesis(wclevMSE_b2, tf)
    hatf_lev_sigma_b2  <- synthesis(wclevSURE_b2, tf)
    fglob1[i,it]=SNR(f, hatf_gl_oracle_b1)
    fglob2[i,it]=SNR(f, hatf_gl_oracle_b2)
    fglsb1[i,it]=SNR(f, hatf_gl_sigma_b1)
    fglsb2[i,it]=SNR(f, hatf_gl_sigma_b2)
    flvob1[i,it]=SNR(f, hatf_lev_oracle_b1)
    flvob2[i,it]=SNR(f, hatf_lev_oracle_b2)
    flvsb1[i,it]=SNR(f, hatf_lev_sigma_b1)
    flvsb2[i,it]=SNR(f, hatf_lev_sigma_b2)
    SNRin[i,it] <- SNR(f, yobs[,i,it])
  }
}
meanSNRin <- apply(SNRin,1,mean)
sdSNRin <- apply(SNRin,1,sd)
meanfglob1 <- apply(fglob1,1,mean)
sdfglob1<- apply(fglob1,1,sd)
meanfglob2 <- apply(fglob2,1,mean)
sdfglob2 <- apply(fglob2,1,sd)
meanfglsb1 <- apply(fglsb1,1,mean)
sdfglsb1<- apply(fglsb1,1,sd)
meanfglsb2 <- apply(fglsb2,1,mean)
sdfglsb2<- apply(fglsb2,1,sd)
meanflvob1 <- apply(flvob1,1,mean)
sdflvob1<- apply(flvob1,1,sd)
meanflvob2 <- apply(flvob2,1,mean)
sdflvob2<- apply(flvob2,1,sd)
meanflvsb1 <- apply(flvsb1,1,mean)
sdflvsb1<- apply(flvsb1,1,sd)
meanflvsb2 <- apply(flvsb2,1,mean)
sdflvsb2<- apply(flvsb2,1,sd)
tf <- Sys.time()

t0W <- Sys.time()
fwiener <- array(0, dim=c(length(sigmal),MC))
fWoracle <- array(0, dim=c(length(sigmal),MC))
for (i in 1:length(sigmal)){
  for (it in 1:MC){
    # Wiener estimators
    ft_f <- evectors%*%f
    ft_y <- evectors%*%yobs[,i,it]
    thres_y <- ft_y*ft_f^2/(ft_f^2+sigmal[i]^2)
    rinf <- sigmal[i]^2*sum(ft_f^2/(ft_f^2+sigmal[i]^2))
    hatf_wiener <- synthesis(thres_y, evectors)

    fwiener[i,it]=SNR(f, hatf_wiener)
    fWoracle[i,it]= 20 * log10(normvec(as.vector(ft_f))/sqrt(rinf) )
  }
}
meanfwiener <- apply(fwiener,1,mean)
sdfwiener <- apply(fwiener,1,sd)
meanfWoracle <- apply(fWoracle,1,mean)
sdfWoracle <- apply(fWoracle,1,sd)
tfW <- Sys.time()

t0trend <- Sys.time()
fto <- array(0, dim=c(length(sigmal),MC))
ftrendS <- array(0, dim=c(length(sigmal),MC))
mint <- array(0, dim=c(length(sigmal),MC))
for (i in 1:length(sigmal)){
  for (it in 1:MC){
    # Trend filtering estimators
    trend <- fusedlasso(yobs[,i,it], graph=grA, maxsteps = maxiter)
    mintrend <- apply(trend$beta-f,2,function(x)(mean(x^2)))
    ftrend <- trend$beta[,which.min(mintrend)]
    erisk_trend <- apply(trend$beta-yobs[,i,it],
                         2,function(x)(sum(x^2)))
    dof_trend <- trend$df
    SUREtrend <- erisk_trend + 2*sigmal[i]^2*dof_trend
    ftrendSURE <- trend$beta[,which.min(SUREtrend)]

    fto[i,it] = SNR(f,ftrend)
    ftrendS[i,it]= SNR(f,ftrendSURE)
    mint[i,it]=which.min(mintrend)

  }
}
meanfto <- apply(fto,1,mean)
sdfto <- apply(fto,1,sd)
meanftrendSURE <- apply(ftrendS,1,mean)
sdftrendSURE <- apply(ftrendS,1,sd)
tftrend <- Sys.time()


res <- data.frame(meanSNRin,
                  meanfglob1,
                  meanfglob2,
                  meanfglsb1,
                  meanfglsb2,
                  meanflvob1,
                  meanflvob2,
                  meanflvsb1,
                  meanflvsb2,
                  meanfwiener,
                  meanfWoracle,
                  meanfto,
                  meanftrendSURE)


res_sd <- data.frame(sdSNRin,
                     sdfglob1,
                     sdfglob2,
                     sdfglsb1,
                     sdfglsb2,
                     sdflvob1,
                     sdflvob2,
                     sdflvsb1,
                     sdflvsb2,
                     sdfwiener,
                     sdfWoracle,
                     sdfto,
                     sdftrendSURE)

print(round(tf-t0,2))
print(round(tftrend-t0trend,2))


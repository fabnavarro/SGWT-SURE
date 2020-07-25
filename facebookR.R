rm(list=ls())
library(gasper)
library(R.matlab)

# From GTF matlab code file ExpFacebook2
load("facebook/facebook.rdata")
A <- full(A)
n <- nrow(A)

t0 <- Sys.time()
L <- laplacian_mat(A)
val1 <- eigendec(L)
evalues <- val1$evalues
evectors <- val1$evectors
lmax <- max(evalues)

b <- 5
kmax <- floor(log(lmax)/log(b)) + 2
tf <- tight_frame(evalues, evectors, b=b)

diagWWt <- colSums(t(tf)^2)
kmax <- nrow(tf)/n-1

#- Comment/Uncomment below

############################################
# Sparse poisson exemple: one relalization #
############################################
# f for the sparse poisson exemple
#f <- read.table("facebook/fFacebook.csv")
#f <- f$V1
#- ten realisation of noisy f for sparse poisson exemple
#ya <- read.table("facebook/yFaceebok.csv",sep=",")
#- sigma value for the sparse poisson exemple
#sigma_list <- c(0.01,0.1,0.3,0.5,0.7,0.9)
#########################################
# Sparse poisson exemple: 5 relaization #
#########################################
#f <- read.table("facebook/fFacebookSparsePoisson.csv")
#f <- f$V1
#ya <- readMat("facebook/yFacebookSparsePoisson.mat")
#########################################
# Dense poisson exemple: 5 relaization #
#########################################
f <- read.table("facebook/fFacebookDensePoisson.csv")
f <- f$V1
ya <- readMat("facebook/yFacebookDensePoisson.mat")
# ###############################################
# # homogeneous random walks: five realizations #
# ###############################################
#f <- read.table("facebook/fFaceebokRandWalkhomo.csv")
#f <- f$V1
#ya <- readMat("facebook/yRandWalkHomo.mat")
###############################################
# inhomogeneous random walks: five realizations #
###############################################
#f <- read.table("facebook/fFacebookRandWalkInomo.csv")
#f <- f$V1
#ya <- readMat("facebook/yRandWalkInomo.mat")

SnR_list=c(30,25,20,15,10,5,0,-5,-10) # in dB
# 10 log (signal_sigma^2/noise_sigma^2)
sigma_list=sqrt(10.^(-SnR_list/10)) #/sqrt(n);

bet <- 2
mse_all <- rep(0,length(sigma_list))
for (is in 1:length(sigma_list)){
  sigma <- sigma_list[is]
  mse_b2 <- rep(0,5)
  for (it in 1:5){
    y <- f + ya$Z[,it,is]
    
    wcn <- analysis(y,tf)
    wcf <- analysis(f,tf)
    
    hatsigma <- sqrt(GVN(y, A, L))
    
    # Level dependent coordinatewise
    lev_thresh_b2 <- list()
    opt_th_MSE <- opt_th_SURE <- opt_th_hatSURE <- rep(0, kmax+1)
    wclevSURE_b2 <- wclevMSE_b2 <- wclevhatSURE_b2 <- rep(0, length(wcn))
    tresh_set_all <- rep(0, length(wcn))
    for (k in 0:kmax){
      indscale <- seq(k*n+1, (k+1)*n)
      wc <- wcn[indscale]
      thresh_set <- sort(abs(wc)/sqrt(diagWWt[indscale]))
      #thresh_set <- sort(abs(wc))
      tresh_set_all[indscale] <- thresh_set
      
      lev_thresh_b2[[k+1]] <- SURE_MSEthresh(wc, 
                                          wcf[indscale], 
                                          thresh_set,
                                          diagWWt[indscale], 
                                          b=bet, 
                                          sigma, 
                                          hatsigma,
                                          policy = "dependent")
      wclevMSE_b2[indscale] <- lev_thresh_b2[[k+1]]$wc[,lev_thresh_b2[[k+1]]$min[1]]
      opt_th_MSE[k+1] <- lev_thresh_b2[[k+1]]$thr[1]
    }
    # Level dependent estimators
    hatf_lev_oracle_b2 <- synthesis(wclevMSE_b2, tf)
    mse_b2[it] <- mean((hatf_lev_oracle_b2-f)^2)
  }
  mse_all[is] = mean(mse_b2)
}
tfin <- Sys.time()

#- Uncomment to save output
#writeMat(mse = mse_all, "densePoissonR.mat")
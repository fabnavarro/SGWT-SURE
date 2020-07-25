rm(list=ls())
library(genlasso)
library(gasper)
source("utils.R")

set.seed(0)
A <- full(pittsburgh$sA)
n <- nrow(A)
f <- as.vector(pittsburgh$f)
sigma <- 0.2

t0 <- Sys.time()
L <- laplacian_mat(A)
val1 <- eigendec(L)
evalues <- val1$evalues
evectors <- val1$evectors
lmax <- max(evalues)
b <- 2
kmax <- floor(log(lmax)/log(b)) + 2
tf <- tight_frame(evalues, evectors)
diagWWt <- colSums(t(tf)^2)
MSElev_oracle_b2 <- array(0, dim = c(10,n))
for (i in 1:10){
  y <- pittsburgh$y[,i]
  #sigma <- sd(f-y[,i]) # sigma is not provided
  #y <- f + rnorm(n,0,sd = sigma)
  wcn <- analysis(y,tf)
  wcf <- analysis(f,tf)
  # Level dependent coordinatewise
  lev_thresh_b2 <- list()
  opt_th_MSE <- rep(0, kmax+1)
  wclevMSE_b2 <- rep(0, length(wcn))
  tresh_set_all <- rep(0, length(wcn))
  for (k in 0:kmax){
    indscale <- seq(k*n+1, (k+1)*n)
    wc <- wcn[indscale]
    thresh_set <- sort(abs(wc)/sqrt(diagWWt[indscale]))
    tresh_set_all[indscale] <- thresh_set

    lev_thresh_b2[[k+1]] <- SURE_MSEthresh(wc,
                                       wcf[indscale],
                                       thresh_set,
                                       diagWWt[indscale],
                                       b=2, 
                                       sigma,
                                       NA,
                                       policy = "dependent")
    wclevMSE_b2[indscale] <- lev_thresh_b2[[k+1]]$wc[,lev_thresh_b2[[k+1]]$min[1]]
    opt_th_MSE[k+1] <- lev_thresh_b2[[k+1]]$thr[1]
  }
  hatf_lev_oracle_b2 <- synthesis(wclevMSE_b2, tf)
  MSElev_oracle_b2[i,] <- hatf_lev_oracle_b2
}
tfin <- Sys.time()

# Trend filtering
t0trend <- Sys.time()
trenditermax <- 2000
MSEtr <- array(0, dim = c(10,n))
mintrenditer <- rep(0,10)
AA <- as(A, "dgCMatrix")
grA = graph_from_adjacency_matrix(adjmatrix=AA,
                                  mode = "undirected")
for (i in 1:10){
  y <- pittsburgh$y[,i]
  trend <- fusedlasso(y, graph=grA, maxsteps = trenditermax)
  mintrend <- apply(trend$beta-f,2,function(x)(mean(x^2)))
  xmintrend <- which.min(mintrend)
  if (xmintrend==trenditermax){
    warning("minimum at the edge,
            increase the number of iterations")
  }
  ftrend <- trend$beta[,xmintrend]
  MSEtr[i,] <- ftrend
}
tftrend <- Sys.time()

SNRb2 <- apply(MSElev_oracle_b2,1,function(x)SNR(f,x))
SNRtrend <- apply(MSEtr,1,function(x)SNR(f,x))
print(paste0("Oracle Trend filtering SNR=",
             round(mean(SNRtrend),2),"dB"))
print(round(tftrend-t0trend,2))
print(paste0("Oracle SGWT beta=2 SNR=",
             round(mean(SNRb2),2),"dB"))
print(round(tfin-t0,2))

rm(list=ls())
library(gasper)
library(foreach)
library(doMC)
library(genlasso)
source("utils.R")

if (as.numeric(R.version$minor)>6){
  RNGkind(sample.kind = "Rounding")
}

f <- as.vector(NYCdata$f)
A <- NYCdata$A

n <- nrow(A)
L <- laplacian_mat(A)
val1 <- eigendec(L)
evalues <- val1$evalues
evectors <- val1$evectors
lmax <- max(evalues)
b <- 30
kmax <- floor(log(lmax)/log(b)) + 2
tf <- tight_frame(evalues, evectors, b=b)
diagWWt <- colSums(t(tf)^2)

modulef <- smoothmodulus(f, A, L)
wcf <- analysis(f, tf)
sigmal <- 2.5

# Sparse Adj for trend filtering
AA <- as(A, "dgCMatrix")
grA <- graph_from_adjacency_matrix(adjmatrix=AA,
                                   mode = "undirected")
maxiter <- 2000

registerDoMC(32)
MC <- 25
# Each process seed must be fixed to ensure reproducibility of results
# (each process have a new random seed by default)
set.seed(0)
rseed <- sample(1:10000, MC)
res <- foreach(i=1:length(sigmal), .combine=cbind) %dopar% {
  foreach(it=1:MC, .combine=c) %dopar%  {
    set.seed(rseed[it])
    noise <- rnorm(n, sd = sigmal[i])
    y <- f + noise
    hatsigma <- sqrt(GVN(y, A, L))
    # noisy wc
    wcn <- analysis(y,tf)
    # level dependent coordinatewise
    lev_thresh_b1 <- lev_thresh_b2 <- list()
    wclevSURE_b2 <- wclevMSE_b2 <- rep(0, length(wcn))
    for (k in 0:kmax){
      indscale <- seq(k*n+1, (k+1)*n)
      wc <- wcn[indscale]
      thresh_set <- sort(abs(wc)/sqrt(diagWWt[indscale]))
      lev_thresh_b2[[k+1]] <- SURE_MSEthresh(wc, 
                                          wcf[indscale], 
                                          thresh_set,
                                          diagWWt[indscale], 
                                          b=2, 
                                          sigmal[i], 
                                          hatsigma,
                                          policy = "dependent")
      wclevMSE_b2[indscale] <- lev_thresh_b2[[k+1]]$wc[,lev_thresh_b2[[k+1]]$min[1]]
      wclevSURE_b2[indscale] <- lev_thresh_b2[[k+1]]$wc[,lev_thresh_b2[[k+1]]$min[2]]
    }
    #level dependent estimators according to each criteria
    hatf_lev_oracle_b2 <- synthesis(wclevMSE_b2, tf)
    hatf_lev_sigma_b2  <- synthesis(wclevSURE_b2, tf)

    # Wiener estimator
    ft_f <- evectors%*%f
    ft_y <- evectors%*%y
    thres_y <- ft_y*ft_f^2/(ft_f^2+sigmal[i]^2)
    hatf_wiener <- synthesis(thres_y, evectors)

    # Trend Filtering
    trend <- fusedlasso(y, graph=grA, maxsteps = maxiter)
    mintrend <- apply(trend$beta-f,2,function(x)(mean(x^2)))
    ftrend <- trend$beta[,which.min(mintrend)]

    data.frame("SNRin"=SNR(f, y),
               "flvob2"=SNR(f, hatf_lev_sigma_b2),
               "fwiener"=SNR(f, hatf_wiener),
               "fto"= SNR(f,ftrend))
  }
}

meanSNRin <- round(mean(unlist(res[seq(1,length(res),4)])),2)
sdSNRin <- round(sd(unlist(res[seq(1,length(res),4)])),2)
meanflvob2 <- round(mean(unlist(res[seq(2,length(res),4)])),2)
sdflvob2 <- round(sd(unlist(res[seq(2,length(res),4)])),2)
meanWiener <- round(mean(unlist(res[seq(3,length(res),4)])),2)
sdWiener <- round(sd(unlist(res[seq(3,length(res),4)])),2)
meantrend <- round(mean(unlist(res[seq(4,length(res),4)])),2)
sdtrend <- round(sd(unlist(res[seq(4,length(res),4)])),2)

print(paste0("Input SNR_in=",
             meanSNRin,"+-",sdSNRin,"dB"))
print(paste0("SURE SNR=",
             meanflvob2,"+-",sdflvob2,"dB"))
print(paste0("Oracle Wiener SNR=",
             meanWiener,"+-",sdWiener,"dB"))
print(paste0("Oracle Trend SNR=",
             meantrend,"+-",sdtrend,"dB"))

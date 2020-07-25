rm(list=ls())
library(gasper)
library(foreach)
library(doMC)
library(genlasso)
source("utils.R")

if (as.numeric(R.version$minor)>6){
  RNGkind(sample.kind = "Rounding")
}

A <- full(pittsburgh$sA)
n <- nrow(A)
L <- laplacian_mat(A)
val1 <- eigendec(L)
evalues <- val1$evalues
evectors <- val1$evectors
lmax <- max(evalues)
b <- 2
kmax <- floor(log(lmax)/log(b)) + 2
tf <- tight_frame(evalues, evectors, b=b)
diagWWt <- colSums(t(tf)^2)

load("pittsburgh/pittsburghGEO.Rdata")
f <- pittsburghGEO$f1
wcf <- analysis(f, tf)

AA <- as(A, "dgCMatrix")
grA <- graph_from_adjacency_matrix(adjmatrix=AA,
                                   mode = "undirected")
sigmal <- c(0.004, 0.005, 0.01)
maxiter <- 2000

registerDoMC(32)
MC <- 10
# Each process seed must be fixed to ensure reproducibility of results
# (each process have a new random seed by default)
set.seed(0)
rseed <- sample(1:10000, MC)
res <- foreach(i=1:length(sigmal), .combine=cbind) %dopar% {
  foreach(it=1:MC, .combine=c) %dopar%  {
    set.seed(rseed[it])
    noise <- rnorm(n, sd = sigmal[i])
    y <- f + noise
    # noisy wc
    wcn <- analysis(y,tf)
    # sigma estimation
    hatsigma <- sqrt(GVN(y, A, L))
    hatsigma_sca <- sig2estim_scale(wcf, wcn, kmax, evalues, n)
    # Level dependent coordinatewise
    lev_thresh_b2 <- list()
    lev_thresh2_b2 <- list()
    opt_th_MSE <- rep(0, kmax+1)
    wclevMSE_b2 <- wclevSURE_b2 <- rep(0, length(wcn))
    wclevhatSURE_b2 <- wclevhatSURE2_b2 <- rep(0, length(wcn))
    tresh_set_all <- rep(0, length(wcn))
    for (k in 0:kmax){
      indscale <- seq(k*n+1, (k+1)*n)
      wc <- wcn[indscale]
      #thresh_set <- sort(abs(wc))
      thresh_set <- sort(abs(wc)/sqrt(diagWWt[indscale]))
      tresh_set_all[indscale] <- thresh_set

      lev_thresh_b2[[k+1]] <- SURE_MSEthresh(wc,
                                         wcf[indscale],
                                         thresh_set,
                                         diagWWt[indscale],
                                         b=2,
                                         sigmal[i],
                                         hatsigma,
                                         policy = "dependent")
      lev_thresh2_b2[[k+1]] <- SURE_MSEthresh(wc,
                                          wcf[indscale],
                                          thresh_set,
                                          diagWWt[indscale],
                                          b=2,
                                          sigmal[i],
                                          hatsigma_sca[kmax+1],
                                          policy = "dependent")

      wclevMSE_b2[indscale] <- lev_thresh_b2[[k+1]]$wc[,lev_thresh_b2[[k+1]]$min[1]]
      wclevSURE_b2[indscale] <- lev_thresh_b2[[k+1]]$wc[,lev_thresh_b2[[k+1]]$min[2]]
      wclevhatSURE_b2[indscale] <- lev_thresh_b2[[k+1]]$wc[,lev_thresh_b2[[k+1]]$min[3]]
      wclevhatSURE2_b2[indscale] <- lev_thresh2_b2[[k+1]]$wc[,lev_thresh2_b2[[k+1]]$min[3]]
      opt_th_MSE[k+1] <- lev_thresh_b2[[k+1]]$thr[1]
    }
    # Level dependent estimators according to each criteria
    hatf_lev_oracle_b2 <- synthesis(wclevMSE_b2, tf)
    hatfSURE <- synthesis(wclevSURE_b2, tf)
    hatfhSURE <- synthesis(wclevhatSURE_b2, tf)
    hatfh2SURE <- synthesis(wclevhatSURE2_b2, tf)

    # Trend filtering
    trend <- fusedlasso(y, graph=grA, maxsteps = maxiter)

    mintrend <- apply(trend$beta-f,2,function(x)(mean(x^2)))
    ftrend <- trend$beta[,which.min(mintrend)]

    erisk_trend <- apply(trend$beta-y,
                         2,function(x)(sum(x^2)))
    dof_trend <- trend$df
    SUREtrend <- erisk_trend + 2*sigmal[i]^2*dof_trend
    hSUREtrend <- erisk_trend + 2*hatsigma^2*dof_trend
    h2SUREtrend <-erisk_trend + 2*hatsigma_sca[kmax+1]^2*dof_trend
    ftrendSURE <- trend$beta[,which.min(SUREtrend)]
    ftrendhSURE <- trend$beta[,which.min(hSUREtrend)]
    ftrendh2SURE <- trend$beta[,which.min(h2SUREtrend)]

    data.frame("SNRin" = SNR(f, y),
               "MSEbeta2" = SNR(f, hatf_lev_oracle_b2),
               "MSEtrend" = SNR(f, ftrend),
               "SUREbeta2" = SNR(f, hatfSURE),
               "SUREtrend" = SNR(f, ftrendSURE),
               "SUREbeta2sig1" = SNR(f, hatfhSURE),
               "SUREtrendsig1" = SNR(f, ftrendhSURE),
               "SUREbeta2sig2" = SNR(f, hatfh2SURE),
               "SUREtrendsig2" = SNR(f, ftrendh2SURE),
               "sig1" = hatsigma,
               "sig2" = hatsigma_sca[kmax+1],
               "mint" = which.min(mintrend))
  }
}

# Uncomment to save the results
#save(res,file = "resPitt.Rdata", version = 2)

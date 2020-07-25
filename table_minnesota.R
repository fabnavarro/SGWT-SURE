rm(list=ls())
library(gasper)
library(foreach)
library(doMC)
library(genlasso)
source("utils.R")

if (as.numeric(R.version$minor)>6){
  RNGkind(sample.kind = "Rounding")
}

x1 <- minnesota$xy[,1]
n <- length(x1)
A <- minnesota$A
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
  #
  maxiter <- 3000
}
# Sigma values for the case f2 eta=0.001 and A^4
if (name == "f2") {
  sigmal <- c(0.001, 0.002, 0.004)
  maxiter <- 2000
}

registerDoMC(cores=2) # 2 by default
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
    # global coordinatewise
    #thresh <- sort(abs(wcn))
    thresh <- sort(abs(wcn)/sqrt(diagWWt))
    gl_thresh_b1 <-  SURE_MSEthresh(wcn, 
                                    wcf, 
                                    thresh, 
                                    diagWWt, 
                                    b=1, 
                                    sigmal[i], 
                                    NA, 
                                    policy = "dependent")
    gl_thresh_b2 <-  SURE_MSEthresh(wcn, 
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
      #thresh_set <- sort(abs(wc))
      thresh_set <- sort(abs(wc)/sqrt(diagWWt[indscale]))
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
    # Wiener estimators
    ft_f <- evectors%*%f
    ft_y <- evectors%*%y
    thres_y <- ft_y*ft_f^2/(ft_f^2+sigmal[i]^2)
    rinf <- sigmal[i]^2*sum(ft_f^2/(ft_f^2+sigmal[i]^2))
    hatf_wiener <- synthesis(thres_y, evectors)
    # Trend filtering estimators
    #- Note: The results for the trend filtering are obtained using
    #        gtf matlab code (that give slightly better results for k=0
    #        and handle the cases k=1 and k=2)
    # trend <- fusedlasso(y, graph=grA, maxsteps = maxiter)
    # mintrend <- apply(trend$beta-f,2,function(x)(mean(x^2)))
    # ftrend <- trend$beta[,which.min(mintrend)]
    # erisk_trend <- apply(trend$beta-y,
    #                      2,function(x)(sum(x^2)))
    # dof_trend <- trend$df
    # SUREtrend <- erisk_trend + 2*sigmal[i]^2*dof_trend
    # ftrendSURE <- trend$beta[,which.min(SUREtrend)]

    data.frame("SNRin"=SNR(f, y),
               "fglob1"=SNR(f, hatf_gl_oracle_b1),
               "fglob2"=SNR(f, hatf_gl_oracle_b2),
               "fglsb1"=SNR(f, hatf_gl_sigma_b1),
               "fglsb2"=SNR(f, hatf_gl_sigma_b2),
               "flvob1"=SNR(f, hatf_lev_oracle_b1),
               "flvob2"=SNR(f, hatf_lev_oracle_b2),
               "flvsb1"=SNR(f, hatf_lev_sigma_b1),
               "flvsb2"=SNR(f, hatf_lev_sigma_b2),
               "fwiener"=SNR(f, hatf_wiener),
               "fWoracle"= 20 * log10(normvec(as.vector(ft_f))/sqrt(rinf) ),
               "fto"= SNR(f,f),
               "ftrendSURE"= SNR(f,f)
               #"mint"=which.min(mintrend)
               )
  }
}

# Uncomment to save the results
#save(res,file = "res0001A4MC10.Rdata", version = 2)
#save(res,file = "res001A2MC10.Rdata", version = 2)

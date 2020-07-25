rm(list=ls())
library(gasper)
library(ggplot2)
library(gridExtra)
source("utils.R")

if (as.numeric(R.version$minor)>6){
  RNGkind(sample.kind = "Rounding")
}

# Set the random seed for reproductibility
set.seed(0)
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

load("minnesota/fminnesota.Rdata")
f <- fminnesota$f1

# Noise structure
id <- diag(rep(1,n))
# Uncomment for the second experiment
#alpha <- 0.5
alpha <- 0.1
gamma <- (id+2*alpha*A+alpha^2*matmult(A,A))
tfc <- matmult(tf,gamma)
sigma <- 0.02
noise <- (id+alpha*A)%*%rnorm(n, sd = sigma)
y <- f + as.vector(noise)

# SURE weights
diagWWtcor <- colSums(t(tfc)*t(tf))
diagWWt <- colSums(t(tf)^2)

wcn <- analysis(y,tf)
wcf <- analysis(f,tf)

# Level dependent coordinatewise
lev_thresh_b2 <- lev_thresh_b2cor <- list()
opt_th_MSEcor <- opt_th_SUREcor <- opt_th_SURE <- rep(0, kmax+1)
wclevSURE_b2cor <- wclevSURE_b2 <- wclevMSE_b2cor <- rep(0, length(wcn))
tresh_set_all <- rep(0, length(wcn))
for (k in 0:kmax){
  indscale <- seq(k*n+1, (k+1)*n)
  wc <- wcn[indscale]
  thresh_set <- sort(abs(wc)/sqrt(diagWWt[indscale]))
  tresh_set_all[indscale] <- thresh_set
  lev_thresh_b2cor[[k+1]] <- SURE_MSEthresh(wc,
                                        wcf[indscale],
                                        thresh_set,
                                        diagWWtcor[indscale],
                                        b=2, 
                                        sigma, 
                                        NA,
                                        policy = "dependent")
  lev_thresh_b2[[k+1]] <- SURE_MSEthresh(wc,
                                     wcf[indscale],
                                     thresh_set,
                                     diagWWt[indscale],
                                     b=2, 
                                     sigma, 
                                     NA,
                                     policy = "dependent")

  wclevMSE_b2cor[indscale] <- lev_thresh_b2cor[[k+1]]$wc[,lev_thresh_b2cor[[k+1]]$min[1]]
  wclevSURE_b2cor[indscale] <- lev_thresh_b2cor[[k+1]]$wc[,lev_thresh_b2cor[[k+1]]$min[2]]
  wclevSURE_b2[indscale] <- lev_thresh_b2[[k+1]]$wc[,lev_thresh_b2[[k+1]]$min[2]]

  opt_th_MSEcor[k+1] <- lev_thresh_b2cor[[k+1]]$thr[1]
  opt_th_SUREcor[k+1] <- lev_thresh_b2cor[[k+1]]$thr[2]
  opt_th_SURE[k+1] <- lev_thresh_b2[[k+1]]$thr[2]
}
# Level dependent estimators
hatf_lev_oracle_b2cor <- synthesis(wclevMSE_b2cor, tf)
hatf_lev_sigma_b2cor  <- synthesis(wclevSURE_b2cor, tf)
hatf_lev_sigma_b2     <- synthesis(wclevSURE_b2, tf)
print(paste0("Input SNR_in=",
             round(SNR(f, y),2),"dB"))
print(paste0("Oracle: SNR=",
             round(SNR(f, hatf_lev_oracle_b2cor),2),"dB"))
print(paste0("SUREcor: SNR=",
             round(SNR(f, hatf_lev_sigma_b2cor),2),"dB"))
print(paste0("SURE: SNR=",
             round(SNR(f, hatf_lev_sigma_b2),2),"dB"))


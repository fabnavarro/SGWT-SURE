rm(list=ls())
library(gasper)
library(ggplot2)
library(genlasso)
library(RColorBrewer)
library(gridExtra)
source("utils.R")

if (as.numeric(R.version$minor)>6){
  RNGkind(sample.kind = "Rounding")
}

set.seed(0)
load("pittsburgh/pittsburghGEO.Rdata")
A <- pittsburghGEO$A

n <- nrow(A)
sigma <- 0.01
f <- pittsburghGEO$f1
y <- f + rnorm(n, 0, sd = sigma)

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
wcn <- analysis(y, tf)
wcf <- analysis(f, tf)
hatsigma <- sqrt(GVN(y, A, L))

# Level dependent coordinatewise
lev_thresh_b2 <- list()
opt_th_MSE <- rep(0, kmax+1)
wclevMSE_b2 <- wclevSURE_b2 <- wclevhatSURE_b2 <- rep(0, length(wcn))
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
                                     sigma,
                                     hatsigma,
                                     policy = "dependent")
  wclevMSE_b2[indscale] <- lev_thresh_b2[[k+1]]$wc[,lev_thresh_b2[[k+1]]$min[1]]
  wclevSURE_b2[indscale] <- lev_thresh_b2[[k+1]]$wc[,lev_thresh_b2[[k+1]]$min[2]]
  wclevhatSURE_b2[indscale] <- lev_thresh_b2[[k+1]]$wc[,lev_thresh_b2[[k+1]]$min[3]]
  opt_th_MSE[k+1] <- lev_thresh_b2[[k+1]]$thr[1]
}
# Level dependent estimators according to each criteria
hatf_lev_oracle_b2 <- synthesis(wclevMSE_b2, tf)
hatfSURE <- synthesis(wclevSURE_b2, tf)
hatfhSURE <- synthesis(wclevhatSURE_b2, tf)
tf <- Sys.time()

# Fused lasso
t0trend <- Sys.time()
trenditermax <- 2000
AA <- as(A, "dgCMatrix")
grA = graph_from_adjacency_matrix(adjmatrix=AA,
                                  mode = "undirected")
trend <- fusedlasso(y,
                    graph=grA,
                    maxsteps = trenditermax)

mintrend <- apply(trend$beta-f,2,function(x)(mean(x^2)))
xmintrend <- which.min(mintrend)
if (xmintrend==trenditermax){
  warning("minimum at the edge, increase the number of iterations")
}
ftrend <- trend$beta[,xmintrend]
erisk_trend <- apply(trend$beta-y,
                     2,function(x)(sum(x^2)))
dof_trend <- trend$df
SUREtrend <- erisk_trend + 2*sigma^2*dof_trend
hSUREtrend <- erisk_trend + 2*hatsigma^2*dof_trend
ftrendSURE <- trend$beta[,which.min(SUREtrend)]
ftrendhSURE <- trend$beta[,which.min(hSUREtrend)]
tftrend <- Sys.time()
print(paste0("Input SNR_in=",
             round(SNR(f, y),2),"dB"))
print(paste0("Oracle Trend filtering SNR=",
             round(SNR(f,ftrend),2),"dB",
             " (",round(tftrend-t0trend,2),"s)"))
print(paste0("Oracle SGWT beta=2 SNR=",
             round(SNR(f,hatf_lev_oracle_b2),2),"dB",
             " (",round(tf-t0,2),"s)"))


# Plot
allegheny_tracts <- pittsburghGEO$geo
allegheny_tracts$f <- f
allegheny_tracts$y <- y
g <- hatfSURE
allegheny_tracts$g <- g
h <- ftrendSURE
allegheny_tracts$h <- h

lim_inf <- min(f,y,g,h)
lim_sup <- max(f,y,g,h)

p1 <- allegheny_tracts %>%
  ggplot() +
  geom_sf(aes(fill = f))+
  scale_fill_gradientn(colours = brewer.pal(n=11,name = "RdBu"),
                       limits = c(lim_inf,lim_sup))+
  theme_void()+
  theme(panel.grid.major = element_line(colour = "white"))

p2 <- allegheny_tracts %>%
  ggplot() +
  geom_sf(aes(fill = y))+
  scale_fill_gradientn(colours = brewer.pal(n=11,name = "RdBu"),
                       limits = c(lim_inf,lim_sup))+
  theme_void()+
  theme(panel.grid.major = element_line(colour = "white"))

p3 <- allegheny_tracts %>%
  ggplot() +
  geom_sf(aes(fill = g))+#, color = NA
  scale_fill_gradientn(colours = brewer.pal(n=11,name = "RdBu"),
                       limits = c(lim_inf,lim_sup))+
  theme_void()+
  theme(legend.title=element_blank())+
  theme(panel.grid.major = element_line(colour = "white"))

p4 <- allegheny_tracts %>%
  ggplot() +
  geom_sf(aes(fill = h))+
  scale_fill_gradientn(colours = brewer.pal(n=11,name = "RdBu"),
                       limits = c(lim_inf,lim_sup))+
  theme_void()+
  theme(legend.title=element_blank())+
  theme(panel.grid.major = element_line(colour = "white"))

grid.arrange(p1,p2,p3,p4)

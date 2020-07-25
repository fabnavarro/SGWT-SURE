#- utils
zetav <- function(x, k, b) {
  if (k == 0) {
    return(gv(x, b))
  } else {
    return(gv(b^(-k) * x, b) - gv(b^(-k + 1) * x, b))
  }
}

gv <- function(x, b) {
  low <- outer(x, 0, "<")
  mid1 <- outer(x, 0, ">=")
  mid2 <- outer(x, 1/b, "<=")
  mid3 <- outer(x, 1/b, ">=")
  mid4 <- outer(x, 1, "<=")
  up <- outer(x, 1, ">")

  gg <- rep(0, length(x))
  gg[low] <- 1
  gg[mid1 & mid2] <- 1
  gg[mid3 & mid4] <- b * x[mid3 & mid4]/(1 - b) + b/(b - 1)
  gg[up] <- 0
  return(gg)
}

sig2estim_scale <- function(wcf, wcn, kmax, evalues, n){
  scaf <- rep(0,kmax+1)
  scay <- rep(0,kmax+1)
  for(j in 1:(kmax+1)) {
    scaf[j]<-sqrt(sum(wcf[(n*(j-1)+1):(n*j)]^2)/sum(zetav(evalues,j-1,b)))
    scay[j]<-sqrt(sum(wcn[(n*(j-1)+1):(n*j)]^2)/sum(zetav(evalues,j-1,b)))
  }
  return(scay)
}

normvec <- function(x) {
  return(sqrt(sum(x^2)))
}

SNR <- function(x, y) {
  v <- 20 * log10(normvec(as.vector(x))/normvec(as.vector(x) - as.vector(y)))
  return(v)
}

blockthresh <- function(y, t, beta) {
  normBlock<-sqrt(sum(y^2))
  x <- max(0, 1 - t^beta/normBlock^beta) * y
  return(x)
}

mse2snr<-function(sumsqx,sumsqy){
  return(10*log10(sumsqx/sumsqy))
}

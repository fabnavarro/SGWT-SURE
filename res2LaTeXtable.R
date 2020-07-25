rm(list=ls())
library(xtable)
library(dplyr)
library(tidyverse)
#- Specify the path to the result data
path <- "minnesota/res001A2MC10.Rdata"  # table 1 left
#path <- "minnesota/res0001A4MC10.Rdata" # table 1 right
#path <- "pittsburgh/resPitt.Rdata" # table 2
load(path)

ind <- unique(rownames(res))
indc <- length(colnames(res))
MC <- 10
df <- dfsd <- matrix(0,
                     ncol=indc,
                     nrow=length(ind))
for (i in 1:length(ind)){
  df[i,] <- colMeans(matrix(unlist(res[which(rownames(res)==ind[i]),]),
                            nrow = MC))
  dfsd[i,] <- apply(matrix(unlist(res[which(rownames(res)==ind[i]),]),
                           nrow = MC),2,sd)
}
if (path=="pittsburgh/resPitt.Rdata"){
  id <- 3
  sig1 <- df[nrow(df)-id+1,]
  sig2 <- df[nrow(df)-id+2,]
  print(paste0("Variance estimation: sig1=",round(sig1,4),
               " sig2=", round(sig2,4)))
}else{
  id <- 2
}
dfout <- round(df[1:(nrow(df)-id),],2)
dfsdout <- round(dfsd[1:(nrow(df)-id),],2)
dfres <- rbind(dfout,dfsdout)
dim(dfres) <- c(nrow(dfout), 2*ncol(dfout))

dfres <- data.frame(dfres)
dfres_sd <- dfres %>%
  mutate_all(as.character) %>%
  unite(A1,c("X1","X2"),sep = "+-") %>%
  unite(A2,c("X3","X4"),sep = "+-") %>%
  unite(A3,c("X5","X6"),sep = "+-")

# Labels Table 1
# 1: SNR_in
# 2: global MSE_b1
# 3: global MSE_b2
# 4: global SURE_b1
# 5: global SURE_b2
# 6: lev MSE_b1
# 7: lev MSE_b2
# 8: lev SURE_b1
# 9: lev SURE_b2
# 10: Wiener
# 11: Wiener oracle

# Labels Table 2
# 1: SNR_in
# 2: lev oracle b2
# 3: oracle trend (genlasso)
# 4: lev SURE_b2
# 5: trend SURE (genlasso)
# 6: lev SURE_b2 sig1
# 7: trend SURE sig1 (genlasso)
# 8: lev SURE_b2 sig2
# 9: trend SURE sig2 (genlasso)

print(xtable(dfres_sd))

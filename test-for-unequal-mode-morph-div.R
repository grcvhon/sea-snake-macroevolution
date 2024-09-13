# R script for test of unequal mode of morphological diversification
# Authors: Vhon Oliver S. Garcia and Simone Blomberg
# Date: 4 February 2022

# Load packages
library(ape)
library(phytools)
library(cluster)
library(MASS)

#### functions ####

# morphometric branch lengths
L_D1 <- function (tree, X) {
  ## Code modified from Liam Revell's phylomorphospace function
  blen <- vector("numeric", length=length(tree$edge.length))
  A <- apply(X, 2, fastAnc, tree = tree)
  aa <- setNames(c(X[tree$tip.label, 1], A[, 1]), c(1:length(tree$tip.label),
                                                    rownames(A)))
  bb <- setNames(c(X[tree$tip.label, 2], A[, 2]), c(1:length(tree$tip.label),
                                                    rownames(A)))
  XX <- matrix(aa[as.character(tree$edge)], nrow(tree$edge), 2)
  YY <- matrix(bb[as.character(tree$edge)], nrow(tree$edge), 2)
  
  for (i in 1:nrow(XX)) {
    blen[i] <- sqrt((XX[i, 1]-XX[i, 2])^2 + (YY[i, 1]-YY[i, 2])^2)
  }
  sum(blen)
}

# volume (here: area since there are only two dimensions)
V_D1 <- function (matrix) {
  area_perm <- ellipsoidhull(matrix) # generates an ellipsoid object
  volume(area_perm) # calculates the volume of the generated ellipsoid object
}

# volume (here: area since there are only two dimensions)
# V_D1 is volume of the hyperellipsoid

#### data ####
tree <- read.tree("data/sea snake MCC tree")
tree_ultra <- chronos(tree)
is.ultrametric(tree_ultra) # TRUE

tree_ultra_AE <- extract.clade(tree_ultra,86)
tree_ultra_hyd <- extract.clade(tree_ultra,51)

morphdata_AE <- read.csv("data/LOG10morphdata_AE.csv")
rownames(morphdata_AE) <- morphdata_AE$X
morphdata_AE <- morphdata_AE[-1]
morphdata_AE
raw_morphdata_AE <- morphdata_AE[,-c(3,4)]
log10_morphdata_AE <- morphdata_AE[,-c(1,2)]

morphdata_hyd <- read.csv("data/LOG10morphdata_hyd.csv")
rownames(morphdata_hyd) <- morphdata_hyd$X
morphdata_hyd <- morphdata_hyd[-1]
morphdata_hyd
raw_morphdata_hyd <- morphdata_hyd[,-c(3,4)]
log10_morphdata_hyd <- morphdata_hyd[,-c(1,2)]

#### analysis ####

# Lineage density 1

# AE (raw)
sumblen_AE <- L_D1(tree_ultra_AE,raw_morphdata_AE) # 1728.821
ellipse_AE <- ellipsoidhull(as.matrix(raw_morphdata_AE))
area_AE <- volume(ellipse_AE) # 176.4902
ldensity1_AE <- sumblen_AE/area_AE # 9.795567

# Hydrophis (raw)
sumblen_hyd <- L_D1(tree_ultra_hyd,raw_morphdata_hyd) # 11529.33
ellipse_hyd <- ellipsoidhull(as.matrix(raw_morphdata_hyd))
area_hyd <- volume(ellipse_hyd) # 4418.76
ldensity1_hyd <- sumblen_hyd/area_hyd # 2.609177

# D1-ratio (raw)
ldensity1_AE/ldensity1_hyd # 3.754275

# Lineage density 2

# AE (raw)
rootblen_AE <- L_D2(tree_ultra_AE,raw_morphdata_AE) # 41.5791
area2_AE <- V_D2(ellipse_AE) # 1394.652
ldensity2_AE <- rootblen_AE/area2_AE # 0.02981325

# Hydrophis (raw)
rootblen_hyd <- L_D2(tree_ultra_hyd,raw_morphdata_hyd) # 107.3747
area2_hyd <- V_D2(ellipse_hyd) # 3332.275
ldensity2_hyd <- rootblen_hyd/area2_hyd # 0.03222264

# D2-ratio (raw)
ldensity2_AE/ldensity2_hyd # 0.9252267

#### analysis LOG10 ####

# Lineage density 1

# AE (LOG10)
sumblen_AElog10 <- L_D1(tree_ultra_AE,log10_morphdata_AE) # 0.7007288
ellipse_AElog10 <- ellipsoidhull(as.matrix(log10_morphdata_AE))
area_AElog10 <- volume(ellipse_AElog10) # 0.01818516
ldensity1_AElog10 <- sumblen_AElog10/area_AElog10 # 38.533

# Hydrophis (LOG10)
sumblen_hydlog10 <- L_D1(tree_ultra_hyd,log10_morphdata_hyd) # 5.143253
ellipse_hydlog10 <- ellipsoidhull(as.matrix(log10_morphdata_hyd))
area_hydlog10 <- volume(ellipse_hydlog10) # 0.238046
ldensity1_hydlog10 <- sumblen_hydlog10/area_hydlog10 # 21.60614

# D1-ratio (LOG10)
ldensity1_AElog10/ldensity1_hydlog10 # 1.783429

##############################################
#### simulations to test for significance ####
##############################################

# additional functions

glsfun <- function (tree, trait) {
  v <- vcv.phylo(tree, corr=TRUE)
  sv <- solve(v)
  ones <- rep(1, length(tree$tip.label))
  root.mean <- solve(ones %*% sv %*% ones) %*% ones %*% sv %*% trait
  sig2 <- 1/(length(tree$tip.label)-1)*(trait-root.mean) %*% sv %*% (trait-root.mean)
  list(root=root.mean, sig2=sig2)
}

# For use: split the array containing simulated tips into two objects 
# to correspond with AE and hyd clades/trees
## Use splitter2 function
splitter2 <- function(arr,r,c,e) {
  split.a <- arr[r,c,e]
}

#### data ####

# tree : loading nexus file containing 500 sample trees
trees500 <- read.nexus("./data/sea snake 500 trees.nex")

# morphological data
morphdata_AEhyd <- read.csv("./data/LOG10morphdata_AE-hyd.csv")
rownames(morphdata_AEhyd) <- morphdata_AEhyd$X
morphdata_AEhyd <- morphdata_AEhyd[-1]
raw_morphdata_AEhyd <- morphdata_AEhyd[,-c(3,4)]
log10_morphdata_AEhyd <- morphdata_AEhyd[,-c(1,2)]

#### prep trees ####

# drop tips non-AE and non-Hydrophis
trees_AEhyd <- drop.tip.multiPhylo(trees500, c("Hydrelaps_darwiniensis","Parahydrophis_mertoni","Ephalophis_greyi"))

# ultrametricise 500 pruned trees
trees_AEhyd_ultra <- lapply(trees_AEhyd, chronos)

# randomly sample 100,000 trees from 500 ultrametricised sample trees
trees_AEhyd_ultra_sample <- sample(trees_AEhyd_ultra, size = 100000, replace = TRUE)

#### raw morphdata simulation ####

# create container objects for for-loop
raw_obs_rndm <- array(dim = c(44,2,100000))
raw_root <- matrix(nrow=100000,ncol=2)
raw_sig2 <- matrix(nrow=100000,ncol=2)
raw_new_tips <- array(dim = c(44,2,100000))

# for loop (raw)
for (t in 1:length(trees_AEhyd_ultra_sample)){
  for (tip in 1:2){
    # obs_rndm: a 44 x 2 x 100,000 array (= 100,000 independent randomisations without replacement of 2 trait columns across 44 taxa)
    raw_obs_rndm[,tip,t] <- sample(raw_morphdata_AEhyd[,tip], replace = FALSE)
    
    # root: a 100,000 x 2 matrix containing the reconstructed trait values for each trait at the root of the 100,000 trees
    raw_root[t,tip] <- (glsfun(trees_AEhyd_ultra_sample[[t]],raw_obs_rndm[,tip,t]))[["root"]]
    
    # root: a 100,000 x 2 matrix containing the reconstructed trait values for each trait at the root of the 100,000 trees
    raw_sig2[t,tip] <- (glsfun(trees_AEhyd_ultra_sample[[t]],raw_obs_rndm[,tip,t]))[["sig2"]]
    
    # new tips: a 44 x 2 x 100,000 array (= 100,000 independent reconstructions of each trait value across 44 taxa using each iteration of root[t,tip])
    raw_new_tips[,tip,t] <- fastBM(trees_AEhyd_ultra_sample[[t]], a = raw_root[t,tip], raw_sig2 = raw_sig2[t,tip])
  }
}

# extract AE from ultrametricised sample trees
trees_AE_ultra_sample <- lapply(trees_AEhyd_ultra_sample, function(x) extract.clade(x,46))

# extract hyd from ultrametricised sample trees
trees_hyd_ultra_sample <- lapply(trees_AEhyd_ultra_sample, function(x) extract.clade(x,54))

# split for AE
AE_raw_new_tips <- splitter2(arr = raw_new_tips, r = c(1:9), c = c(1:2), e = c(1:100000))
# Complete array with colnames and rownames
colnames(AE_raw_new_tips) <- c("rel_girth","max_tot_len") # cols are the traits
for (m in dim(AE_raw_new_tips)[1]){
  rownames(AE_raw_new_tips) <- trees_AE_ultra_sample[[m]]$tip.label
}

# split for hyd
hyd_raw_new_tips <- splitter2(arr = raw_new_tips, r = c(10:44), c = c(1:2), e = c(1:100000))
# Complete array with colnames and rownames
colnames(hyd_raw_new_tips) <- c("rel_girth","max_tot_len")
for (g in dim(hyd_raw_new_tips)[1]){
  rownames(hyd_raw_new_tips) <- trees_hyd_ultra_sample[[g]]$tip.label
}

### null distribution lineage density 1 ratios

raw_D1sim_AE <- numeric(length(trees_AEhyd_ultra_sample))
raw_D1sim_hyd <- numeric(length(trees_AEhyd_ultra_sample))

### calculate raw_D1 for AE
for (o in 1:length(trees_AEhyd_ultra_sample)) {
  raw_D1sim_AE[o] <- L_D1(trees_AE_ultra_sample[[o]],AE_raw_new_tips[,,o]/V_D1(AE_raw_new_tips[,,o]))
  print(o-1+1)
}

### calculate raw_D1 for hyd
for (o in 1:length(trees_AEhyd_ultra_sample)) {
  raw_D1sim_hyd[o] <- L_D1(trees_hyd_ultra_sample[[o]],hyd_raw_new_tips[,,o]/V_D1(hyd_raw_new_tips[,,o]))
  print(o-1+1)
}

hist(raw_D1sim_AE/raw_D1sim_hyd)
abline(v=3.75,col="red",lwd=2)
(sum((raw_D1sim_AE/raw_D1sim_hyd)>=3.754275))/100000 # 0; observed is outside of simulated distribution


#### log10 morphdata simulation ####

# create container objects for for-loop
log10_obs_rndm <- array(dim = c(44,2,100000))
log10_root <- matrix(nrow=100000,ncol=2)
log10_sig2 <- matrix(nrow=100000,ncol=2)
log10_new_tips <- array(dim = c(44,2,100000))

# for loop (raw)
for (t in 1:length(trees_AEhyd_ultra_sample)){
  for (tip in 1:2){
    # obs_rndm: a 44 x 2 x 100,000 array (= 100,000 independent randomisations without replacement of 2 trait columns across 44 taxa)
    log10_obs_rndm[,tip,t] <- sample(log10_morphdata_AEhyd[,tip], replace = FALSE)
    
    # root: a 100,000 x 2 matrix containing the reconstructed trait values for each trait at the root of the 100,000 trees
    log10_root[t,tip] <- (glsfun(trees_AEhyd_ultra_sample[[t]],log10_obs_rndm[,tip,t]))[["root"]]
    
    # root: a 100,000 x 2 matrix containing the reconstructed trait values for each trait at the root of the 100,000 trees
    log10_sig2[t,tip] <- (glsfun(trees_AEhyd_ultra_sample[[t]],log10_obs_rndm[,tip,t]))[["sig2"]]
    
    # new tips: a 44 x 2 x 100,000 array (= 100,000 independent reconstructions of each trait value across 44 taxa using each iteration of root[t,tip])
    log10_new_tips[,tip,t] <- fastBM(trees_AEhyd_ultra_sample[[t]], a = log10_root[t,tip], log10_sig2 = log10_sig2[t,tip])
  }
}

# extract AE from ultrametricised sample trees
log10_trees_AE_ultra_sample <- lapply(trees_AEhyd_ultra_sample, function(x) extract.clade(x,46))

# extract hyd from ultrametricised sample trees
log10_trees_hyd_ultra_sample <- lapply(trees_AEhyd_ultra_sample, function(x) extract.clade(x,54))

# split for AE
AE_log10_new_tips <- splitter2(arr = log10_new_tips, r = c(1:9), c = c(1:2), e = c(1:100000))
# Complete array with colnames and rownames
colnames(AE_log10_new_tips) <- c("rel_girth","max_tot_len") # cols are the traits
for (m in dim(AE_log10_new_tips)[1]){
  rownames(AE_log10_new_tips) <- log10_trees_AE_ultra_sample[[m]]$tip.label
}

# split for hyd
hyd_log10_new_tips <- splitter2(arr = raw_new_tips, r = c(10:44), c = c(1:2), e = c(1:100000))
# Complete array with colnames and rownames
colnames(hyd_log10_new_tips) <- c("rel_girth","max_tot_len")
for (g in dim(hyd_log10_new_tips)[1]){
  rownames(hyd_log10_new_tips) <- log10_trees_hyd_ultra_sample[[g]]$tip.label
}

### null distribution lineage density 1 ratios

log10_D1sim_AE <- numeric(length(trees_AEhyd_ultra_sample))
log10_D1sim_hyd <- numeric(length(trees_AEhyd_ultra_sample))

### calculate log10_D1 for AE
for (o in 1:length(trees_AEhyd_ultra_sample)) {
  log10_D1sim_AE[o] <- L_D1(log10_trees_AE_ultra_sample[[o]],AE_log10_new_tips[,,o]/V_D1(AE_log10_new_tips[,,o]))
  print(o-1+1) # print progress
}

### calculate log10_D1 for hyd
for (o in 1:length(trees_AEhyd_ultra_sample)) {
  log10_D1sim_hyd[o] <- L_D1(log10_trees_hyd_ultra_sample[[o]],hyd_log10_new_tips[,,o]/V_D1(hyd_log10_new_tips[,,o]))
  print(o-1+1) # print progress
}

hist(log10_D1sim_AE/log10_D1sim_hyd)
abline(v=1.78,col="red",lwd=2)
(sum((log10_D1sim_AE/log10_D1sim_hyd)>=1.78))/100000 # 0.00268
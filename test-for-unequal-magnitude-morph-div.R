# R script for test of unequal magnitude of morphological diversification
# Authors: Vhon Oliver S. Garcia and Simone Blomberg
# Date: 4 February 2022

# Load packages
library(ape)
library(phytools)
library(cluster)
library(nlme)
library(MASS)

### Magnitude = unequal magnitudes of morphological change per phylogenetic branch ###

# function takes a phylo object for tree and a matrix for X
mean.pms3 <- function (tree, X) {
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
  mean(blen)
}

# Load tree (Using sea snake phylogeny from Sherratt et al. 2018)
tree <- read.tree("./data/sea snake MCC tree")
plot(tree, cex=0.5)
nodelabels(cex=0.6)

# Ultrametricise the tree (scaling the branch lengths such that basal node is t=1)
tree_ultra <- chronos(tree)
is.ultrametric(tree_ultra)
plot(tree_ultra,cex=0.5)

# Extract AE clade from tree based on node number
ext_AE_ultra <- extract.clade(tree_ultra,86)
plot(ext_AE_ultra)

# Load trait dataset with AE-only data
trait_data_AE <- read.csv("./data/trait_data_AE.csv") 
rownames(trait_data_AE) <- trait_data_AE$X
trait_data_AE <- trait_data_AE[,-1]
morph_trait_data_AE <- trait_data_AE[,-3]

# Plot the phylomorphospace
phylomorphospace(ext_AE_ultra,morph_trait_data_AE[,c(2,1)],label = "off")

# Extract Hydrophis clade
ext_hyd <- extract.clade(tree_ultra,51)
is.ultrametric(ext_hyd)
plot(ext_hyd, cex=0.5)

# Load trait dataset with hyd-only data
trait_data_hyd <- read.csv("./data/trait_data_hyd.csv")
rownames(trait_data_hyd) <- trait_data_hyd$X
trait_data_hyd <- trait_data_hyd[,-1]
morph_trait_data_hyd <- trait_data_hyd[,-3]

### Calculate mean morphometric branch lengths (observed) ###

# Aipysurus-Emydocephalus
mean.pms3(ext_AE_ultra,morph_trait_data_AE)
# Hydrophis
mean.pms3(ext_hyd,morph_trait_data_hyd)

# Calcualte observed M-ratio (observed)
(mean.pms3(ext_hyd,morph_trait_data_hyd))/(mean.pms3(ext_AE_ultra,morph_trait_data_AE))

##############################################
#### simulations to test for significance ####
##############################################

# updated code from Simone
glsfun <- function (tree, trait) {
  v <- vcv.phylo(tree, corr=TRUE)
  sv <- solve(v)
  ones <- rep(1, length(tree$tip.label))
  root.mean <- solve(ones %*% sv %*% ones) %*% ones %*% sv %*% trait
  sig2 <- 1/(length(tree$tip.label)-1)*(trait-root.mean) %*% sv %*% (trait-root.mean)
  list(root=root.mean, sig2=sig2)
}

#### Load tree data ####
# tree : loading nexus file containing 500 sample trees
trees500 <- read.nexus("./data/sea snake 500 trees.nex")

# Drop taxa not included in analyses (i.e. non-AE and non-Hydrophis)
aehyd_tree <- drop.tip.multiPhylo(trees500, c("Hydrelaps_darwiniensis","Parahydrophis_mertoni","Ephalophis_greyi"))

# Ultrametricise the 500 aehyd_tree
lapply_um_aehyd <- lapply(aehyd_tree, chronos)

# Randomly sample 100,000 trees from the 500 ultrametricised sample trees
um_aehyd_treesamp <- sample(lapply_um_aehyd, size = 100000, replace = TRUE)

### Load morphological trait data ###
trait <- read.csv("./data/trait_data_AE-hyd.csv")
rownames(trait) <- trait$X
trait <- trait[,-1]
aehyd_morph <- trait[,-3]

### Simulation of `new tip trait values` under Brownian Motion model of evolution ###

# create container objects for for-loop
obs_rndm <- array(dim = c(44,2,100000))
root <- matrix(nrow=100000,ncol=2)
sig2 <- matrix(nrow=100000,ncol=2)
new_tips <- array(dim = c(44,2,100000))

# for loop
for (t in 1:length(um_aehyd_treesamp)){
  for (tip in 1:2){
    # obs_rndm: a 44 x 2 x 100,000 array (= 100,000 independent randomisations without replacement of 2 trait columns across 44 taxa)
    obs_rndm[,tip,t] <- sample(aehyd_morph[,tip], replace = FALSE)
    
    # root: a 100,000 x 2 matrix containing the reconstructed trait values for each trait at the root of the 100,000 trees
    root[t,tip] <- (glsfun(um_aehyd_treesamp[[t]],obs_rndm[,tip,t]))[["root"]]
    
    # root: a 100,000 x 2 matrix containing the reconstructed trait values for each trait at the root of the 100,000 trees
    sig2[t,tip] <- (glsfun(um_aehyd_treesamp[[t]],obs_rndm[,tip,t]))[["sig2"]]
    
    # new tips: a 44 x 2 x 100,000 array (= 100,000 independent reconstructions of each trait value across 44 taxa using each iteration of root[t,tip])
    new_tips[,tip,t] <- fastBM(um_aehyd_treesamp[[t]], a = root[t,tip], sig2 = sig2[t,tip])
  }
}

### Prepare input objects for analyses ###

# Extract AE clade from ultrametricised trees
lapply_ext_ae <- lapply(um_aehyd_treesamp, function(x) extract.clade(x,46)) # 46 is node leading to AE clade

# Extract hyd clade from ultrametricised trees
lapply_ext_hyd <- lapply(um_aehyd_treesamp, function(x) extract.clade(x,54)) # 54 is node leading to hyd clade

# Split the array containing simulated tips into two objects to correspond with AE and hyd clades/trees

## Use splitter2 function
splitter2 <- function(arr,r,c,e) {
  split.a <- arr[r,c,e]
}

### Split for AE
ae_new_tips <- splitter2(arr = new_tips, r = c(1:9), c = c(1:2), e = c(1:100000))

#### Complete array with colnames and rownames
colnames(ae_new_tips) <- c("rel_girth","max_tot_len") # cols are the traits
for (m in dim(ae_new_tips)[1]){
  rownames(ae_new_tips) <- lapply_ext_ae[[m]]$tip.label
}

### Split for hyd
hyd_new_tips <- splitter2(arr = new_tips, r = c(10:44), c = c(1:2), e = c(1:100000))

#### Complete array with colnames and rownames
colnames(hyd_new_tips) <- c("rel_girth","max_tot_len")
for (g in dim(hyd_new_tips)[1]){
  rownames(hyd_new_tips) <- lapply_ext_hyd[[g]]$tip.label
}

Msim_ae <- vector("list",length(um_aehyd_treesamp))
Msim_hyd <- vector("list",length(um_aehyd_treesamp))

for (o in 1:length(um_aehyd_treesamp)){
  Msim_ae[o] <- mean.pms3(tree = lapply_ext_ae[[o]], ae_new_tips[,,o])
}

for (o in 1:length(um_aehyd_treesamp)){
  Msim_hyd[o] <- mean.pms3(tree = lapply_ext_hyd[[o]], hyd_new_tips[,,o])
}

### [[PLOT]] NULL DISTRIBUTION M-RATIOS ###

hist((unlist(Msim_hyd))/(unlist(Msim_ae)))
abline(v=1.57, col = "red", lwd = 2)
sum((unlist(Msim_hyd))/(unlist(Msim_ae))>=1.57)/100000 # 0.10371

##########################
# Log10 transformed data
##########################

### Magnitude = unequal magnitudes of morphological change per phylogenetic branch ###

# function takes a phylo object for tree and a matrix for X
mean.pms3 <- function (tree, X) {
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
  mean(blen)
}

# Load tree (Using sea snake phylogeny from Sherratt et al. 2018)
tree <- read.tree("./data/sea snake MCC tree")
plot(tree, cex=0.5)
nodelabels(cex=0.6)

# Ultrametricise the tree (scaling the branch lengths such that basal node is t=1)
tree_ultra <- chronos(tree)
is.ultrametric(tree_ultra)
plot(tree_ultra,cex=0.5)

# Extract AE clade from tree based on node number
ext_AE_ultra <- extract.clade(tree_ultra,86)
plot(ext_AE_ultra)

# Load trait dataset with AE-only data
log10morph_data_AE <- read.csv("./data/LOG10morphdata_AE.csv") 
rownames(log10morph_data_AE) <- log10morph_data_AE$X
log10morph_data_AE <- log10morph_data_AE[,-1]
log10morph_data_AE <- log10morph_data_AE[,-c(1,2)]

# Plot the phylomorphospace
phylomorphospace(ext_AE_ultra,log10morph_data_AE[,c(2,1)],label = "off")

# Extract Hydrophis clade
ext_hyd <- extract.clade(tree_ultra,51)
is.ultrametric(ext_hyd)
plot(ext_hyd, cex=0.5)

# Load trait dataset with hyd-only data
log10morph_data_hyd <- read.csv("./data/LOG10morphdata_hyd.csv")
rownames(log10morph_data_hyd) <- log10morph_data_hyd$X
log10morph_data_hyd <- log10morph_data_hyd[,-1]
log10morph_data_hyd <- log10morph_data_hyd[,-c(1,2)]

phylomorphospace(ext_hyd,log10morph_data_hyd[,c(2,1)],label = "off")

### Calculate mean morphometric branch lengths (observed) ###

# Aipysurus-Emydocephalus
mean.pms3(ext_AE_ultra,log10morph_data_AE) # 0.04379555
# Hydrophis
mean.pms3(ext_hyd,log10morph_data_hyd) # 0.07563608

# Calcualte observed M-ratio (observed)
(mean.pms3(ext_hyd,log10morph_data_hyd))/(mean.pms3(ext_AE_ultra,log10morph_data_AE)) # 1.727027

##############################################
#### simulations to test for significance ####
##############################################

# updated code from Simone
glsfun <- function (tree, trait) {
  v <- vcv.phylo(tree, corr=TRUE)
  sv <- solve(v)
  ones <- rep(1, length(tree$tip.label))
  root.mean <- solve(ones %*% sv %*% ones) %*% ones %*% sv %*% trait
  sig2 <- 1/(length(tree$tip.label)-1)*(trait-root.mean) %*% sv %*% (trait-root.mean)
  list(root=root.mean, sig2=sig2)
}

#### Load tree data ####
# tree : loading nexus file containing 500 sample trees
trees500 <- read.nexus("./data/sea snake 500 trees.nex")

# Drop taxa not included in analyses (i.e. non-AE and non-Hydrophis)
aehyd_tree <- drop.tip.multiPhylo(trees500, c("Hydrelaps_darwiniensis","Parahydrophis_mertoni","Ephalophis_greyi"))

# Ultrametricise the 500 aehyd_tree
lapply_um_aehyd <- lapply(aehyd_tree, chronos)

# Randomly sample 100,000 trees from the 500 ultrametricised sample trees
um_aehyd_treesamp <- sample(lapply_um_aehyd, size = 100000, replace = TRUE)

### Load morphological trait data ###
log10trait <- read.csv("./data/LOG10morphdata_AE-hyd.csv")
rownames(log10trait) <- log10trait$X
log10trait <- log10trait[,-1]
log10_aehyd_morph <- log10trait[,-c(1,2)]

### Simulation of `new tip trait values` under Brownian Motion model of evolution ###

# create container objects for for-loop
log10M_obs_rndm <- array(dim = c(44,2,100000))
log10M_root <- matrix(nrow=100000,ncol=2)
log10M_sig2 <- matrix(nrow=100000,ncol=2)
log10M_new_tips <- array(dim = c(44,2,100000))

# for loop
for (t in 1:length(um_aehyd_treesamp)){
  for (tip in 1:2){
    # obs_rndm: a 44 x 2 x 100,000 array (= 100,000 independent randomisations without replacement of 2 trait columns across 44 taxa)
    log10M_obs_rndm[,tip,t] <- sample(log10_aehyd_morph[,tip], replace = FALSE)
    
    # root: a 100,000 x 2 matrix containing the reconstructed trait values for each trait at the root of the 100,000 trees
    log10M_root[t,tip] <- (glsfun(um_aehyd_treesamp[[t]],log10M_obs_rndm[,tip,t]))[["root"]]
    
    # root: a 100,000 x 2 matrix containing the reconstructed trait values for each trait at the root of the 100,000 trees
    log10M_sig2[t,tip] <- (glsfun(um_aehyd_treesamp[[t]],log10M_obs_rndm[,tip,t]))[["sig2"]]
    
    # new tips: a 44 x 2 x 100,000 array (= 100,000 independent reconstructions of each trait value across 44 taxa using each iteration of root[t,tip])
    log10M_new_tips[,tip,t] <- fastBM(um_aehyd_treesamp[[t]], a = log10M_root[t,tip], log10M_sig2 = log10M_sig2[t,tip])
    print(t-1+1)
  }
}

### Prepare input objects for analyses ###

# Extract AE clade from ultrametricised trees
log10M_lapply_ext_ae <- lapply(um_aehyd_treesamp, function(x) extract.clade(x,46)) # 46 is node leading to AE clade

# Extract hyd clade from ultrametricised trees
log10M_lapply_ext_hyd <- lapply(um_aehyd_treesamp, function(x) extract.clade(x,54)) # 54 is node leading to hyd clade

# Split the array containing simulated tips into two objects to correspond with AE and hyd clades/trees

## Use splitter2 function
splitter2 <- function(arr,r,c,e) {
  split.a <- arr[r,c,e]
}

### Split for AE
log10M_ae_new_tips <- splitter2(arr = log10M_new_tips, r = c(1:9), c = c(1:2), e = c(1:100000))

#### Complete array with colnames and rownames
colnames(log10M_ae_new_tips) <- c("rel_girth","max_tot_len") # cols are the traits
for (m in dim(log10M_ae_new_tips)[1]){
  rownames(log10M_ae_new_tips) <- lapply_ext_ae[[m]]$tip.label
}

### Split for hyd
log10M_hyd_new_tips <- splitter2(arr = log10M_new_tips, r = c(10:44), c = c(1:2), e = c(1:100000))

#### Complete array with colnames and rownames
colnames(log10M_hyd_new_tips) <- c("rel_girth","max_tot_len")
for (g in dim(log10M_hyd_new_tips)[1]){
  rownames(log10M_hyd_new_tips) <- lapply_ext_hyd[[g]]$tip.label
}

log10M_sim_ae <- vector("list",length(um_aehyd_treesamp))
log10M_sim_hyd <- vector("list",length(um_aehyd_treesamp))

for (o in 1:length(um_aehyd_treesamp)){
  log10M_sim_ae[o] <- mean.pms3(tree = lapply_ext_ae[[o]], log10M_ae_new_tips[,,o])
  print(o-1+1)
}

for (o in 1:length(um_aehyd_treesamp)){
  log10M_sim_hyd[o] <- mean.pms3(tree = lapply_ext_hyd[[o]], log10M_hyd_new_tips[,,o])
  print(o-1+1)
}

### [[PLOT]] NULL DISTRIBUTION M-RATIOS ###

hist((unlist(log10M_sim_hyd))/(unlist(log10M_sim_ae)))
abline(v=1.72, col = "red", lwd = 2)
sum((unlist(Msim_hyd))/(unlist(Msim_ae))>=1.72)/100000 # 0.06634


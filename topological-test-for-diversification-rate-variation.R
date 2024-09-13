# R script for testing significance of observed clade size disparity
# Author: Vhon Oliver S. Garcia
# Date: 11 January 2022

# Load packages
library(ape)
library(TreeSim)

# Simulate trees under equal-rates Markov model using `sim.bd.taxa()`
# ntaxa = number of total tips (i.e. total taxa of two clades being compared)
# numsimtree = number of simulated trees to generate null distribution
dat <- sim.bd.taxa(n = ntaxa, numbsim = numsimtree, lambda = 1, mu = 0)

# Wrapper function
disparity.multiPhylo <- function(dat, low_sp) {
  # dat = object containing simulated trees
  # low_sp = lower number of the clade split (e.g. 9 for a 9/34 clade split)
  message("-Reading trees block...")
  trees_block <- dat
  message("-Trees block read succesfully...")
  message("-Performing balance on trees block...")
  balance <- lapply(trees_block, balance)
  message("-Balance function complete...")
  balance_specific <- lapply(balance, function(x) x[1,1:2])
  bal_spec_unlist <- unlist(balance_specific)
  bal_spec_df <- as.data.frame(matrix(bal_spec_unlist, ncol = 2, byrow = T))
  bal_spec_min <- apply(bal_spec_df, 1, FUN = min)
  bal_spec_min_df <- as.data.frame(matrix(bal_spec_min, ncol = 1, byrow = T))
  message("Complete!")
  print(sum(bal_spec_min_df$V1<=low_sp)/length(bal_spec_min_df$V1))
}

# Example using Pyron et al. 2013 (split: 6/21)
pyron <- sim.bd.taxa(n = 27, numbsim = 100000, lambda = 1, mu = 0)
disparity.multiPhylo(pyron)

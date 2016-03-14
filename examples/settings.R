rm(list = setdiff(ls(), lsf.str()))

library(DatABEL)

n   <- 2000         # number of individuals
l   <- 3            # number of covariates
int <- 3
r   <- 2            # nr of genetic predictors
nr_snps <- 100000   # number of snps
m   <- r *  nr_snps # total nr of genetic data columns
t   <- 10000        # number of traits
var <- 0.05

set.seed(1001)

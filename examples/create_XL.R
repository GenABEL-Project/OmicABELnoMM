source("settings.R")

## Remove any existing output files
fn <- "XL.fvi"
if (file.exists(fn)) file.remove(fn)
fn <- "XL.fvd"
if (file.exists(fn)) file.remove(fn)

## Create the new files
## dim(XL) = (nr_covar + 1) * nr_indiv
XL <- matrix(rnorm((l+1)*n), ncol=(l+1))

## Add some NAs
for (i in 1:(n*(l+1))) {
    if(sample(1:100, 1) > 95){
        XL[i] <- 0/0
    }
}

## First column should be ones (intercept)
for (i in 1:n) {
    XL[i] <- 1
}

colnames(XL) <- c("intercept", paste0("cov",1:l))
rownames(XL) <- paste0("ind", 1:n)
XL_db <- matrix2databel(XL, filename="XL", type="FLOAT")

#XL[1:n,1:(l+1)]
#XL

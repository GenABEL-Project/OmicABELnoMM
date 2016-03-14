source("settings.R")

## Remove any existing output files
fn <- "XR.fvi"
if (file.exists(fn)) file.remove(fn)
fn <- "XR.fvd"
if (file.exists(fn)) file.remove(fn)

## Load the XL data
XL <- as(databel("XL.fvi"), "matrix")

## Load the Y data
Y <- as(databel("Y.fvi"), "matrix")

## Create the new files
## dim(XR) = (nr_gen_predictors * nr_snps) * nr_indiv
XR <- matrix(rnorm(m*n), ncol=m)

## LCK: commented the following code because it doesn't work for the
## current settings (the yIdx variable gets too large) and I'm not
## sure what Alvaro was trying to do here.
## for (i in 1 + r * (0:((m-2)/r)) ) {
##     ##print(i)
##     yIdx <- ceiling(i/r)
##     ##print(i)
##     ##print(yIdx)
##     for (j in 1:n)
##     {
##         ## Set the SNP value for individual j, SNP i to
##         XR[j, i] <- Y[j, yIdx]
##         ## Introduce a bit of correlation with the covariates into
##         ## XR[j, i], for the first genetic predictor (r=1)
##         for (k in 1:l)
##         {
##             XR[j, i] <- XR[j, i] - XL[j, k] * 0.01
##         }
##         ## Add correlation with the other r columns
##         for (k in 1:(r-1))
##         {
##             XR[j, i] <- XR[j, i] - XR[j, i+k] * 0.01
##         }
##         ##XR[j,i]=XR[j,i]/2.8888
##         ##XR[j,i] = XR[j,i]*runif(1, 1.0-var, 1.0)
##     }
##}

## Add some NAs
for (i in 1:(n*m)){
    if(sample(1:100,1) > 90) {
        XR[i] <- 0/0
    }
}

colnames(XR) <- paste0("miss", 1:m)
cnames <- vector(m)
for (i in 1:(m/r)) {
    for (j in 1:r) {
        cnames[(i-1)*r+(j)] <- paste0("snp",
                                            paste(i, j, sep="_")
                                            )
    }
}
colnames(XR) <- cnames

rownames(XR) <- paste0("ind", 1:n)
XR_db <- matrix2databel(XR, filename="XR", type="FLOAT")

source("settings.R")

## Remove any existing output files
fn <- "Y.fvi"
if (file.exists(fn)) file.remove(fn)
fn <- "Y.fvd"
if (file.exists(fn)) file.remove(fn)

## Create the new files
## dim(Y) = traits * nr_indiv
Y <- matrix(rnorm(t*n), ncol=t)

## Add some NAs
for (i in 1:(t*n)) {
    if(sample(1:100,1) > 90) {
        Y[i] <- 0/0
    }
}

colnames(Y) <- paste0("ph", 1:t)
rownames(Y) <- paste0("ind", 1:n)
Y_db <- matrix2databel(Y, filename="Y", type="FLOAT")

#Y[1:n,1:t]
#Y

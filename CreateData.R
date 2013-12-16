library(DatABEL)

n = 10 # number of individuals
l = 2    # number of covariates
m = 10 # number of snps
t = 5  # number of traits

set.seed(1001)
runif(3)

XL <- matrix(rnorm((l+1)*n),ncol=(l+1)) # first column should be ones (intercept)
for(i in 1:n*l){ if(sample(1:10,1) > 5) XL[i]=0/0}

colnames(XL) <- c("intercept", paste("cov",1:l,sep=""))
rownames(XL) <- paste("ind",1:n,sep="")
XL_db <- matrix2databel(XL,filename="XL",type="FLOAT")

#XL[1:n,1:(l+1)]
XL

#prob= TRUE or FALSE
prob=FALSE
if (!prob){
	XR <- matrix(rnorm(m*n),ncol=m)
	for(i in 1:n*m){ if(sample(1:10,1) > 5) XR[i]=0/0}

	colnames(XR) <- paste("snp",1:m,sep="")
	rownames(XR) <- paste("ind",1:n,sep="")
	XR_db <- matrix2databel(XR,filename="XR",type="FLOAT")
} else{
	XR <- matrix(runif(2*m*n),ncol=2*m)
	for(i in 1:n*m){ if(sample(1:10,1) > 5) XR[i]=0/0}
	colnames(XR) <- paste("snp",(1:(2*m)),sep="")
	colnames(XR)[seq(from=1,to=2*m,by=2)] <- paste("snp",1:m,"_11",sep="")
	colnames(XR)[seq(from=2,to=2*m,by=2)] <- paste("snp",1:m,"_01",sep="")
	rownames(XR) <- paste("ind",1:n,sep="")
	tmp11=XR[,seq(from=1,to=2*m,by=2)]
	tmp01_00=1-tmp11
	tmp=matrix(runif(m*n),ncol=m)
	tmp01=tmp01_00*tmp
	XR[,seq(from=2,to=2*m,by=2)]=tmp01
	XR_db <- matrix2databel(XR,filename="XR",type="FLOAT")
}

#XR[1:n,1:m]
XR

Y <- matrix(rnorm(t*n),ncol=t)
for(i in 1:n*t){ if(sample(1:10,1) > 5) Y[i]=0/0}
colnames(Y) <- paste("ph",1:t,sep="")
rownames(Y) <- paste("ind",1:n,sep="")
Y_db <- matrix2databel(Y,filename="Y",type="FLOAT")

#Y[1:n,1:t]
Y

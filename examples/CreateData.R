rm(list = setdiff(ls(), lsf.str()))

library(DatABEL)

n = 1000 # number of individuals
l = 4    # number of covariates+1 for intercept
r = 2
m = r*1000 # number of snps
t = 1000  # number of traits


set.seed(1001)
runif(3)

XL <- matrix(rnorm((l+1)*n),ncol=(l+1)) # first column should be ones (intercept)
for(i in 1:(n*(l+1))){ if(sample(1:100,1) > 95){XL[i]=0/0} }
for(i in 1:n){ XL[i]=1}

colnames(XL) <- c("intercept", paste("cov",1:l,sep=""))
rownames(XL) <- paste("ind",1:n,sep="")
XL_db <- matrix2databel(XL,filename="XL",type="FLOAT")

#XL[1:n,1:(l+1)]
#XL

Y <- matrix(rnorm(t*n),ncol=t)
for(i in 1:(n*t)){ if(sample(1:100,1) > 85) {Y[i]=0/0} }
colnames(Y) <- paste("ph",1:t,sep="")
rownames(Y) <- paste("ind",1:n,sep="")
Y_db <- matrix2databel(Y,filename="Y",type="FLOAT")

#Y[1:n,1:t]
#Y

#prob= TRUE or FALSE
prob=FALSE
if (!prob){
	XR <- matrix(rnorm(m*n),ncol=m)

	
	for(i in 1 + r*(0:((m-2)/r)) )
	{ 
		yIdx=sample(1:t, 1)
		#print(i)
		#print(yIdx)
		for(j in 1:n)
		{ 
			XR[j,i]=Y[j,yIdx]
			for(k in 1:l)
			{
				XR[j,i]=XR[j,i]-XL[j,k]*0.01
			}
			for(k in 1:(r-1))
			{
				XR[j,i]=XR[j,i]-XR[j,i+k]*0.01
			} 
			#XR[j,i]=XR[j,i]/2.8888
		}
	}

	for(i in 1:(n*m)){ if(sample(1:100,1) > 85) XR[i]=0/0}

	

	colnames(XR) <- paste("snp",1:m,sep="")
	rownames(XR) <- paste("ind",1:n,sep="")
	XR_db <- matrix2databel(XR,filename="XR",type="FLOAT")
} else{
	XR <- matrix(runif(2*m*n),ncol=2*m)
	for(i in 1:(n*m)){ if(sample(1:100,1) > 85) {XR[i]=0/0}}
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
#XR

XR[is.nan(XR)] <- 0
for(i in 1:r )
{ 
	Avg=sum(XR[,i])/(n-sum(0==XR[,i]))
	XR[,i][ XR[,i] == 0 ] <- Avg
}


b=c()
for(k in 1:t )
{
	for(i in 1 + r*(0:((m-2)/r)) )
	{ 
		Xtemp=cbind(XL,XR[,i])
		for(j in 1:(r-1))
		{
			Xtemp=cbind(Xtemp,XR[,i+j])
		}
		
		to_remove=c()
		idx=1
		for(j in 1:n)
		{
			if(sum(is.nan(XL[j,])) > 0 || is.nan(Y[j,k]))
			{
				to_remove[idx]=j
				idx=idx+1
			}
		}
		X=Xtemp[-c(to_remove),]
		y=Y[-c(to_remove),k]

		S=base::t(X)%*%X
		Xy = base::t(X)%*%y

		btemp=solve(S,Xy)
		#print(y)
		#print(X)
		#print(btemp)
		b=c(b,btemp)
		

		

	} 
}

to.write = file("bpre.fvd", "wb")
writeBin(b, to.write,size=4)
close(to.write)




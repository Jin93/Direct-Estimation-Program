##### Realdata_Direct.R
====================================================================================
Description: Modified Estimation and analysis of Differential Networks of FASD data 
             by directly using FASD data matrix
====================================================================================

====================================================================================
####### 1. dpmchangelanbdaloss.R #######
## based on Dave Zhao's "dpm.R" package
====================================================================================


dpm <- function(X1,X0,
                lambda=NULL,nlambda=100,lambda.min.ratio=NULL,
                rho=NULL,shrink=NULL,prec=0.0001,max.ite=100,
                correlation=FALSE,perturb=FALSE,
                tuning=c("none","aic","bic","cv"),folds=5)
{
    ###calculate covariance matrices, calculate lambdas, etc###
    if(ncol(X1)!=ncol(X0))
        {
            cat("X1 and X0 need to have the same number of columns.\n");
            return(NULL);
        }
    n1 <- nrow(X1); n0 <- nrow(X0);
    ## the number of parameters is p(p+1)/2
    p <- ncol(X1); d <- p*(p+1)/2;
    maxdf <- max(n1,n0,d);
    
    ## construct kronecker product, first perturb individual
    ## matrices
    if(correlation){ S1 <- cor(X1); S0 <- cor(X0); } else
        { S1 <- cov(X1)*(1-1/n1); S0 <- cov(X0)*(1-1/n0); }
    if(is.logical(perturb))
        {
            if(perturb)
                {
                    ## same perturbation as the clime software
                    eigvals1 <- eigen(S1,only.values=TRUE)$values;
                    eigvals0 <- eigen(S0,only.values=TRUE)$values;
                    perturb1 <- max(max(eigvals1)-p*min(eigvals1),0)/(p-1);
                    perturb0 <- max(max(eigvals0)-p*min(eigvals0),0)/(p-1);
                } else { perturb1 <- 0; perturb0 <- 0; }
        }
    S <- kronecker(S0+diag(p)*perturb0,S1+diag(p)*perturb1);

    gc(reset=TRUE);
    
    if(is.null(rho)){ rho <- sqrt(d); }
    if(is.null(shrink)){ shrink <- 1.5; }
    
    ## create smaller problem with p(p+1)/2 rows and columns
    ind <- matrix(1:p^2,nrow=p);
    lowtri <- lower.tri(ind,diag=TRUE);
    ## sum appropriate columns
    S[,ind[lowtri]] <- S[,ind[lowtri]]+S[,t(ind)[lowtri]];
    S[,diag(ind)] <- S[,diag(ind)]/2;
    S <- S[,ind[lowtri]];
    ## sum appropriate rows
    S[ind[lowtri],] <- S[ind[lowtri],]+S[t(ind)[lowtri],];
    S[diag(ind),] <- S[diag(ind),]/2;
    S <- S[ind[lowtri],];
    ## multiply entries of e by 2 to be consistent.
    e <- as.vector(S1-S0);
    e[ind[lowtri]] <- 2*e[ind[lowtri]];
    e[diag(ind)] <- e[diag(ind)]/2;
    e <- e[ind[lowtri]];

    gc(reset=TRUE);
    
    ## need to pass indices of the diagonal elements to the C
    ## fxn so that it can divide the corresponding lambdas by 2
    diags <- diag(p)[lowtri];
    nlambda=100;########################################  nlambda=100
    ## determine lambda based on S and e
    if(!is.null(lambda)){ nlambda <- length(lambda); }
    if(is.null(lambda))
        {
            if(is.null(lambda.min.ratio)){ lambda.min.ratio <- 0.004; }
            ##lambda.max <- min(max(S-diag(diag(S))),
            ##                -min(S-diag(diag(S))));
            lambda.max <- max(abs(e));
            lambda.min <- lambda.min.ratio*lambda.max;
            lambda <- exp(seq(log(lambda.max),log(lambda.min),
                              length=nlambda));
        }

    ## call C function
    gamma <- 1/rho;
    lambda <- lambda - shrink*prec;
    nlambda = length(lambda);
    ret = vector("list", nlambda);
    diff = .C("dpm", as.double(S), diff=double(d*nlambda),
        as.integer(d), as.integer(diags),as.double(lambda),
        as.integer(nlambda), as.double(gamma), as.integer(max.ite),
        as.double(prec), as.double(e))$diff;
    ## DON'T PERTURB S1-S0
    for(i in 1:nlambda)
    {
        ret[[i]] <- matrix(NA,nrow=p,ncol=p);
        ret[[i]][ind[lowtri]] <- diff[((i-1)*d+1):(i*d)];
        ret[[i]][t(ind)[lowtri]] <- diff[((i-1)*d+1):(i*d)];
    }
    
    rm(list=c("S","diff","ind","lowtri")); gc(reset=TRUE);
    
    ## ==============================================================
    ## run tuning
    ## ==============================================================
    opt <- switch(tuning[1], ## default is "none"
                  none=NA,
                  cv=dpm.cv(X1,X0,
                      lambda,
                      rho,shrink,prec,max.ite,
                      correlation,perturb,
                      folds),
                  aic=dpm.ic(S1,S0,ret,n1+n0,2),
                  bic=dpm.ic(S1,S0,ret,n1+n0,log(n1+n0)));
    if(!is.na(opt[1])){ names(opt) <- c("max","1","L1","sp","F","nc"); }
    
    return(list(dpm=ret,lambda=lambda,nlambda=nlambda,opt=opt));
}

## **************************************************************
## function to calculate loss
## return a vector of different types of losses
## D=estimated difference matrix
## S1,S0=true, validation set, or training set covariances
## **************************************************************
loss <- function(D,S1,S0)
{
    err <- S1%*%D%*%S0-S1+S0;
    return(c(max(abs(err)), ## max
             sum(abs(err)), ## l1, element-wise
             max(apply(err,1,function(r){ sum(abs(r)) })), ## matrix L1
             svd(err,nu=0,nv=0)$d[1], ## spectral
             sqrt(sum(err^2)), ## frobenius
             sum(svd(err,nu=0,nv=0)$d))); ## nuclear
}

## **************************************************************
## tuning methods
## **************************************************************
## ==============================================================
## cv
## ==============================================================
dpm.cv <- function(X1,X0,
                   lambda=NULL,
                   rho=NULL,shrink=NULL,prec=0.0001,max.ite=100,
                   correlation=FALSE,perturb=FALSE,
                   folds=5)
{
  if(ncol(X1)!=ncol(X0))
  {
    cat("X1 and X0 need to have the same number of columns.\n");
    return(NULL);
  }
  n1 <- nrow(X1); n0 <- nrow(X0);
  p <- ncol(X1); d <- p*(p+1)/2;
  if(is.null(rho)){ rho <- sqrt(d); }
  if(is.null(shrink)){ shrink <- 1.5; }
  
  ind1 <- sample(1:n1,n1,replace=FALSE);
  ind0 <- sample(1:n0,n0,replace=FALSE);
  losses <- array(NA,c(folds,6,length(lambda)));
  cat("CV fold:");
  for(i in 1:folds)
  {
    cat("",i);
    test1 <- ind1[((i-1)*n1/folds+1):(i*n1/folds)];
    test0 <- ind0[((i-1)*n0/folds+1):(i*n0/folds)];
    fit <- dpm(X1[-test1,],X0[-test0,],lambda,lenth(lambda),NULL,
               rho,shrink,prec,max.ite,correlation,perturb,tuning="none");
    ## don't need to perturb the estimated cov/cor from test set
    if(correlation)
    { S1.test <- cor(X1[test1,]); S0.test <- cor(X0[test0,]); } else
    {
      S1.test <- cov(X1[test1,])*(1-1/length(test1));
      S0.test <- cov(X0[test0,])*(1-1/length(test0));
    }

    losses[i,,] <- sapply(fit$dpm,loss,S1.test,S0.test);
  }
  cat("\n");
  opt <- apply(apply(losses,c(2,3),mean),1,which.min);
  return(opt);
}

## ==============================================================
## ic (information criteria)
## ==============================================================
dpm.ic <- function(S1,S0,ret,n,penalty)
{
    lowtri <- which(lower.tri(ret[[1]],diag=TRUE));
    df <- sapply(ret,function(x){ sum(x[lowtri]!=0); });
    ic <- scale(n*sapply(ret,loss,S1,S0),
                center=-penalty*df,scale=FALSE);
    return(apply(ic,1,which.min));
}


====================================================================================
####### 2. directgenerate.R #######
## DEA.R on Real Data
====================================================================================
## **************************************************************
## DEA.R on Real Data
## **************************************************************

rm(list=ls());
dyn.load("dpm.so");
source("dpmchangelambda.R");
library(MASS);
library(ROCR);


n=180;
n1=n*25;
n0=n*30;
## **************************************************************
## Read Data
## **************************************************************

data=read.csv("/home/merganser/jinjin/data1.csv",header=T)
group=read.csv("/home/merganser/jinjin/group.csv",header=T)
data=as.matrix(data)
group=as.matrix(group)
ns=55
nt=180
p=ncol(data)
nlambda=100;
 
ret = vector("list", ns);
dd=matrix(0,nrow(data),ncol(data))
for (i in 1:ns)
{
ret[[i]]=matrix(NA,nrow=nt,ncol=p)
ret[[i]]=data[(nt*(i-1)+1):(nt*i),]

for (j in 1:nt)
{for (k in 3:p)
ret[[i]][j,k]=ret[[i]][j,k]-mean(ret[[i]][,k]);
}
}
############## Or Standardize by subtracting each entries of FASD matrix by the total mean #########
#m=c(1:p);
#for (j in 1:nt)
#{for (k in 3:p)
#{m[k]=mean(data[,k]);
#ret[[i]][j,k]=ret[[i]][j,k]-m[k];
#}
#}
#}

numcontrol=0;
for (i in 1:ns)
{ if (group[i,4]==1)
{numcontrol=numcontrol+1}
}
numFASD=0;
for (i in 1:ns)
{ if (group[i,4]==2)
{numFASD=numFASD+1}
}

X1.t=matrix(0,numcontrol*nt,p-2)
X0.t=matrix(0,numFASD*nt,p-2)

i=1;
for (j in 1:ns)
{
if (group[j,4]==1)
{
for (k in 1:ns)
{ if(ret[[k]][1,1]==group[j,1])
{
X1.t[(nt*(i-1)+1):(nt*i),]=ret[[k]][,(3):(p)];
i=i+1;
}
}
}
}

i=1;
for (j in 1:ns)
{
if (group[j,4]==2)
{
for (k in 1:ns)
{ if(ret[[k]][1,1]==group[j,1])
{
X0.t[(nt*(i-1)+1):(nt*i),]=ret[[k]][,(3):(p)];
i=i+1;
}
}
}
}
#####Calculation#####
setwd("/home/merganser/jinjin")
fit.aic <- dpm(X1.t,X0.t,nlambda=100,tuning="aic");
fit.cv <- dpm(X1.t,X0.t,nlambda=100,tuning="cv",folds=3);
setwd("/home/merganser/jinjin/TestDE/Direct/Result/Direct100")
save(fit.aic,fit.cv, file="directlambda.RData")

====================================================================================
####### 3. directanalysis.R #######
====================================================================================

library(MASS);
p=74;
n=180;
n1=n*25;
n0=n*30;
nlambda=10;
q=fit.aic[[4]][[5]]
th=1e-05;

setwd("/home/merganser/jinjin/TestDE/direct/Result")
load("directlambda.RData")

m=fit.aic[[1]][[q]];
x=matrix(0,1,p*p);
y=matrix(0,1,p*p);
k=1;
for (i in (1:(p)))
{for (j in (1:(i)))
{ if (m[i,j]>=th) {x[1,k]=j;y[1,k]=i;k=k+1;}
}
}
pdf("plotNonzeroDE.pdf")
plot(x,y,type = "p",cex=.5)
dev.off()

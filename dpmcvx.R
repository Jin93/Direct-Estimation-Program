##### dpmcvx.R
========================================================================================
Description: Modified Estimation and analysis of Differential Networks of Simulation 
             based on FASD Data Calling CVX package from Matlab in R and example 
             usage. R code was modified based on code from the dpm.R 
             (https://github.com/sdzhao/dpm/blob/master/dpm.R).
             Use SICE as true precision matrices to generate data and compare dpmcvx.R
             and dpm.R
========================================================================================

========================================================================================
####### 1. dpmcvxgenerate.R #######
## To generate data with precision matrices obtained by SICE simulation and use dpmcvx.R
   and dpm.R to compare the results. 
========================================================================================

setwd("/home/merganser/jinjin");
rm(list=ls());
dyn.load("dpm.so");
source("dpm.R");
library(MASS);

setwd("/home/merganser/jinjin/result")
load("newdata2.RData")
O1=out1$O1;
O0=out1$O0;
D0=O0-O1;
p=74;
n=180;
n1=n*25;
n0=n*30;

X1.t<- mvrnorm(n1,rep(0,p),solve(O1));
X0.t<- mvrnorm(n0,rep(0,p),solve(O0));

setwd("/home/merganser/jinjin")
fit.aic <- dpm(X1.t,X0.t,nlambda=10,tuning="aic");
k=fit.aic[[4]][[6]];
est = fit.aic[[1]][[4]];

write.matrix(format(est, scientific=FALSE), 
               file = paste("/home/merganser/jinjin/Jindpm/comparecode/result", "initialdpm.csv", sep="/"), sep=",")
               
========================================================================================
####### 2. dpmcvx.R #######
## based on Dave Zhao's "dpm.R" package
========================================================================================
rm(list=ls());
dyn.load("dpm.so");
source("dpm.R");
library(MASS);
library(ROCR);
p=74;
n=180;
n1=n*25;
n0=n*30;

dpmcvx<- function(X1,X0,
                lambda=NULL,nlambda=10,lambda.min.ratio=NULL,
                rho=NULL,shrink=NULL,prec=0.001,max.ite=100,
                correlation=FALSE,perturb=FALSE,
                tuning=c("none","aic","bic","cv"),folds=5)
{
lambda=NULL;
nlambda=10;
lambda.min.ratio=NULL;
rho=NULL;
shrink=NULL;
prec=0.001;
max.ite=100;
correlation=FALSE;
perturb=FALSE;
tuning=c("none","aic","bic","cv");
folds=5;
    ## ==========================================================
    ## calculate covariance matrices, calculate lambdas, etc
    ## ==========================================================
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
    
    ## determine lambda based on S and e
    if(!is.null(lambda)){ nlambda <- length(lambda); }
    if(is.null(lambda))
        {
            if(is.null(lambda.min.ratio)){ lambda.min.ratio <- 0.04; }
            ##lambda.max <- min(max(S-diag(diag(S))),
            ##                -min(S-diag(diag(S))));
            lambda.max <- max(abs(e));
            lambda.min <- lambda.min.ratio*lambda.max;
            lambda <- exp(seq(log(lambda.max),log(lambda.min),
                              length=nlambda));
        }
    gamma <- 1/rho;
    lambda <- lambda - shrink*prec;
    nlambda = length(lambda);
    ## ===============================================================================
    ## Calling Matlab and use CVX package to solve constrained optimization problem
    ## ===============================================================================
  write.matrix(format(S, scientific=FALSE), 
               file = paste("/home/merganser/jinjin/Jindpm/CVX", "S.csv", sep="/"), sep=",")
  write.matrix(format(e, scientific=FALSE), 
               file = paste("/home/merganser/jinjin/Jindpm/CVX", "e.csv", sep="/"), sep=",")
  write.matrix(format(lambda, scientific=FALSE), 
               file = paste("/home/merganser/jinjin/Jindpm/CVX", "lambda.csv", sep="/"), sep=",")
setwd("/home/merganser/jinjin/Jindpm/CVX")

#matlab code reads in our variable x, creates two variables y and z, 
#and write z in a csv file
matlab.lines <- c(
"S=csvread('/home/merganser/jinjin/Jindpm/CVX/S.csv');",
"e=csvread('/home/merganser/jinjin/Jindpm/CVX/e.csv');",
"lambda=csvread('/home/merganser/jinjin/Jindpm/CVX/lambda.csv');",
"n1 = 4500; n2=5400; p = 74; d=p*(p+1)/2; nlambda=10;",
"beta=zeros(nlambda*d,1)",
"for i=1:nlambda",
"cvx_begin",
    "variable x(d)",
   " minimize( norm(x,1) )",
   " subject to",
        "norm( S*x-e, 1 ) <= lambda(i)",
"cvx_end",
"beta(((i-1)*d+1):(i*d),1)=x",
"end ",
"save /home/merganser/jinjin/Jindpm/CVX/beta1.mat beta",
)

#create a MATLAB script containing all the commands in matlab.lines
writeLines(matlab.lines, con="/home/merganser/jinjin/Jindpm/CVX/myscript.m")
 
#run MATLAB script
system("matlab -nodisplay -r \"run('/home/merganser/jinjin/Jindpm/CVX/myscript.m'); exit\"")

setwd("/home/merganser/jinjin/Jindpm/R.matlab_3.1.1/R.matlab")
library(R.matlab)
path <- system.file("mat-files", package="R.matlab")
pathname <- file.path("/home/merganser/jinjin/Jindpm/CVX", "beta1.mat")
data <- readMat(pathname)
diff=data$beta

ret = vector("list", nlambda);
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
    
    fit=(list(dpm=ret,lambda=lambda,nlambda=nlambda,opt=opt));
return(fit);
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
dpm.cv <- function(X1,X0,
                   lambda=NULL,
                   rho=NULL,shrink=NULL,prec=0.001,max.ite=100,
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
dpm.ic <- function(S1,S0,ret,n,penalty)
{
    lowtri <- which(lower.tri(ret[[1]],diag=TRUE));
    df <- sapply(ret,function(x){ sum(x[lowtri]!=0); });
    ic <- scale(n*sapply(ret,loss,S1,S0),
                center=-penalty*df,scale=FALSE);
    return(apply(ic,1,which.min));
}

setwd("/home/merganser/jinjin/result")
load("newdata2.RData")
O1=out1$O1;
O0=out1$O0;
D0=O0-O1;
p=74;
n=180;
n1=n*25;
n0=n*30;
X1<- mvrnorm(n1,rep(0,p),solve(O1));
X0<- mvrnorm(n0,rep(0,p),solve(O0));

setwd("/home/merganser/jinjin")
fitcvx.aic<- dpmcvx(X1,X0,nlambda=10,tuning="aic");
fitcvx.cv<- dpmcvx(X1,X0,nlambda=10,tuning="cv",folds=3);
k=fitcvx.aic[[4]][[5]];
est = fit.aic[[1]][[k]]                # To calculate estimation of Difference between two inverse covariance metrices #

estimation=c(1:p*p);
true=c(1:p*p);
k=1;
for (i in (1:(p)))
{for (j in (1:(p)))
{estimation[k]=abs(est[i,j]);    #####Modified by the user#####
true[k]=abbs(D0[i,j]);
k=k+1;
}
}

for (i in (1:(p*p)))
{if (abs(true[i])>1e-04){true[i]=1;}       ##### Threshold Modified by the user #####
if (abs(true[i])<=1e-04){true[i]=0;}
if (abs(estimation[i])>1e-04){estimation[i]=1;}
if (abs(estimation[i])<=1e-04){estimation[i]=0;}
}

pred <- prediction(estimation,true)
perf <- performance(pred,"tpr","fpr")
setwd("/home/merganser/jinjin/jindpm/comparecode/result")
pdf("dpm1.pdf")
plot(perf,colorize=TRUE)
dev.off()

output= NA
output= list(estimation=estimation; true=true)
save(output, file="dpm1.RData")

##### function dpm.cv and dpm.ic are the same as the original ones #####
##### Then we can compare estimation from dpmcvx.R with estimation from dpm.R (Need different tuning parameters to obtain same results)

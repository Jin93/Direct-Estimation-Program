##### SimulationSICE.R
====================================================================================
Description: Modified Estimation and analysis of Differential Networks of Data with 
             Theoretical Covariance Matrix obtained by SICE of FASD data
====================================================================================

====================================================================================
####### 1. SICEgenerate.R #######
## based on Chen's SRC program
====================================================================================
rm(list=ls());
dyn.load("dpm.so");
source("dpm.R");
library(MASS);
library(ROCR);

setwd("~/")
dat = read.csv("data1.csv",sep = ",")
grp = read.csv("group.csv",sep = ",")
sdata = rep(grp$Group_Binary, each = 180)
sdata1 = sdata[1:(27*180)]
##1:27 15'1' 12'2' 28:55 16'1' 12'2'
Y = as.matrix(dat[1:(27*180),-c(1,2)])
for(i in 1:27){
  temp = apply(Y[(1+(i-1)*180):(i*180),], 2, mean)
  Y[(1+(i-1)*180):(i*180),] = t(t(Y[(1+(i-1)*180):(i*180),])-temp)
}
library(MASS)
library(mvtnorm)
a = Sys.time()
setwd("~/src") # set working directory

# The C and R functions are in C and R directory respectively, and both are under src directory
# Use the makefile to compile the C files
# Double check the file paths in the makefile
source("R/MGGM.R")   # the R wrapper for C functions
dyn.load('lib/c_funcs.so')  # load the external C functions

fmvnorm = function(p, covmat_inv, Y, mu){ # Y n*p, Y[j,]; mu L*p, mu[l,]
#mvnormden = -0.5 * p * log(2*pi) + 0.5 * log(det(covmat_inv )) -0.5 * t(Y-mu) %*% covmat_inv %*% (Y-mu)
mvnormden = -p * log(2*pi) + log(det(covmat_inv )) - drop(t(Y-mu) %*% covmat_inv %*% (Y-mu))
return(0.5*mvnormden) # return the log likelihood
}

emcv1 = function(Y, mu, ncluster, covinv, z){ 
  # EM algorithm for cross-validation with known mu, covinv, pie_i
  # initialize clustering using kmeans
  n = nrow(Y)
  p = ncol(Y)
  cvlogliksum = 0 
  for(j in 1:n){
    l = which(z[,j]==1) # find which class it belongs to
      cvlogliksum = cvlogliksum +  fmvnorm(p, covinv[,((l-1)*p+1):(l*p)], Y[j,],mu[l,])
  }
  output = list(ploglik = cvlogliksum) 
}

n = nrow(Y)
p = ncol(Y)
ncluster=2
L=ncluster
Lambda1.vec <- log(p)*c(10,5, 1, .5, 0.1, 0.05,0.01, 0.005, 0.001, 0.0001, 1e-5)  #lasso penlaty
Lambda2.vec <- log(p)*c(10,5, 1, .5, 0.1, 0.05,0.01, 0.005, 0.001, 0.0001, 1e-5) 
#tau <- c(0.001,0.01,0.1,1) #thresholding parameter 
tau = 0.01
graph_complete = matrix(0,2,L*(L-1)/2)
for (l1 in 1:(L-1)){ #dimension of graph needs to be changed when varying ncluster = L
  graph_complete[,(L*(l1-1)-(l1-1)*l1/2+1):(L*l1-l1*(l1+1)/2)] = rbind(rep(l1,L-l1),(l1+1):L)
}
graph <- graph_complete - 1

mu = matrix(nrow = ncluster, ncol = p)
S_bar <- matrix(0,p,L*p) 
for (l in 1:L){
  S_bar[,((l-1)*p+1):(l*p)] = var(Y[sdata1==l,])
  mu[l,] = mean(Y[sdata1==l,])
}
z =matrix(0, ncol = n, nrow =L)
for(j in 1:n){
  for(l in 1:L){
    if(sdata1[j] == l) z[l,j] =1
  }
}
nn = as.vector(table(sdata1))
covmat_inv1 = list()
target = -Inf
posdef = T
  for(i in 1: length(Lambda1.vec)){
    for(j in 1:length(Lambda2.vec)){
      sol_path <- MGGM.path(S_bar, nn, Lambda1.vec[i], Lambda2.vec[j], graph, tau)
      out = sol_path$sol_convex[,,1,1]
      for(l in 1:L){
        covmat_inv1[[l]] = out[,((l-1)*p+1):(l*p)]
        if(det(covmat_inv1[[l]]) <0) posdef = F
      }
      #newtarget = obj_L0(p, L, covmat_inv1, S_bar, nn, graph, Lambda1.vec[i], Lambda2.vec[j], tau[k]) 
      #newtarget = loglikcalc (Y, pie, mu, covmat_inv1)
      pie = nn/n
      newtarget = emcv1(Y, mu, ncluster, out, z)$ploglik
      #cat("tau=", tau[k],"lambda_1=",Lambda1.vec[i], "lambda_2=", Lambda2.vec[j], "loglik=", newtarget,"\n")
      if(posdef ==F){
        cat("lambda_1=",Lambda1.vec[i], "lambda_2=", Lambda2.vec[j], "loglik=", newtarget,"\n")
        } else {
        if (newtarget > target){
          target = newtarget
          inv_sigma1 = covmat_inv1
          targeti = i
          targetj = j
        }
      }
    }
  }
target
Lambda1.vec[targeti]
Lambda2.vec[targetj]
targeti
targetj

out = NA
out = list(ploglik = target, lambda1 = Lambda1.vec[targeti], lambda2 = Lambda2.vec[targetj],
           inv_sigma= inv_sigma1)
Sys.time()-a
setwd("/home/merganser/jinjin/result")
save(out, file="Olddata2.RData")
out1 = NA
out1=list(O1=out$inv_sigma[[1]], O0= out$inv_sigma[[2]])
save(out1, file = "newdata2.RData")

====================================================================================
####### 2. DEgenerate.R #######
## based on Chen's SRC program
====================================================================================
rm(list=ls());
dyn.load("dpm.so");
source("dpmchangelambda.R");
library(MASS);
library(ROCR);

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
fit.aic <- dpm(X1.t,X0.t,nlambda=100,tuning="aic");
fit.cv <- dpm(X1.t,X0.t,nlambda=100,tuning="cv",folds=3);
setwd("/home/merganser/jinjin/TestDE/lambda/Result")
save(fit.aic,fit.cv, file="reallambda.RData")

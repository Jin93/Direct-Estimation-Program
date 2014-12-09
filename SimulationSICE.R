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

====================================================================================
####### 3. SimulationAnalysis.R #######
## ROC, Calculation of TPN etc.
====================================================================================
library(MASS);
library(ROCR);
p=74;
n=180;
n1=n*25;
n0=n*30;
nlambda=10;
setwd("/home/merganser/jinjin/result")
load("newdata2.RData")
O1=out1$O1;
O0=out1$O0;
D0=O0-O1;
setwd("/home/merganser/jinjin/TestDE/ROC")
load("fitgenerate.RData")

i=fit.aic[[4]][[5]];
j=fit.cv[[4]][[5]];
for (i in (1:10))           #####################absolute value!####################
{fit.aic[[1]][[i]]=abs(fit.aic[[1]][[i]])
fit.cv[[1]][[j]]=abs(fit.cv[[1]][[j]])
}
D0=abs(D0);

############AIC############
esaic=matrix(0,p,p*nlambda)
list=list()
k=c(1:6);
for (i in (1:6))
{k[i]=fit.aic[[4]][[i]];
esaic[,(p*(i-1)+1):(p*i)] = fit.aic[[1]][[k[i]]]
}

eaic=list();
eaic$a1=esaic[,(p*(1-1)+1):(p*1)];
eaic$a2=esaic[,(p*(2-1)+1):(p*2)];
eaic$a3=esaic[,(p*(3-1)+1):(p*3)];
eaic$a4=esaic[,(p*(4-1)+1):(p*4)];
eaic$a5=esaic[,(p*(5-1)+1):(p*5)];
eaic$a6=esaic[,(p*(6-1)+1):(p*6)];
############CV##############
escv=matrix(0,p,p*nlambda)
k=c(1:6);
for (i in (1:6))
{k[i]=fit.cv[[4]][[i]];
escv[,(p*(i-1)+1):(p*i)] = fit.cv[[1]][[k[i]]]}
ecv=list();
ecv$a1=escv[,(p*(1-1)+1):(p*1)];
ecv$a2=escv[,(p*(2-1)+1):(p*2)];
ecv$a3=escv[,(p*(3-1)+1):(p*3)];
ecv$a4=escv[,(p*(4-1)+1):(p*4)];
ecv$a5=escv[,(p*(5-1)+1):(p*5)];
ecv$a6=escv[,(p*(6-1)+1):(p*6)];

####################   ROC   ###################
setwd("/home/merganser/jinjin/TestDE/ROC")
library(ROCR);
######1######
q=5;  ###specified by the user###
estimation=c(1:p*p);
true=c(1:p*p);
k=1;
for (i in (1:(p)))
{for (j in (1:(p)))
{estimation[k]=eaic[[q]][i,j];
true[k]=D0[i,j];
k=k+1;
}
}
for (i in (1:(p*p)))
{if (abs(true[i])>1e-04){true[i]=1;}
if (abs(true[i])<=1e-04){true[i]=0;}
if (abs(estimation[i])>1e-04){estimation[i]=1;}
if (abs(estimation[i])<=1e-04){estimation[i]=0;}
}
pred <- prediction(estimation,true)
perf <- performance(pred,"tpr","fpr")
pdf("ROCAIC1.pdf")
plot(perf,colorize=TRUE)
dev.off()

############AIC############
tpup=0;tpdown=0;
tnup=0;tndown=0;
tdup=0;tddown=0;
tndup=0;tnddown=0;
th=1.00000e-03  ### Specified by the user ###
#######################Loss Function
q=5; ###specified by the user###
for (i in 1:(p))
{
for (j in 1:(p))
{
if ((eaic[[q]][i,j]>th)&&(D0[i,j]>th)){tpup=tpup+1;}
if (D0[i,j]>th) {tpdown=tpdown+1;}
if ((eaic[[q]][i,j]<=th)&&(D0[i,j]<=th)) {tnup=tnup+1;}
if (D0[i,j]<=th){tndown=tndown+1;}
if ((D0[i,j]>th)&&(eaic[[q]][i,j]>th)){tdup=tdup+1;}
if (eaic[[q]][i,j]>th){tddown=tddown+1;}
if ((eaic[[q]][i,j]<=th)&&(D0[i,j]<=th)) {tndup=tndup+1;}
if (eaic[[q]][i,j]<=th){tnddown=tnddown+1;}
}
}

tp=tpup/tpdown;
tn=tnup/tndown;
td=tdup/tddown;
tnd=tndup/tnddown;

###### Frobenius norm estimation accuracy:######
Fnorm=0;
for (i in 1:(p))
{
for (j in 1:(p))
{
Fnorm=Fnorm+(eaic[[q]][i,j]-D0[i,j])^2;
}}
Fnorm=sqrt(Fnorm);

setwd("/home/merganser/jinjin/TestDE/Simulation/Result")
save(tp,tn,td,tnd,Fnorm, file="truepositiverate_etc.RData")

######################   CV   #######################
q=fit.aic[[4]][[5]];  
tpup2=0;tpdown2=0;
tnup2=0;tndown2=0;
tdup2=0;tddown2=0;
tndup2=0;tnddown2=0;
th=1.000000e-05  ### Specified by the user ###
for (i in 1:(p))
{
for (j in 1:(p))
{
if ((ecv[[q]][i,j]>th)&&(D0[i,j]>th)) {tpup2=tpup2+1;}
if (D0[i,j]>th) {tpdown2=tpdown2+1;}
if ((ecv[[q]][i,j]<=th)&&(D0[i,j]<=th)) {tnup2=tnup2+1;}
if (D0[i,j]<=th){tndown2=tndown2+1;}
if ((D0[i,j]>th)&&(ecv[[q]][i,j]>th)){tdup2=tdup2+1;}
if (ecv[[q]][i,j]>th){tddown2=tddown2+1;}
if ((ecv[[q]][i,j]<=th)&&(D0[i,j]<=th)) {tndup2=tndup2+1;}
if (ecv[[q]][i,j]<=th){tnddown2=tnddown2+1;}
}
}

tp2=tpup2/tpdown2;
tn2=tnup2/tndown2;
td2=tdup2/tddown2;
tnd2=tndup2/tnddown2;

###### Frobenius norm estimation accuracy:######
Fnorm2=0;
for (i in 1:(p))
{
for (j in 1:(p))
{
Fnorm2=Fnorm2+(ecv[[q]][i,j]-D0[i,j])^2;
}}
Fnorm2=sqrt(Fnorm2);

setwd("/home/merganser/jinjin/TestDE/Simulation/Result")
save(tp2,tn2,td2,tnd2,Fnorm2,Fnorm2, file="truepositiverate2_etc.RData")



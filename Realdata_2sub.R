##### Realdata_2sub.R
====================================================================================
Description: Modified Estimation and analysis of Differential Networks of FASD data 
             by dividing each of the two groups into two subgroups
====================================================================================

====================================================================================
####### 1. dpmchangelanbdaloss.R #######
## Same as the one in RealdataDirect.R
====================================================================================

====================================================================================
####### 2. 2subgenerate.R #######
====================================================================================

rm(list=ls());
dyn.load("dpm.so");
source("dpmchangelambda.R");
library(MASS);
library(ROCR);
## **************************************************************
## Read Data
## **************************************************************
data=read.csv("/home/merganser/jinjin/data1.csv",header=T)
group=read.csv("/home/merganser/jinjin/group.csv",header=T)
data=as.matrix(data)
group=as.matrix(group)
ns=55;
nt=180;
p=ncol(data);
nlambda=10;

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
########## Or Standardize by subtracting each entries of FASD matrix by the total mean #########
#m=c(1:p);
#for (j in 1:nt)
#{for (k in 3:p)
#{m[k]=mean(data[,k]);
#ret[[i]][j,k]=ret[[i]][j,k]-m[k];
#}}}
###############

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
}}}}

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
}}}}

X11.t=matrix(0,(floor(numcontrol/2))*nt,p-2)
X12.t=matrix(0,(ceiling(numcontrol/2))*nt,p-2)
X01.t=matrix(0,(floor(numFASD/2))*nt,p-2)
X02.t=matrix(0,(ceiling(numFASD/2))*nt,p-2)

for (i in 1:(floor(numcontrol/2)))
{X11.t[(nt*(i-1)+1):(nt*i),]=X1.t[(nt*(i-1)+1):(nt*i),];
}
for (i in 1:(ceiling(numcontrol/2)))
{X12.t[(nt*(i-1)+1):(nt*i),]=X1.t[(nt*(floor(numcontrol/2)+i-1)+1):(nt*(floor(numcontrol/2)+i)),];
}
for (i in 1:(floor(numFASD/2)))
{X01.t[(nt*(i-1)+1):(nt*i),]=X0.t[(nt*(i-1)+1):(nt*i),];
}
for (i in 1:(ceiling(numFASD/2)))
{X02.t[(nt*(i-1)+1):(nt*i),]=X0.t[(nt*(floor(numFASD/2)+i-1)+1):(nt*(floor(numFASD/2)+i)),];
}
##### Calculation #####
fit.aic1 <- dpm(X11.t,X01.t,nlambda=10,tuning="aic");
fit.cv1 <- dpm(X11.t,X01.t,nlambda=10,tuning="cv",folds=3);
fit.aic2 <- dpm(X12.t,X02.t,nlambda=10,tuning="aic");
fit.cv2 <- dpm(X12.t,X02.t,nlambda=10,tuning="cv",folds=3);
setwd("/home/merganser/jinjin/TestDE/2sub/Result")
save(fit.aic1,fit.cv1,fit.aic2,fit.cv2, file="2sublambda2.RData")

====================================================================================
####### 3. 2subanalysis.R #######
====================================================================================

setwd("/home/merganser/jinjin/TestDE/2sub/Result")
load("2sublambda.RData")
p=74;
q=fit.aic1[[4]][[5]];###72###
r=fit.aic2[[4]][[5]];###70###
#########Decide the threshold (th) (Can be run outside this program)#########
#no1=0;
#for (i in 1:(p))
#{for (j in 1:(p))
#{if (s2[i,j]>2.000000e-04){no1=no1+1;}
#}}
#no1;
#########
th=5e-04; ###specify by user###
m=fit.aic1[[1]][[q]];
mm=fit.aic2[[1]][[r]];
x1=matrix(0,1,p*p);y1=matrix(0,1,p*p);
x2=matrix(0,1,p*p);y2=matrix(0,1,p*p);
x3=matrix(0,1,p*p);y3=matrix(0,1,p*p);
k1=1;k2=1;
for (i in (1:(p)))
{for (j in (1:(i)))
{ if (m[i,j]>=th) {x1[1,k1]=j;y1[1,k1]=i;k1=k1+1;}
  if (mm[i,j]>=th) {x2[1,k2]=j;y2[1,k2]=i;k2=k2+1;}
}
}
n1=0;n2=0;
for (i in 1:(p*p))
{
if (x1[1,i]>0) {n1=n1+1;}
if (x2[1,i]>0) {n2=n2+1;}
}
k3=1;
for (i in (1:(n1)))
{for (j in (1:(n2)))
{ if ((x1[1,i]==x2[1,j])&(y1[1,i]==y2[1,j])&(x1[1,i]!=0)&(x2[1,j]!=0)) {x3[1,k3]=x1[1,i];y3[1,k3]=y1[1,i];k3=k3+1;}
}
}
pdf("plotNonzeroDE2sub.pdf")
par(mfrow=c(2,2))
plot(x1,y1,type = "p",cex=.35,pch=21,col="red", xlim=c(0,75),ylim=c(0,75))
plot(x2,y2,type = "p",cex=.35,pch=21,col="green",xlim=c(0,75),ylim=c(0,75))
plot(x3,y3,type = "p",cex=.4,pch=21,col="black",bg="black",xlim=c(0,75),ylim=c(0,75))
plot(x1,y1,type = "p",cex=.35,pch=21,col="red",xlim=c(0,75),ylim=c(0,75))
par(new=TRUE)
plot(x2,y2,type = "p",cex=.35,pch=21,col="green",xlim=c(0,75),ylim=c(0,75))
par(new=TRUE)
plot(x3,y3,type = "p",cex=.4,pch=21,col="black",bg="black",xlim=c(0,75),ylim=c(0,75))
dev.off()

pdf("Combine.pdf")
plot(x1,y1,type = "p",cex=.35,pch=21,col="red",xlim=c(0,75),ylim=c(0,75))
par(new=TRUE)
plot(x2,y2,type = "p",cex=.35,pch=21,col="green",xlim=c(0,75),ylim=c(0,75))
par(new=TRUE)
plot(x3,y3,type = "p",cex=.4,pch=21,col="black",bg="black",xlim=c(0,75),ylim=c(0,75))
dev.off()

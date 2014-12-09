######Direct-Estimation-Program
=========================
################################################################
README, 12.08.2014
Program of Direct estimation of differential networks
Jin Jin
################################################################
==============================================================
######Description
==============================================================
R program for Direct Estimation, its modifications and applications on FASD data. 
R code for Direct Estimation was modified based on code from the dpm.R by 
Sihai D. Zhao, T. Tony Cai, and Hongzhe Li (https://github.com/sdzhao/dpm/blob/master/dpm.R).
  
==============================================================
#####Contents
==============================================================
######Realdata_Direct.R
- Description: Modified (tuning parameters, loss function, etc.) Estimation and analysis 
               of Differential Networks of FASD data by directly using FASD data matrix.
######Realdata_2sub.R
- Description: Modified Estimation and analysis of Differential Networks of FASD data 
               by dividing each of the two groups into two subgroups.
######SimulationSICE.R
- Description: Modified Estimation and analysis of Differential Networks of Data with 
               Theoretical Covariance Matrix obtained by SICE of FASD data.
######Dpmcvx.R
- Description: Modified Estimation and analysis of Differential Networks of Simulation 
               based on FASD Data Calling CVX package from Matlab in R. R code was modified 
               based on code from the dpm.R (https://github.com/sdzhao/dpm/blob/master/dpm.R).
               Use SICE as true precision matrices to generate data and compare dpmcvx.R
               and dpm.R.
All of the programs need compilation when used.

### Direct-Estimation-Program
=========================
################################################################
### README, 12.08.2014
### Program of Direct estimation of differential networks
### Jin Jin
################################################################
==============================================================
Description
==============================================================
R program for Direct Estimation, its modifications and applications on FASD data. R code for Direct Estimation was modified based on code from the dpm.R by Sihai D. Zhao, T. Tony Cai, and Hongzhe Li (https://github.com/sdzhao/dpm/blob/master/dpm.R).

==============================================================
Contents
==============================================================
dpm.c
- C code
- Linux: compile it using R CMD SHLIB dpm.c
- Windows: see http://mcglinn.web.unc.edu/blog/linking-c-with-r-in-windows/

dpm.R
- R wrapper function for dpm.c

example.R
- example usage



## ##############################################################
## README, 5.24.2014
## Direct estimation of differential networks
## Sihai D. Zhao, T. Tony Cai, and Hongzhe Li
## ##############################################################
## ==============================================================
## Description
## ==============================================================
C implementation of methods, along with R wrapper function and example usage. C code was modified based on code from the flare R package (http://cran.r-project.org/web/packages/flare/index.html).

## ==============================================================
## Contents
## ==============================================================
dpm.c
- C code
- Linux: compile it using R CMD SHLIB dpm.c
- Windows: see http://mcglinn.web.unc.edu/blog/linking-c-with-r-in-windows/

dpm.R
- R wrapper function for dpm.c

example.R
- example usage

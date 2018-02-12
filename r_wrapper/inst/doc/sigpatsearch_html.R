## ----include = FALSE-----------------------------------------------------
library(sigpatsearch)

## ------------------------------------------------------------------------
library(sigpatsearch)

## ------------------------------------------------------------------------
sig_fais1 <- sigpatsearch(method='fais')

## ------------------------------------------------------------------------
print(sig_fais1)

## ------------------------------------------------------------------------
sig_fais2 <- sigpatsearch(method='FAIS')
sig_fais3 <- sigpatsearch(method='FaIs')
sig_fais4 <- sigpatsearch(use_intervals=T) 
sig_fais5 <- sigpatsearch(use_intervals=T, use_covariate=F) #default use_covariate is False
sig_fais6 <- sigpatsearch(use_combinations=F)

## ------------------------------------------------------------------------
sig_fastcmh1 <- sigpatsearch(method='fastcmh')
sig_fastcmh2 <- sigpatsearch(method='FASTCMH')
sig_fastcmh3 <- sigpatsearch(method='FastCMH')
sig_fastcmh4 <- sigpatsearch(use_intervals=T, use_covariate=T)
sig_fastcmh5 <- sigpatsearch(use_combinations=F, use_covariate=T)
print(sig_fastcmh5)

## ------------------------------------------------------------------------
sig_facs1 <- sigpatsearch(method='facs')
sig_facs2 <- sigpatsearch(method='FACS')
sig_facs3 <- sigpatsearch(method='FaCs')
sig_facs4 <- sigpatsearch(use_intervals=F, use_covariate=T)
sig_facs5 <- sigpatsearch(use_intervals=F, use_covariate=F)
sig_facs6 <- sigpatsearch(use_combinations=T, use_covariate=T)
sig_facs7 <- sigpatsearch(use_combinations=T, use_covariate=F)
print(sig_facs7)

## ------------------------------------------------------------------------
sig_fais8 <- sigpatsearch(method='fais')
sig_fais8$set_alpha(0.001)
sig_fais8$set_lmax(10)
print(sig_fais8)


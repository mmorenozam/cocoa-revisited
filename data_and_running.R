#######################################################################
#                                                                     #
# R Script for computation of posterior predictions of the manuscript #
# "Exploring cocoa bean fermentation processes by kinetic modelling"  #
#                                                                     #
# The following script contains data from 28 fermentation trials. For #
# a more detailed description of datasets, see Supplementary Material #
# Section 1, and Supplementary Table S1.                              #
#                                                                     #
# This code works in prompt shell of Linux OS. For its execution it   #
# needs seven arguments:                                              #
# d = dataset (integer from 1 to 28):                                 #
#     1:  ghhp1           15: ghhp7                                   #
#     2:  dowb1           16: brwb3                                   #
#     3:  cipb1           17: brwb4                                   #
#     4:  brwb1           18: brst2                                   #
#     5:  brwb2           19: brwb5                                   #
#     6:  ecpt1           20: brwb6                                   #
#     7:  ecpt2           21: brwb7                                   #
#     8:  ecwb1           22: mywb3                                   #
#     9:  ecwb2           23: brwb8                                   #
#     10: brpb1           24: hnwb1                                   #
#     11: brst1           25: ecpb1                                   #
#     12: mywb1           26: niwb1                                   #
#     13: mywb2           27: niwb2                                   #
#     14: ghhp6           28: ecpb2                                   #
#                                                                     #
# mi = model iteration (integer), refers to model iteration code, for #
#      instance, MI(1,2,3), here is called as 123                     #
#                                                                     #
# The following are Stan parameters. For further understanding, see   #
# reference manual for Stan. Between parenthesis are noted the values #
# used in our paper                                                   #
# ad = adapt_delta parameter (0.995 to 0.999)                         #
#                                                                     #
# st = step_size parameter (0.1)                                      #
#                                                                     #
# mt = max_treedepth parameter (10)                                   #
#                                                                     #
# s  = seed parameter (141085)                                        #
#                                                                     #
# ni = number of iterations (3000)                                    #
#                                                                     #
# Hence, for running this script you must follow the following form:  #
#                                                                     #
# $ Rscript post_sims1.R d mi ad st mt s ni                           #
#                                                                     #
# Once the script is done running, an .Rsave file is saved with the   #
# dataset coding as name containing the Stan output.                  #
#                                                                     #
# Code written by Mauricio Moreno-Zambrano                            #
# 2021                                                                #
#                                                                     #
#######################################################################


setwd('/home') # directory name, change accordingly if needed
rm(list = ls())
args = commandArgs(TRUE)

library(rstan)

rstan_options(auto_write = TRUE)

source("trial_data.R")

# datasets ordered as mentioned above, now defined in trial_data.R

trial <- get_trial(as.numeric(args[1]))
T <- trial$T
x0 <- trial$x0
t0 <- trial$t0
ts <- trial$ts
x <- trial$x
lbl <- trial$lbl
x1 <- trial$x1
scl <- trial$scl


#depending on chosen MI, initial values are randomly generated for MCMC-NUTS sampler

if (args[2]==0){
  ini = function(){
    list(mu1=abs(rnorm(1,0.5,0.2)),mu2=abs(rnorm(1,0.5,0.2)),
         mu3=abs(rnorm(1,0.5,0.2)),mu4=abs(rnorm(1,0.5,0.2)),
         mu5=abs(rnorm(1,0.5,0.2)),ks1=abs(rnorm(1,0.5,0.2)),
         ks2=abs(rnorm(1,0.5,0.2)),ks3=abs(rnorm(1,0.5,0.2)),
         ks4=abs(rnorm(1,0.5,0.2)),ks5=abs(rnorm(1,0.5,0.2)),
         k1=abs(rnorm(1,0.5,0.2)),k2=abs(rnorm(1,0.5,0.2)),
         k3=abs(rnorm(1,0.5,0.2)),yc1=abs(rnorm(1,0.5,0.2)),
         yc2=abs(rnorm(1,0.5,0.2)),yc3=abs(rnorm(1,0.5,0.2)),
         yc4=abs(rnorm(1,0.5,0.2)),yc5=abs(rnorm(1,0.5,0.2)),
         yc6=abs(rnorm(1,0.5,0.2)),yc7=abs(rnorm(1,0.5,0.2)),
         yc8=abs(rnorm(1,0.5,0.2)),yc9=abs(rnorm(1,0.5,0.2)),
         yc10=abs(rnorm(1,0.5,0.2)),yc11=abs(rnorm(1,0.5,0.2)),
         sigma=abs(rnorm(1,0.5,0.2)))
  }
}
if (args[2]==1){
  ini = function(){
    list(mu1=abs(rnorm(1,0.5,0.2)),mu2=abs(rnorm(1,0.5,0.2)),
         mu3=abs(rnorm(1,0.5,0.2)),mu4=abs(rnorm(1,0.5,0.2)),
         mu5=abs(rnorm(1,0.5,0.2)),ks1=abs(rnorm(1,0.5,0.2)),
         ks2=abs(rnorm(1,0.5,0.2)),ks3=abs(rnorm(1,0.5,0.2)),
         ks4=abs(rnorm(1,0.5,0.2)),ks5=abs(rnorm(1,0.5,0.2)),
         k1=abs(rnorm(1,0.5,0.2)),k2=abs(rnorm(1,0.5,0.2)),
         k3=abs(rnorm(1,0.5,0.2)),yc1=abs(rnorm(1,0.5,0.2)),
         yc2=abs(rnorm(1,0.5,0.2)),yc3=abs(rnorm(1,0.5,0.2)),
         yc4=abs(rnorm(1,0.5,0.2)),yc5=abs(rnorm(1,0.5,0.2)),
         yc6=abs(rnorm(1,0.5,0.2)),yc7=abs(rnorm(1,0.5,0.2)),
         yc8=abs(rnorm(1,0.5,0.2)),yc9=abs(rnorm(1,0.5,0.2)),
         yc10=abs(rnorm(1,0.5,0.2)),yc11=abs(rnorm(1,0.5,0.2)),
         ev1=abs(rnorm(1,0.5,0.2)),ev2=abs(rnorm(1,0.5,0.2)),
         ev3=abs(rnorm(1,0.5,0.2)),sigma=abs(rnorm(1,0.5,0.2)))
  }
}
if (args[2]==2){
  ini = function(){
    list(mu1=abs(rnorm(1,0.5,0.2)),mu2=abs(rnorm(1,0.5,0.2)),
         mu3=abs(rnorm(1,0.5,0.2)),mu4=abs(rnorm(1,0.5,0.2)),
         mu5=abs(rnorm(1,0.5,0.2)),mu6=abs(rnorm(1,0.5,0.2)),
         ks1=abs(rnorm(1,0.5,0.2)),ks2=abs(rnorm(1,0.5,0.2)),
         ks3=abs(rnorm(1,0.5,0.2)),ks4=abs(rnorm(1,0.5,0.2)),
         ks5=abs(rnorm(1,0.5,0.2)),ks6=abs(rnorm(1,0.5,0.2)),
         k1=abs(rnorm(1,0.5,0.2)),k2=abs(rnorm(1,0.5,0.2)),
         k3=abs(rnorm(1,0.5,0.2)),yc1=abs(rnorm(1,0.5,0.2)),
         yc2=abs(rnorm(1,0.5,0.2)),yc3=abs(rnorm(1,0.5,0.2)),
         yc4=abs(rnorm(1,0.5,0.2)),yc5=abs(rnorm(1,0.5,0.2)),
         yc6=abs(rnorm(1,0.5,0.2)),yc7=abs(rnorm(1,0.5,0.2)),
         yc8=abs(rnorm(1,0.5,0.2)),yc9=abs(rnorm(1,0.5,0.2)),
         yc10=abs(rnorm(1,0.5,0.2)),yc11=abs(rnorm(1,0.5,0.2)),
         yc12=abs(rnorm(1,0.5,0.2)),yc13=abs(rnorm(1,0.5,0.2)),
         yc14=abs(rnorm(1,0.5,0.2)),yc15=abs(rnorm(1,0.5,0.2)),
         yc16=abs(rnorm(1,0.5,0.2)),sigma=abs(rnorm(1,0.5,0.2)))
  }
}
if (args[2]==3){
  ini = function(){
    list(mu1=abs(rnorm(1,0.5,0.2)),mu2=abs(rnorm(1,0.5,0.2)),
         mu3=abs(rnorm(1,0.5,0.2)),mu4=abs(rnorm(1,0.5,0.2)),
         mu5=abs(rnorm(1,0.5,0.2)),ks1=abs(rnorm(1,0.5,0.2)),
         ks2=abs(rnorm(1,0.5,0.2)),ks3=abs(rnorm(1,0.5,0.2)),
         ks4=abs(rnorm(1,0.5,0.2)),ks5=abs(rnorm(1,0.5,0.2)),
         k1=abs(rnorm(1,0.5,0.2)),k2=abs(rnorm(1,0.5,0.2)),
         k3=abs(rnorm(1,0.5,0.2)),yc1=abs(rnorm(1,0.5,0.2)),
         yc2=abs(rnorm(1,0.5,0.2)),yc3=abs(rnorm(1,0.5,0.2)),
         yc4=abs(rnorm(1,0.5,0.2)),yc5=abs(rnorm(1,0.5,0.2)),
         yc6=abs(rnorm(1,0.5,0.2)),yc7=abs(rnorm(1,0.5,0.2)),
         yc8=abs(rnorm(1,0.5,0.2)),yc9=abs(rnorm(1,0.5,0.2)),
         yc10=abs(rnorm(1,0.5,0.2)),yc11=abs(rnorm(1,0.5,0.2)),
         yc12=abs(rnorm(1,0.5,0.2)),yc13=abs(rnorm(1,0.5,0.2)),
         sigma=abs(rnorm(1,0.5,0.2)))
  }
}

if (args[2]==4){
  ini = function(){
    list(mu1=abs(rnorm(1,0.5,0.2)),mu2=abs(rnorm(1,0.5,0.2)),
         mu3=abs(rnorm(1,0.5,0.2)),mu4=abs(rnorm(1,0.5,0.2)),
         mu5=abs(rnorm(1,0.5,0.2)),mu6=abs(rnorm(1,0.5,0.2)),
         ks1=abs(rnorm(1,0.5,0.2)),ks2=abs(rnorm(1,0.5,0.2)),
         ks3=abs(rnorm(1,0.5,0.2)),ks4=abs(rnorm(1,0.5,0.2)),
         ks5=abs(rnorm(1,0.5,0.2)),ks6=abs(rnorm(1,0.5,0.2)),
         k1=abs(rnorm(1,0.5,0.2)),k2=abs(rnorm(1,0.5,0.2)),
         k3=abs(rnorm(1,0.5,0.2)),yc1=abs(rnorm(1,0.5,0.2)),
         yc2=abs(rnorm(1,0.5,0.2)),yc3=abs(rnorm(1,0.5,0.2)),
         yc4=abs(rnorm(1,0.5,0.2)),yc5=abs(rnorm(1,0.5,0.2)),
         yc6=abs(rnorm(1,0.5,0.2)),yc7=abs(rnorm(1,0.5,0.2)),
         yc8=abs(rnorm(1,0.5,0.2)),yc9=abs(rnorm(1,0.5,0.2)),
         yc10=abs(rnorm(1,0.5,0.2)),yc11=abs(rnorm(1,0.5,0.2)),
         yc12=abs(rnorm(1,0.5,0.2)),yc13=abs(rnorm(1,0.5,0.2)),
         sigma=abs(rnorm(1,0.5,0.2)))
  }
}

if (args[2]==5){
  ini = function(){
    list(mu1=abs(rnorm(1,0.5,0.2)),mu2=abs(rnorm(1,0.5,0.2)),
         mu3=abs(rnorm(1,0.5,0.2)),mu4=abs(rnorm(1,0.5,0.2)),
         mu5=abs(rnorm(1,0.5,0.2)),mu6=abs(rnorm(1,0.5,0.2)),
         ks1=abs(rnorm(1,0.5,0.2)),ks2=abs(rnorm(1,0.5,0.2)),
         ks3=abs(rnorm(1,0.5,0.2)),ks4=abs(rnorm(1,0.5,0.2)),
         ks5=abs(rnorm(1,0.5,0.2)),ks6=abs(rnorm(1,0.5,0.2)),
         k1=abs(rnorm(1,0.5,0.2)),k2=abs(rnorm(1,0.5,0.2)),
         k3=abs(rnorm(1,0.5,0.2)),yc1=abs(rnorm(1,0.5,0.2)),
         yc2=abs(rnorm(1,0.5,0.2)),yc3=abs(rnorm(1,0.5,0.2)),
         yc4=abs(rnorm(1,0.5,0.2)),yc5=abs(rnorm(1,0.5,0.2)),
         yc6=abs(rnorm(1,0.5,0.2)),yc7=abs(rnorm(1,0.5,0.2)),
         yc8=abs(rnorm(1,0.5,0.2)),yc9=abs(rnorm(1,0.5,0.2)),
         yc10=abs(rnorm(1,0.5,0.2)),yc11=abs(rnorm(1,0.5,0.2)),
         yc12=abs(rnorm(1,0.5,0.2)),sigma=abs(rnorm(1,0.5,0.2)))
  }
}
if (args[2]==12){
  ini = function(){
    list(mu1=abs(rnorm(1,0.5,0.2)),mu2=abs(rnorm(1,0.5,0.2)),
         mu3=abs(rnorm(1,0.5,0.2)),mu4=abs(rnorm(1,0.5,0.2)),
         mu5=abs(rnorm(1,0.5,0.2)),mu6=abs(rnorm(1,0.5,0.2)),
         ks1=abs(rnorm(1,0.5,0.2)),ks2=abs(rnorm(1,0.5,0.2)),
         ks3=abs(rnorm(1,0.5,0.2)),ks4=abs(rnorm(1,0.5,0.2)),
         ks5=abs(rnorm(1,0.5,0.2)),ks6=abs(rnorm(1,0.5,0.2)),
         k1=abs(rnorm(1,0.5,0.2)),k2=abs(rnorm(1,0.5,0.2)),
         k3=abs(rnorm(1,0.5,0.2)),yc1=abs(rnorm(1,0.5,0.2)),
         yc2=abs(rnorm(1,0.5,0.2)),yc3=abs(rnorm(1,0.5,0.2)),
         yc4=abs(rnorm(1,0.5,0.2)),yc5=abs(rnorm(1,0.5,0.2)),
         yc6=abs(rnorm(1,0.5,0.2)),yc7=abs(rnorm(1,0.5,0.2)),
         yc8=abs(rnorm(1,0.5,0.2)),yc9=abs(rnorm(1,0.5,0.2)),
         yc10=abs(rnorm(1,0.5,0.2)),yc11=abs(rnorm(1,0.5,0.2)),
         yc12=abs(rnorm(1,0.5,0.2)),yc13=abs(rnorm(1,0.5,0.2)),
         yc14=abs(rnorm(1,0.5,0.2)),yc15=abs(rnorm(1,0.5,0.2)),
         yc16=abs(rnorm(1,0.5,0.2)),sigma=abs(rnorm(1,0.5,0.2)),
         ev1=abs(rnorm(1,0.5,0.2)),ev2=abs(rnorm(1,0.5,0.2)),
         ev3=abs(rnorm(1,0.5,0.2)))
  }
}
if (args[2]==13){
  ini = function(){
    list(mu1=abs(rnorm(1,0.5,0.2)),mu2=abs(rnorm(1,0.5,0.2)),
         mu3=abs(rnorm(1,0.5,0.2)),mu4=abs(rnorm(1,0.5,0.2)),
         mu5=abs(rnorm(1,0.5,0.2)),ks1=abs(rnorm(1,0.5,0.2)),
         ks2=abs(rnorm(1,0.5,0.2)),ks3=abs(rnorm(1,0.5,0.2)),
         ks4=abs(rnorm(1,0.5,0.2)),ks5=abs(rnorm(1,0.5,0.2)),
         k1=abs(rnorm(1,0.5,0.2)),k2=abs(rnorm(1,0.5,0.2)),
         k3=abs(rnorm(1,0.5,0.2)),yc1=abs(rnorm(1,0.5,0.2)),
         yc2=abs(rnorm(1,0.5,0.2)),yc3=abs(rnorm(1,0.5,0.2)),
         yc4=abs(rnorm(1,0.5,0.2)),yc5=abs(rnorm(1,0.5,0.2)),
         yc6=abs(rnorm(1,0.5,0.2)),yc7=abs(rnorm(1,0.5,0.2)),
         yc8=abs(rnorm(1,0.5,0.2)),yc9=abs(rnorm(1,0.5,0.2)),
         yc10=abs(rnorm(1,0.5,0.2)),yc11=abs(rnorm(1,0.5,0.2)),
         yc12=abs(rnorm(1,0.5,0.2)),yc13=abs(rnorm(1,0.5,0.2)),
         ev1=abs(rnorm(1,0.5,0.2)),ev2=abs(rnorm(1,0.5,0.2)),
         ev3=abs(rnorm(1,0.5,0.2)),sigma=abs(rnorm(1,0.5,0.2)))
  }
}
if (args[2]==14){
  ini = function(){
    list(mu1=abs(rnorm(1,0.5,0.2)),mu2=abs(rnorm(1,0.5,0.2)),
         mu3=abs(rnorm(1,0.5,0.2)),mu4=abs(rnorm(1,0.5,0.2)),
         mu5=abs(rnorm(1,0.5,0.2)),mu6=abs(rnorm(1,0.5,0.2)),
         ks1=abs(rnorm(1,0.5,0.2)),ks2=abs(rnorm(1,0.5,0.2)),
         ks3=abs(rnorm(1,0.5,0.2)),ks4=abs(rnorm(1,0.5,0.2)),
         ks5=abs(rnorm(1,0.5,0.2)),ks6=abs(rnorm(1,0.5,0.2)),
         k1=abs(rnorm(1,0.5,0.2)),k2=abs(rnorm(1,0.5,0.2)),
         k3=abs(rnorm(1,0.5,0.2)),yc1=abs(rnorm(1,0.5,0.2)),
         yc2=abs(rnorm(1,0.5,0.2)),yc3=abs(rnorm(1,0.5,0.2)),
         yc4=abs(rnorm(1,0.5,0.2)),yc5=abs(rnorm(1,0.5,0.2)),
         yc6=abs(rnorm(1,0.5,0.2)),yc7=abs(rnorm(1,0.5,0.2)),
         yc8=abs(rnorm(1,0.5,0.2)),yc9=abs(rnorm(1,0.5,0.2)),
         yc10=abs(rnorm(1,0.5,0.2)),yc11=abs(rnorm(1,0.5,0.2)),
         yc12=abs(rnorm(1,0.5,0.2)),yc13=abs(rnorm(1,0.5,0.2)),
         ev1=abs(rnorm(1,0.5,0.2)),ev2=abs(rnorm(1,0.5,0.2)),
         ev3=abs(rnorm(1,0.5,0.2)),sigma=abs(rnorm(1,0.5,0.2)))
  }
}
if (args[2]==15){
  ini = function(){
    list(mu1=abs(rnorm(1,0.5,0.2)),mu2=abs(rnorm(1,0.5,0.2)),
         mu3=abs(rnorm(1,0.5,0.2)),mu4=abs(rnorm(1,0.5,0.2)),
         mu5=abs(rnorm(1,0.5,0.2)),mu6=abs(rnorm(1,0.5,0.2)),
         ks1=abs(rnorm(1,0.5,0.2)),ks2=abs(rnorm(1,0.5,0.2)),
         ks3=abs(rnorm(1,0.5,0.2)),ks4=abs(rnorm(1,0.5,0.2)),
         ks5=abs(rnorm(1,0.5,0.2)),ks6=abs(rnorm(1,0.5,0.2)),
         k1=abs(rnorm(1,0.5,0.2)),k2=abs(rnorm(1,0.5,0.2)),
         k3=abs(rnorm(1,0.5,0.2)),yc1=abs(rnorm(1,0.5,0.2)),
         yc2=abs(rnorm(1,0.5,0.2)),yc3=abs(rnorm(1,0.5,0.2)),
         yc4=abs(rnorm(1,0.5,0.2)),yc5=abs(rnorm(1,0.5,0.2)),
         yc6=abs(rnorm(1,0.5,0.2)),yc7=abs(rnorm(1,0.5,0.2)),
         yc8=abs(rnorm(1,0.5,0.2)),yc9=abs(rnorm(1,0.5,0.2)),
         yc10=abs(rnorm(1,0.5,0.2)),yc11=abs(rnorm(1,0.5,0.2)),
         yc12=abs(rnorm(1,0.5,0.2)),ev1=abs(rnorm(1,0.5,0.2)),
         ev2=abs(rnorm(1,0.5,0.2)),ev3=abs(rnorm(1,0.5,0.2)),
         sigma=abs(rnorm(1,0.5,0.2)))
  }
}
if (args[2]==23){
  ini = function(){
    list(mu1=abs(rnorm(1,0.5,0.2)),mu2=abs(rnorm(1,0.5,0.2)),
         mu3=abs(rnorm(1,0.5,0.2)),mu4=abs(rnorm(1,0.5,0.2)),
         mu5=abs(rnorm(1,0.5,0.2)),mu6=abs(rnorm(1,0.5,0.2)),
         ks1=abs(rnorm(1,0.5,0.2)),ks2=abs(rnorm(1,0.5,0.2)),
         ks3=abs(rnorm(1,0.5,0.2)),ks4=abs(rnorm(1,0.5,0.2)),
         ks5=abs(rnorm(1,0.5,0.2)),ks6=abs(rnorm(1,0.5,0.2)),
         k1=abs(rnorm(1,0.5,0.2)),k2=abs(rnorm(1,0.5,0.2)),
         k3=abs(rnorm(1,0.5,0.2)),yc1=abs(rnorm(1,0.5,0.2)),
         yc2=abs(rnorm(1,0.5,0.2)),yc3=abs(rnorm(1,0.5,0.2)),
         yc4=abs(rnorm(1,0.5,0.2)),yc5=abs(rnorm(1,0.5,0.2)),
         yc6=abs(rnorm(1,0.5,0.2)),yc7=abs(rnorm(1,0.5,0.2)),
         yc8=abs(rnorm(1,0.5,0.2)),yc9=abs(rnorm(1,0.5,0.2)),
         yc10=abs(rnorm(1,0.5,0.2)),yc11=abs(rnorm(1,0.5,0.2)),
         yc12=abs(rnorm(1,0.5,0.2)),yc13=abs(rnorm(1,0.5,0.2)),
         yc14=abs(rnorm(1,0.5,0.2)),yc15=abs(rnorm(1,0.5,0.2)),
         yc16=abs(rnorm(1,0.5,0.2)),yc17=abs(rnorm(1,0.5,0.2)),
         yc18=abs(rnorm(1,0.5,0.2)),sigma=abs(rnorm(1,0.5,0.2)))
  }
}
if (args[2]==24){
  ini = function(){
    list(mu1=abs(rnorm(1,0.5,0.2)),mu2=abs(rnorm(1,0.5,0.2)),
         mu3=abs(rnorm(1,0.5,0.2)),mu4=abs(rnorm(1,0.5,0.2)),
         mu5=abs(rnorm(1,0.5,0.2)),mu6=abs(rnorm(1,0.5,0.2)),
         mu7=abs(rnorm(1,0.5,0.2)),ks1=abs(rnorm(1,0.5,0.2)),
         ks2=abs(rnorm(1,0.5,0.2)),ks3=abs(rnorm(1,0.5,0.2)),
         ks4=abs(rnorm(1,0.5,0.2)),ks5=abs(rnorm(1,0.5,0.2)),
         ks6=abs(rnorm(1,0.5,0.2)),ks7=abs(rnorm(1,0.5,0.2)),
         k1=abs(rnorm(1,0.5,0.2)),k2=abs(rnorm(1,0.5,0.2)),
         k3=abs(rnorm(1,0.5,0.2)),yc1=abs(rnorm(1,0.5,0.2)),
         yc2=abs(rnorm(1,0.5,0.2)),yc3=abs(rnorm(1,0.5,0.2)),
         yc4=abs(rnorm(1,0.5,0.2)),yc5=abs(rnorm(1,0.5,0.2)),
         yc6=abs(rnorm(1,0.5,0.2)),yc7=abs(rnorm(1,0.5,0.2)),
         yc8=abs(rnorm(1,0.5,0.2)),yc9=abs(rnorm(1,0.5,0.2)),
         yc10=abs(rnorm(1,0.5,0.2)),yc11=abs(rnorm(1,0.5,0.2)),
         yc12=abs(rnorm(1,0.5,0.2)),yc13=abs(rnorm(1,0.5,0.2)),
         yc14=abs(rnorm(1,0.5,0.2)),yc15=abs(rnorm(1,0.5,0.2)),
         yc16=abs(rnorm(1,0.5,0.2)),yc17=abs(rnorm(1,0.5,0.2)),
         yc18=abs(rnorm(1,0.5,0.2)),sigma=abs(rnorm(1,0.5,0.2)))
  }
}
if (args[2]==25){
  ini = function(){
    list(mu1=abs(rnorm(1,0.5,0.2)),mu2=abs(rnorm(1,0.5,0.2)),
         mu3=abs(rnorm(1,0.5,0.2)),mu4=abs(rnorm(1,0.5,0.2)),
         mu5=abs(rnorm(1,0.5,0.2)),mu6=abs(rnorm(1,0.5,0.2)),
         mu7=abs(rnorm(1,0.5,0.2)),ks1=abs(rnorm(1,0.5,0.2)),
         ks2=abs(rnorm(1,0.5,0.2)),ks3=abs(rnorm(1,0.5,0.2)),
         ks4=abs(rnorm(1,0.5,0.2)),ks5=abs(rnorm(1,0.5,0.2)),
         ks6=abs(rnorm(1,0.5,0.2)),ks7=abs(rnorm(1,0.5,0.2)),
         k1=abs(rnorm(1,0.5,0.2)),k2=abs(rnorm(1,0.5,0.2)),
         k3=abs(rnorm(1,0.5,0.2)),yc1=abs(rnorm(1,0.5,0.2)),
         yc2=abs(rnorm(1,0.5,0.2)),yc3=abs(rnorm(1,0.5,0.2)),
         yc4=abs(rnorm(1,0.5,0.2)),yc5=abs(rnorm(1,0.5,0.2)),
         yc6=abs(rnorm(1,0.5,0.2)),yc7=abs(rnorm(1,0.5,0.2)),
         yc8=abs(rnorm(1,0.5,0.2)),yc9=abs(rnorm(1,0.5,0.2)),
         yc10=abs(rnorm(1,0.5,0.2)),yc11=abs(rnorm(1,0.5,0.2)),
         yc12=abs(rnorm(1,0.5,0.2)),yc13=abs(rnorm(1,0.5,0.2)),
         yc14=abs(rnorm(1,0.5,0.2)),yc15=abs(rnorm(1,0.5,0.2)),
         yc16=abs(rnorm(1,0.5,0.2)),yc17=abs(rnorm(1,0.5,0.2)),
         sigma=abs(rnorm(1,0.5,0.2)))
  }
}
if (args[2]==34){
  ini = function(){
    list(mu1=abs(rnorm(1,0.5,0.2)),mu2=abs(rnorm(1,0.5,0.2)),
         mu3=abs(rnorm(1,0.5,0.2)),mu4=abs(rnorm(1,0.5,0.2)),
         mu5=abs(rnorm(1,0.5,0.2)),mu6=abs(rnorm(1,0.5,0.2)),
         ks1=abs(rnorm(1,0.5,0.2)),ks2=abs(rnorm(1,0.5,0.2)),
         ks3=abs(rnorm(1,0.5,0.2)),ks4=abs(rnorm(1,0.5,0.2)),
         ks5=abs(rnorm(1,0.5,0.2)),ks6=abs(rnorm(1,0.5,0.2)),
         k1=abs(rnorm(1,0.5,0.2)),k2=abs(rnorm(1,0.5,0.2)),
         k3=abs(rnorm(1,0.5,0.2)),yc1=abs(rnorm(1,0.5,0.2)),
         yc2=abs(rnorm(1,0.5,0.2)),yc3=abs(rnorm(1,0.5,0.2)),
         yc4=abs(rnorm(1,0.5,0.2)),yc5=abs(rnorm(1,0.5,0.2)),
         yc6=abs(rnorm(1,0.5,0.2)),yc7=abs(rnorm(1,0.5,0.2)),
         yc8=abs(rnorm(1,0.5,0.2)),yc9=abs(rnorm(1,0.5,0.2)),
         yc10=abs(rnorm(1,0.5,0.2)),yc11=abs(rnorm(1,0.5,0.2)),
         yc12=abs(rnorm(1,0.5,0.2)),yc13=abs(rnorm(1,0.5,0.2)),
         yc14=abs(rnorm(1,0.5,0.2)),yc15=abs(rnorm(1,0.5,0.2)),
         sigma=abs(rnorm(1,0.5,0.2)))
  }
}
if (args[2]==35){
  ini = function(){
    list(mu1=abs(rnorm(1,0.5,0.2)),mu2=abs(rnorm(1,0.5,0.2)),
         mu3=abs(rnorm(1,0.5,0.2)),mu4=abs(rnorm(1,0.5,0.2)),
         mu5=abs(rnorm(1,0.5,0.2)),mu6=abs(rnorm(1,0.5,0.2)),
         ks1=abs(rnorm(1,0.5,0.2)),ks2=abs(rnorm(1,0.5,0.2)),
         ks3=abs(rnorm(1,0.5,0.2)),ks4=abs(rnorm(1,0.5,0.2)),
         ks5=abs(rnorm(1,0.5,0.2)),ks6=abs(rnorm(1,0.5,0.2)),
         k1=abs(rnorm(1,0.5,0.2)),k2=abs(rnorm(1,0.5,0.2)),
         k3=abs(rnorm(1,0.5,0.2)),yc1=abs(rnorm(1,0.5,0.2)),
         yc2=abs(rnorm(1,0.5,0.2)),yc3=abs(rnorm(1,0.5,0.2)),
         yc4=abs(rnorm(1,0.5,0.2)),yc5=abs(rnorm(1,0.5,0.2)),
         yc6=abs(rnorm(1,0.5,0.2)),yc7=abs(rnorm(1,0.5,0.2)),
         yc8=abs(rnorm(1,0.5,0.2)),yc9=abs(rnorm(1,0.5,0.2)),
         yc10=abs(rnorm(1,0.5,0.2)),yc11=abs(rnorm(1,0.5,0.2)),
         yc12=abs(rnorm(1,0.5,0.2)),yc13=abs(rnorm(1,0.5,0.2)),
         yc14=abs(rnorm(1,0.5,0.2)),sigma=abs(rnorm(1,0.5,0.2)))
  }
}
if (args[2]==45){
  ini = function(){
    list(mu1=abs(rnorm(1,0.5,0.2)),mu2=abs(rnorm(1,0.5,0.2)),
         mu3=abs(rnorm(1,0.5,0.2)),mu4=abs(rnorm(1,0.5,0.2)),
         mu5=abs(rnorm(1,0.5,0.2)),mu6=abs(rnorm(1,0.5,0.2)),
         mu7=abs(rnorm(1,0.5,0.2)),ks1=abs(rnorm(1,0.5,0.2)),
         ks2=abs(rnorm(1,0.5,0.2)),ks3=abs(rnorm(1,0.5,0.2)),
         ks4=abs(rnorm(1,0.5,0.2)),ks5=abs(rnorm(1,0.5,0.2)),
         ks6=abs(rnorm(1,0.5,0.2)),ks7=abs(rnorm(1,0.5,0.2)),
         k1=abs(rnorm(1,0.5,0.2)),k2=abs(rnorm(1,0.5,0.2)),
         k3=abs(rnorm(1,0.5,0.2)),yc1=abs(rnorm(1,0.5,0.2)),
         yc2=abs(rnorm(1,0.5,0.2)),yc3=abs(rnorm(1,0.5,0.2)),
         yc4=abs(rnorm(1,0.5,0.2)),yc5=abs(rnorm(1,0.5,0.2)),
         yc6=abs(rnorm(1,0.5,0.2)),yc7=abs(rnorm(1,0.5,0.2)),
         yc8=abs(rnorm(1,0.5,0.2)),yc9=abs(rnorm(1,0.5,0.2)),
         yc10=abs(rnorm(1,0.5,0.2)),yc11=abs(rnorm(1,0.5,0.2)),
         yc12=abs(rnorm(1,0.5,0.2)),yc13=abs(rnorm(1,0.5,0.2)),
         yc14=abs(rnorm(1,0.5,0.2)),sigma=abs(rnorm(1,0.5,0.2)))
  }
}
if (args[2]==123){
  ini = function(){
    list(mu1=abs(rnorm(1,0.5,0.2)),mu2=abs(rnorm(1,0.5,0.2)),
         mu3=abs(rnorm(1,0.5,0.2)),mu4=abs(rnorm(1,0.5,0.2)),
         mu5=abs(rnorm(1,0.5,0.2)),mu6=abs(rnorm(1,0.5,0.2)),
         ks1=abs(rnorm(1,0.5,0.2)),ks2=abs(rnorm(1,0.5,0.2)),
         ks3=abs(rnorm(1,0.5,0.2)),ks4=abs(rnorm(1,0.5,0.2)),
         ks5=abs(rnorm(1,0.5,0.2)),ks6=abs(rnorm(1,0.5,0.2)),
         k1=abs(rnorm(1,0.5,0.2)),k2=abs(rnorm(1,0.5,0.2)),
         k3=abs(rnorm(1,0.5,0.2)),yc1=abs(rnorm(1,0.5,0.2)),
         yc2=abs(rnorm(1,0.5,0.2)),yc3=abs(rnorm(1,0.5,0.2)),
         yc4=abs(rnorm(1,0.5,0.2)),yc5=abs(rnorm(1,0.5,0.2)),
         yc6=abs(rnorm(1,0.5,0.2)),yc7=abs(rnorm(1,0.5,0.2)),
         yc8=abs(rnorm(1,0.5,0.2)),yc9=abs(rnorm(1,0.5,0.2)),
         yc10=abs(rnorm(1,0.5,0.2)),yc11=abs(rnorm(1,0.5,0.2)),
         yc12=abs(rnorm(1,0.5,0.2)),yc13=abs(rnorm(1,0.5,0.2)),
         yc14=abs(rnorm(1,0.5,0.2)),yc15=abs(rnorm(1,0.5,0.2)),
         yc16=abs(rnorm(1,0.5,0.2)),y17=abs(rnorm(1,0.5,0.2)),
         yc18=abs(rnorm(1,0.5,0.2)),ev1=abs(rnorm(1,0.5,0.2)),
         ev2=abs(rnorm(1,0.5,0.2)),ev3=abs(rnorm(1,0.5,0.2)),
         sigma=abs(rnorm(1,0.5,0.2)))
  }
}
if (args[2]==124){
  ini = function(){
    list(mu1=abs(rnorm(1,0.5,0.2)),mu2=abs(rnorm(1,0.5,0.2)),
         mu3=abs(rnorm(1,0.5,0.2)),mu4=abs(rnorm(1,0.5,0.2)),
         mu5=abs(rnorm(1,0.5,0.2)),mu6=abs(rnorm(1,0.5,0.2)),
         mu7=abs(rnorm(1,0.5,0.2)),ks1=abs(rnorm(1,0.5,0.2)),
         ks2=abs(rnorm(1,0.5,0.2)),ks3=abs(rnorm(1,0.5,0.2)),
         ks4=abs(rnorm(1,0.5,0.2)),ks5=abs(rnorm(1,0.5,0.2)),
         ks6=abs(rnorm(1,0.5,0.2)),ks7=abs(rnorm(1,0.5,0.2)),
         k1=abs(rnorm(1,0.5,0.2)),k2=abs(rnorm(1,0.5,0.2)),
         k3=abs(rnorm(1,0.5,0.2)),yc1=abs(rnorm(1,0.5,0.2)),
         yc2=abs(rnorm(1,0.5,0.2)),yc3=abs(rnorm(1,0.5,0.2)),
         yc4=abs(rnorm(1,0.5,0.2)),yc5=abs(rnorm(1,0.5,0.2)),
         yc6=abs(rnorm(1,0.5,0.2)),yc7=abs(rnorm(1,0.5,0.2)),
         yc8=abs(rnorm(1,0.5,0.2)),yc9=abs(rnorm(1,0.5,0.2)),
         yc10=abs(rnorm(1,0.5,0.2)),yc11=abs(rnorm(1,0.5,0.2)),
         yc12=abs(rnorm(1,0.5,0.2)),yc13=abs(rnorm(1,0.5,0.2)),
         yc14=abs(rnorm(1,0.5,0.2)),yc15=abs(rnorm(1,0.5,0.2)),
         yc16=abs(rnorm(1,0.5,0.2)),yc17=abs(rnorm(1,0.5,0.2)),
         yc18=abs(rnorm(1,0.5,0.2)),ev1=abs(rnorm(1,0.5,0.2)),
         ev2=abs(rnorm(1,0.5,0.2)),ev3=abs(rnorm(1,0.5,0.2)),
         sigma=abs(rnorm(1,0.5,0.2)))
  }
}

if (args[2]==125){
  ini = function(){
    list(mu1=abs(rnorm(1,0.5,0.2)),mu2=abs(rnorm(1,0.5,0.2)),
         mu3=abs(rnorm(1,0.5,0.2)),mu4=abs(rnorm(1,0.5,0.2)),
         mu5=abs(rnorm(1,0.5,0.2)),mu6=abs(rnorm(1,0.5,0.2)),
         mu7=abs(rnorm(1,0.5,0.2)),ks1=abs(rnorm(1,0.5,0.2)),
         ks2=abs(rnorm(1,0.5,0.2)),ks3=abs(rnorm(1,0.5,0.2)),
         ks4=abs(rnorm(1,0.5,0.2)),ks5=abs(rnorm(1,0.5,0.2)),
         ks6=abs(rnorm(1,0.5,0.2)),ks7=abs(rnorm(1,0.5,0.2)),
         k1=abs(rnorm(1,0.5,0.2)),k2=abs(rnorm(1,0.5,0.2)),
         k3=abs(rnorm(1,0.5,0.2)),yc1=abs(rnorm(1,0.5,0.2)),
         yc2=abs(rnorm(1,0.5,0.2)),yc3=abs(rnorm(1,0.5,0.2)),
         yc4=abs(rnorm(1,0.5,0.2)),yc5=abs(rnorm(1,0.5,0.2)),
         yc6=abs(rnorm(1,0.5,0.2)),yc7=abs(rnorm(1,0.5,0.2)),
         yc8=abs(rnorm(1,0.5,0.2)),yc9=abs(rnorm(1,0.5,0.2)),
         yc10=abs(rnorm(1,0.5,0.2)),yc11=abs(rnorm(1,0.5,0.2)),
         yc12=abs(rnorm(1,0.5,0.2)),yc13=abs(rnorm(1,0.5,0.2)),
         yc14=abs(rnorm(1,0.5,0.2)),yc15=abs(rnorm(1,0.5,0.2)),
         yc16=abs(rnorm(1,0.5,0.2)),yc17=abs(rnorm(1,0.5,0.2)),
         ev1=abs(rnorm(1,0.5,0.2)),ev2=abs(rnorm(1,0.5,0.2)),
         ev3=abs(rnorm(1,0.5,0.2)),sigma=abs(rnorm(1,0.5,0.2)))
  }
}
if (args[2]==134){
  ini = function(){
    list(mu1=abs(rnorm(1,0.5,0.2)),mu2=abs(rnorm(1,0.5,0.2)),
         mu3=abs(rnorm(1,0.5,0.2)),mu4=abs(rnorm(1,0.5,0.2)),
         mu5=abs(rnorm(1,0.5,0.2)),mu6=abs(rnorm(1,0.5,0.2)),
         ks1=abs(rnorm(1,0.5,0.2)),ks2=abs(rnorm(1,0.5,0.2)),
         ks3=abs(rnorm(1,0.5,0.2)),ks4=abs(rnorm(1,0.5,0.2)),
         ks5=abs(rnorm(1,0.5,0.2)),ks6=abs(rnorm(1,0.5,0.2)),
         k1=abs(rnorm(1,0.5,0.2)),k2=abs(rnorm(1,0.5,0.2)),
         k3=abs(rnorm(1,0.5,0.2)),yc1=abs(rnorm(1,0.5,0.2)),
         yc2=abs(rnorm(1,0.5,0.2)),yc3=abs(rnorm(1,0.5,0.2)),
         yc4=abs(rnorm(1,0.5,0.2)),yc5=abs(rnorm(1,0.5,0.2)),
         yc6=abs(rnorm(1,0.5,0.2)),yc7=abs(rnorm(1,0.5,0.2)),
         yc8=abs(rnorm(1,0.5,0.2)),yc9=abs(rnorm(1,0.5,0.2)),
         yc10=abs(rnorm(1,0.5,0.2)),yc11=abs(rnorm(1,0.5,0.2)),
         yc12=abs(rnorm(1,0.5,0.2)),yc13=abs(rnorm(1,0.5,0.2)),
         yc14=abs(rnorm(1,0.5,0.2)),yc15=abs(rnorm(1,0.5,0.2)),
         ev1=abs(rnorm(1,0.5,0.2)),ev2=abs(rnorm(1,0.5,0.2)),
         ev3=abs(rnorm(1,0.5,0.2)),sigma=abs(rnorm(1,0.5,0.2)))
  }
}
if (args[2]==135){
  ini = function(){
    list(mu1=abs(rnorm(1,0.5,0.2)),mu2=abs(rnorm(1,0.5,0.2)),
         mu3=abs(rnorm(1,0.5,0.2)),mu4=abs(rnorm(1,0.5,0.2)),
         mu5=abs(rnorm(1,0.5,0.2)),mu6=abs(rnorm(1,0.5,0.2)),
         ks1=abs(rnorm(1,0.5,0.2)),ks2=abs(rnorm(1,0.5,0.2)),
         ks3=abs(rnorm(1,0.5,0.2)),ks4=abs(rnorm(1,0.5,0.2)),
         ks5=abs(rnorm(1,0.5,0.2)),ks6=abs(rnorm(1,0.5,0.2)),
         k1=abs(rnorm(1,0.5,0.2)),k2=abs(rnorm(1,0.5,0.2)),
         k3=abs(rnorm(1,0.5,0.2)),yc1=abs(rnorm(1,0.5,0.2)),
         yc2=abs(rnorm(1,0.5,0.2)),yc3=abs(rnorm(1,0.5,0.2)),
         yc4=abs(rnorm(1,0.5,0.2)),yc5=abs(rnorm(1,0.5,0.2)),
         yc6=abs(rnorm(1,0.5,0.2)),yc7=abs(rnorm(1,0.5,0.2)),
         yc8=abs(rnorm(1,0.5,0.2)),yc9=abs(rnorm(1,0.5,0.2)),
         yc10=abs(rnorm(1,0.5,0.2)),yc11=abs(rnorm(1,0.5,0.2)),
         yc12=abs(rnorm(1,0.5,0.2)),yc13=abs(rnorm(1,0.5,0.2)),
         yc14=abs(rnorm(1,0.5,0.2)),ev1=abs(rnorm(1,0.5,0.2)),
         ev2=abs(rnorm(1,0.5,0.2)),ev3=abs(rnorm(1,0.5,0.2)),
         sigma=abs(rnorm(1,0.5,0.2)))
  }
}
if (args[2]==145){
  ini = function(){
    list(mu1=abs(rnorm(1,0.5,0.2)),mu2=abs(rnorm(1,0.5,0.2)),
         mu3=abs(rnorm(1,0.5,0.2)),mu4=abs(rnorm(1,0.5,0.2)),
         mu5=abs(rnorm(1,0.5,0.2)),mu6=abs(rnorm(1,0.5,0.2)),
         mu7=abs(rnorm(1,0.5,0.2)),ks1=abs(rnorm(1,0.5,0.2)),
         ks2=abs(rnorm(1,0.5,0.2)),ks3=abs(rnorm(1,0.5,0.2)),
         ks4=abs(rnorm(1,0.5,0.2)),ks5=abs(rnorm(1,0.5,0.2)),
         ks6=abs(rnorm(1,0.5,0.2)),ks7=abs(rnorm(1,0.5,0.2)),
         k1=abs(rnorm(1,0.5,0.2)),k2=abs(rnorm(1,0.5,0.2)),
         k3=abs(rnorm(1,0.5,0.2)),yc1=abs(rnorm(1,0.5,0.2)),
         yc2=abs(rnorm(1,0.5,0.2)),yc3=abs(rnorm(1,0.5,0.2)),
         yc4=abs(rnorm(1,0.5,0.2)),yc5=abs(rnorm(1,0.5,0.2)),
         yc6=abs(rnorm(1,0.5,0.2)),yc7=abs(rnorm(1,0.5,0.2)),
         yc8=abs(rnorm(1,0.5,0.2)),yc9=abs(rnorm(1,0.5,0.2)),
         yc10=abs(rnorm(1,0.5,0.2)),yc11=abs(rnorm(1,0.5,0.2)),
         yc13=abs(rnorm(1,0.5,0.2)),yc14=abs(rnorm(1,0.5,0.2)),
         ev1=abs(rnorm(1,0.5,0.2)),ev2=abs(rnorm(1,0.5,0.2)),
         ev3=abs(rnorm(1,0.5,0.2)),sigma=abs(rnorm(1,0.5,0.2)))
  }
}
if (args[2]==234){
  ini = function(){
    list(mu1=abs(rnorm(1,0.5,0.2)),mu2=abs(rnorm(1,0.5,0.2)),
         mu3=abs(rnorm(1,0.5,0.2)),mu4=abs(rnorm(1,0.5,0.2)),
         mu5=abs(rnorm(1,0.5,0.2)),mu6=abs(rnorm(1,0.5,0.2)),
         mu7=abs(rnorm(1,0.5,0.2)),ks1=abs(rnorm(1,0.5,0.2)),
         ks2=abs(rnorm(1,0.5,0.2)),ks3=abs(rnorm(1,0.5,0.2)),
         ks4=abs(rnorm(1,0.5,0.2)),ks5=abs(rnorm(1,0.5,0.2)),
         ks6=abs(rnorm(1,0.5,0.2)),ks7=abs(rnorm(1,0.5,0.2)),
         k1=abs(rnorm(1,0.5,0.2)),k2=abs(rnorm(1,0.5,0.2)),
         k3=abs(rnorm(1,0.5,0.2)),yc1=abs(rnorm(1,0.5,0.2)),
         yc2=abs(rnorm(1,0.5,0.2)),yc3=abs(rnorm(1,0.5,0.2)),
         yc4=abs(rnorm(1,0.5,0.2)),yc5=abs(rnorm(1,0.5,0.2)),
         yc6=abs(rnorm(1,0.5,0.2)),yc7=abs(rnorm(1,0.5,0.2)),
         yc8=abs(rnorm(1,0.5,0.2)),yc9=abs(rnorm(1,0.5,0.2)),
         yc10=abs(rnorm(1,0.5,0.2)),yc11=abs(rnorm(1,0.5,0.2)),
         yc12=abs(rnorm(1,0.5,0.2)),yc13=abs(rnorm(1,0.5,0.2)),
         yc14=abs(rnorm(1,0.5,0.2)),yc15=abs(rnorm(1,0.5,0.2)),
         yc16=abs(rnorm(1,0.5,0.2)),yc17=abs(rnorm(1,0.5,0.2)),
         yc18=abs(rnorm(1,0.5,0.2)),yc19=abs(rnorm(1,0.5,0.2)),
         yc20=abs(rnorm(1,0.5,0.2)),sigma=abs(rnorm(1,0.5,0.2)))
  }
}
if (args[2]==235){
  ini = function(){
    list(mu1=abs(rnorm(1,0.5,0.2)),mu2=abs(rnorm(1,0.5,0.2)),
         mu3=abs(rnorm(1,0.5,0.2)),mu4=abs(rnorm(1,0.5,0.2)),
         mu5=abs(rnorm(1,0.5,0.2)),mu6=abs(rnorm(1,0.5,0.2)),
         mu7=abs(rnorm(1,0.5,0.2)),ks1=abs(rnorm(1,0.5,0.2)),
         ks2=abs(rnorm(1,0.5,0.2)),ks3=abs(rnorm(1,0.5,0.2)),
         ks4=abs(rnorm(1,0.5,0.2)),ks5=abs(rnorm(1,0.5,0.2)),
         ks6=abs(rnorm(1,0.5,0.2)),ks7=abs(rnorm(1,0.5,0.2)),
         k1=abs(rnorm(1,0.5,0.2)),k2=abs(rnorm(1,0.5,0.2)),
         k3=abs(rnorm(1,0.5,0.2)),yc1=abs(rnorm(1,0.5,0.2)),
         yc2=abs(rnorm(1,0.5,0.2)),yc3=abs(rnorm(1,0.5,0.2)),
         yc4=abs(rnorm(1,0.5,0.2)),yc5=abs(rnorm(1,0.5,0.2)),
         yc6=abs(rnorm(1,0.5,0.2)),yc7=abs(rnorm(1,0.5,0.2)),
         yc8=abs(rnorm(1,0.5,0.2)),yc9=abs(rnorm(1,0.5,0.2)),
         yc10=abs(rnorm(1,0.5,0.2)),yc11=abs(rnorm(1,0.5,0.2)),
         yc12=abs(rnorm(1,0.5,0.2)),yc13=abs(rnorm(1,0.5,0.2)),
         yc14=abs(rnorm(1,0.5,0.2)),yc15=abs(rnorm(1,0.5,0.2)),
         yc16=abs(rnorm(1,0.5,0.2)),yc17=abs(rnorm(1,0.5,0.2)),
         yc18=abs(rnorm(1,0.5,0.2)),yc19=abs(rnorm(1,0.5,0.2)),
         sigma=abs(rnorm(1,0.5,0.2)))
  }
}
if (args[2]==245){
  ini = function(){
    list(mu1=abs(rnorm(1,0.5,0.2)),mu2=abs(rnorm(1,0.5,0.2)),
         mu3=abs(rnorm(1,0.5,0.2)),mu4=abs(rnorm(1,0.5,0.2)),
         mu5=abs(rnorm(1,0.5,0.2)),mu6=abs(rnorm(1,0.5,0.2)),
         mu7=abs(rnorm(1,0.5,0.2)),mu8=abs(rnorm(1,0.5,0.2)),
         ks1=abs(rnorm(1,0.5,0.2)),ks2=abs(rnorm(1,0.5,0.2)),
         ks3=abs(rnorm(1,0.5,0.2)),ks4=abs(rnorm(1,0.5,0.2)),
         ks5=abs(rnorm(1,0.5,0.2)),ks6=abs(rnorm(1,0.5,0.2)),
         ks7=abs(rnorm(1,0.5,0.2)),ks8=abs(rnorm(1,0.5,0.2)),
         k1=abs(rnorm(1,0.5,0.2)),k2=abs(rnorm(1,0.5,0.2)),
         k3=abs(rnorm(1,0.5,0.2)),yc1=abs(rnorm(1,0.5,0.2)),
         yc2=abs(rnorm(1,0.5,0.2)),yc3=abs(rnorm(1,0.5,0.2)),
         yc4=abs(rnorm(1,0.5,0.2)),yc5=abs(rnorm(1,0.5,0.2)),
         yc6=abs(rnorm(1,0.5,0.2)),yc7=abs(rnorm(1,0.5,0.2)),
         yc8=abs(rnorm(1,0.5,0.2)),yc9=abs(rnorm(1,0.5,0.2)),
         yc10=abs(rnorm(1,0.5,0.2)),yc11=abs(rnorm(1,0.5,0.2)),
         yc12=abs(rnorm(1,0.5,0.2)),yc13=abs(rnorm(1,0.5,0.2)),
         yc14=abs(rnorm(1,0.5,0.2)),yc15=abs(rnorm(1,0.5,0.2)),
         yc16=abs(rnorm(1,0.5,0.2)),yc17=abs(rnorm(1,0.5,0.2)),
         yc18=abs(rnorm(1,0.5,0.2)),yc19=abs(rnorm(1,0.5,0.2)),
         sigma=abs(rnorm(1,0.5,0.2)))
  }
}
if (args[2]==345){
  ini = function(){
    list(mu1=abs(rnorm(1,0.5,0.2)),mu2=abs(rnorm(1,0.5,0.2)),
         mu3=abs(rnorm(1,0.5,0.2)),mu4=abs(rnorm(1,0.5,0.2)),
         mu5=abs(rnorm(1,0.5,0.2)),mu6=abs(rnorm(1,0.5,0.2)),
         mu7=abs(rnorm(1,0.5,0.2)),ks1=abs(rnorm(1,0.5,0.2)),
         ks2=abs(rnorm(1,0.5,0.2)),ks3=abs(rnorm(1,0.5,0.2)),
         ks4=abs(rnorm(1,0.5,0.2)),ks5=abs(rnorm(1,0.5,0.2)),
         ks6=abs(rnorm(1,0.5,0.2)),ks7=abs(rnorm(1,0.5,0.2)),
         k1=abs(rnorm(1,0.5,0.2)),k2=abs(rnorm(1,0.5,0.2)),
         k3=abs(rnorm(1,0.5,0.2)),yc1=abs(rnorm(1,0.5,0.2)),
         yc2=abs(rnorm(1,0.5,0.2)),yc3=abs(rnorm(1,0.5,0.2)),
         yc4=abs(rnorm(1,0.5,0.2)),yc5=abs(rnorm(1,0.5,0.2)),
         yc6=abs(rnorm(1,0.5,0.2)),yc7=abs(rnorm(1,0.5,0.2)),
         yc8=abs(rnorm(1,0.5,0.2)),yc9=abs(rnorm(1,0.5,0.2)),
         yc10=abs(rnorm(1,0.5,0.2)),yc11=abs(rnorm(1,0.5,0.2)),
         yc12=abs(rnorm(1,0.5,0.2)),yc13=abs(rnorm(1,0.5,0.2)),
         yc14=abs(rnorm(1,0.5,0.2)),yc15=abs(rnorm(1,0.5,0.2)),
         yc16=abs(rnorm(1,0.5,0.2)),sigma=abs(rnorm(1,0.5,0.2)))
  }
}
if (args[2]==1234){
  ini = function(){
    list(mu1=abs(rnorm(1,0.5,0.2)),mu2=abs(rnorm(1,0.5,0.2)),
         mu3=abs(rnorm(1,0.5,0.2)),mu4=abs(rnorm(1,0.5,0.2)),
         mu5=abs(rnorm(1,0.5,0.2)),mu6=abs(rnorm(1,0.5,0.2)),
         mu7=abs(rnorm(1,0.5,0.2)),ks1=abs(rnorm(1,0.5,0.2)),
         ks2=abs(rnorm(1,0.5,0.2)),ks3=abs(rnorm(1,0.5,0.2)),
         ks4=abs(rnorm(1,0.5,0.2)),ks5=abs(rnorm(1,0.5,0.2)),
         ks6=abs(rnorm(1,0.5,0.2)),ks7=abs(rnorm(1,0.5,0.2)),
         k1=abs(rnorm(1,0.5,0.2)),k2=abs(rnorm(1,0.5,0.2)),
         k3=abs(rnorm(1,0.5,0.2)),yc1=abs(rnorm(1,0.5,0.2)),
         yc2=abs(rnorm(1,0.5,0.2)),yc3=abs(rnorm(1,0.5,0.2)),
         yc4=abs(rnorm(1,0.5,0.2)),yc5=abs(rnorm(1,0.5,0.2)),
         yc6=abs(rnorm(1,0.5,0.2)),yc7=abs(rnorm(1,0.5,0.2)),
         yc8=abs(rnorm(1,0.5,0.2)),yc9=abs(rnorm(1,0.5,0.2)),
         yc10=abs(rnorm(1,0.5,0.2)),yc11=abs(rnorm(1,0.5,0.2)),
         yc12=abs(rnorm(1,0.5,0.2)),yc13=abs(rnorm(1,0.5,0.2)),
         yc14=abs(rnorm(1,0.5,0.2)),yc15=abs(rnorm(1,0.5,0.2)),
         yc16=abs(rnorm(1,0.5,0.2)),yc17=abs(rnorm(1,0.5,0.2)),
         yc18=abs(rnorm(1,0.5,0.2)),yc19=abs(rnorm(1,0.5,0.2)),
         yc20=abs(rnorm(1,0.5,0.2)),ev1=abs(rnorm(1,0.5,0.2)),
         ev2=abs(rnorm(1,0.5,0.2)),ev3=abs(rnorm(1,0.5,0.2)),
         sigma=abs(rnorm(1,0.5,0.2)))
  }
}
if (args[2]==1235){
  ini = function(){
    list(mu1=abs(rnorm(1,0.5,0.2)),mu2=abs(rnorm(1,0.5,0.2)),
         mu3=abs(rnorm(1,0.5,0.2)),mu4=abs(rnorm(1,0.5,0.2)),
         mu5=abs(rnorm(1,0.5,0.2)),mu6=abs(rnorm(1,0.5,0.2)),
         mu7=abs(rnorm(1,0.5,0.2)),ks1=abs(rnorm(1,0.5,0.2)),
         ks2=abs(rnorm(1,0.5,0.2)),ks3=abs(rnorm(1,0.5,0.2)),
         ks4=abs(rnorm(1,0.5,0.2)),ks5=abs(rnorm(1,0.5,0.2)),
         ks6=abs(rnorm(1,0.5,0.2)),ks7=abs(rnorm(1,0.5,0.2)),
         k1=abs(rnorm(1,0.5,0.2)),k2=abs(rnorm(1,0.5,0.2)),
         k3=abs(rnorm(1,0.5,0.2)),yc1=abs(rnorm(1,0.5,0.2)),
         yc2=abs(rnorm(1,0.5,0.2)),yc3=abs(rnorm(1,0.5,0.2)),
         yc4=abs(rnorm(1,0.5,0.2)),yc5=abs(rnorm(1,0.5,0.2)),
         yc6=abs(rnorm(1,0.5,0.2)),yc7=abs(rnorm(1,0.5,0.2)),
         yc8=abs(rnorm(1,0.5,0.2)),yc9=abs(rnorm(1,0.5,0.2)),
         yc10=abs(rnorm(1,0.5,0.2)),yc11=abs(rnorm(1,0.5,0.2)),
         yc12=abs(rnorm(1,0.5,0.2)),yc13=abs(rnorm(1,0.5,0.2)),
         yc14=abs(rnorm(1,0.5,0.2)),yc15=abs(rnorm(1,0.5,0.2)),
         yc16=abs(rnorm(1,0.5,0.2)),yc17=abs(rnorm(1,0.5,0.2)),
         yc18=abs(rnorm(1,0.5,0.2)),yc19=abs(rnorm(1,0.5,0.2)),
         ev1=abs(rnorm(1,0.5,0.2)),ev2=abs(rnorm(1,0.5,0.2)),
         ev3=abs(rnorm(1,0.5,0.2)),sigma=abs(rnorm(1,0.5,0.2)))
  }
}
if (args[2]==1245){
  ini = function(){
    list(mu1=abs(rnorm(1,0.5,0.2)),mu2=abs(rnorm(1,0.5,0.2)),
         mu3=abs(rnorm(1,0.5,0.2)),mu4=abs(rnorm(1,0.5,0.2)),
         mu5=abs(rnorm(1,0.5,0.2)),mu6=abs(rnorm(1,0.5,0.2)),
         mu7=abs(rnorm(1,0.5,0.2)),mu8=abs(rnorm(1,0.5,0.2)),
         ks1=abs(rnorm(1,0.5,0.2)),ks2=abs(rnorm(1,0.5,0.2)),
         ks3=abs(rnorm(1,0.5,0.2)),ks4=abs(rnorm(1,0.5,0.2)),
         ks5=abs(rnorm(1,0.5,0.2)),ks6=abs(rnorm(1,0.5,0.2)),
         ks7=abs(rnorm(1,0.5,0.2)),ks8=abs(rnorm(1,0.5,0.2)),
         k1=abs(rnorm(1,0.5,0.2)),k2=abs(rnorm(1,0.5,0.2)),
         k3=abs(rnorm(1,0.5,0.2)),yc1=abs(rnorm(1,0.5,0.2)),
         yc2=abs(rnorm(1,0.5,0.2)),yc3=abs(rnorm(1,0.5,0.2)),
         yc4=abs(rnorm(1,0.5,0.2)),yc5=abs(rnorm(1,0.5,0.2)),
         yc6=abs(rnorm(1,0.5,0.2)),yc7=abs(rnorm(1,0.5,0.2)),
         yc8=abs(rnorm(1,0.5,0.2)),yc9=abs(rnorm(1,0.5,0.2)),
         yc10=abs(rnorm(1,0.5,0.2)),yc11=abs(rnorm(1,0.5,0.2)),
         yc12=abs(rnorm(1,0.5,0.2)),yc13=abs(rnorm(1,0.5,0.2)),
         yc14=abs(rnorm(1,0.5,0.2)),yc15=abs(rnorm(1,0.5,0.2)),
         yc16=abs(rnorm(1,0.5,0.2)),yc17=abs(rnorm(1,0.5,0.2)),
         yc18=abs(rnorm(1,0.5,0.2)),yc19=abs(rnorm(1,0.5,0.2)),
         ev1=abs(rnorm(1,0.5,0.2)),ev2=abs(rnorm(1,0.5,0.2)),
         ev3=abs(rnorm(1,0.5,0.2)),sigma=abs(rnorm(1,0.5,0.2)))
  }
}

if (args[2]==1345){
  ini = function(){
    list(mu1=abs(rnorm(1,0.5,0.2)),mu2=abs(rnorm(1,0.5,0.2)),
         mu3=abs(rnorm(1,0.5,0.2)),mu4=abs(rnorm(1,0.5,0.2)),
         mu5=abs(rnorm(1,0.5,0.2)),mu6=abs(rnorm(1,0.5,0.2)),
         mu7=abs(rnorm(1,0.5,0.2)),ks1=abs(rnorm(1,0.5,0.2)),
         ks2=abs(rnorm(1,0.5,0.2)),ks3=abs(rnorm(1,0.5,0.2)),
         ks4=abs(rnorm(1,0.5,0.2)),ks5=abs(rnorm(1,0.5,0.2)),
         ks6=abs(rnorm(1,0.5,0.2)),ks7=abs(rnorm(1,0.5,0.2)),
         k1=abs(rnorm(1,0.5,0.2)),k2=abs(rnorm(1,0.5,0.2)),
         k3=abs(rnorm(1,0.5,0.2)),yc1=abs(rnorm(1,0.5,0.2)),
         yc2=abs(rnorm(1,0.5,0.2)),yc3=abs(rnorm(1,0.5,0.2)),
         yc4=abs(rnorm(1,0.5,0.2)),yc5=abs(rnorm(1,0.5,0.2)),
         yc6=abs(rnorm(1,0.5,0.2)),yc7=abs(rnorm(1,0.5,0.2)),
         yc8=abs(rnorm(1,0.5,0.2)),yc9=abs(rnorm(1,0.5,0.2)),
         yc10=abs(rnorm(1,0.5,0.2)),yc11=abs(rnorm(1,0.5,0.2)),
         yc12=abs(rnorm(1,0.5,0.2)),yc13=abs(rnorm(1,0.5,0.2)),
         yc14=abs(rnorm(1,0.5,0.2)),yc15=abs(rnorm(1,0.5,0.2)),
         yc16=abs(rnorm(1,0.5,0.2)),ev1=abs(rnorm(1,0.5,0.2)),
         ev2=abs(rnorm(1,0.5,0.2)),ev3=abs(rnorm(1,0.5,0.2)),
         sigma=abs(rnorm(1,0.5,0.2)))
  }
}

if (args[2]==2345){
  ini = function(){
    list(mu1=abs(rnorm(1,0.5,0.2)),mu2=abs(rnorm(1,0.5,0.2)),
         mu3=abs(rnorm(1,0.5,0.2)),mu4=abs(rnorm(1,0.5,0.2)),
         mu5=abs(rnorm(1,0.5,0.2)),mu6=abs(rnorm(1,0.5,0.2)),
         mu7=abs(rnorm(1,0.5,0.2)),mu8=abs(rnorm(1,0.5,0.2)),
         ks1=abs(rnorm(1,0.5,0.2)),ks2=abs(rnorm(1,0.5,0.2)),
         ks3=abs(rnorm(1,0.5,0.2)),ks4=abs(rnorm(1,0.5,0.2)),
         ks5=abs(rnorm(1,0.5,0.2)),ks6=abs(rnorm(1,0.5,0.2)),
         ks7=abs(rnorm(1,0.5,0.2)),ks8=abs(rnorm(1,0.5,0.2)),
         k1=abs(rnorm(1,0.5,0.2)),k2=abs(rnorm(1,0.5,0.2)),
         k3=abs(rnorm(1,0.5,0.2)),yc1=abs(rnorm(1,0.5,0.2)),
         yc2=abs(rnorm(1,0.5,0.2)),yc3=abs(rnorm(1,0.5,0.2)),
         yc4=abs(rnorm(1,0.5,0.2)),yc5=abs(rnorm(1,0.5,0.2)),
         yc6=abs(rnorm(1,0.5,0.2)),yc7=abs(rnorm(1,0.5,0.2)),
         yc8=abs(rnorm(1,0.5,0.2)),yc9=abs(rnorm(1,0.5,0.2)),
         yc10=abs(rnorm(1,0.5,0.2)),yc11=abs(rnorm(1,0.5,0.2)),
         yc12=abs(rnorm(1,0.5,0.2)),yc13=abs(rnorm(1,0.5,0.2)),
         yc14=abs(rnorm(1,0.5,0.2)),yc15=abs(rnorm(1,0.5,0.2)),
         yc16=abs(rnorm(1,0.5,0.2)),yc17=abs(rnorm(1,0.5,0.2)),
         yc18=abs(rnorm(1,0.5,0.2)),yc19=abs(rnorm(1,0.5,0.2)),
         yc20=abs(rnorm(1,0.5,0.2)),yc21=abs(rnorm(1,0.5,0.2)),
         sigma=abs(rnorm(1,0.5,0.2)))
  }
}

if (args[2]==12345){
  ini = function(){
    list(mu1=abs(rnorm(1,0.5,0.2)),mu2=abs(rnorm(1,0.5,0.2)),
         mu3=abs(rnorm(1,0.5,0.2)),mu4=abs(rnorm(1,0.5,0.2)),
         mu5=abs(rnorm(1,0.5,0.2)),mu6=abs(rnorm(1,0.5,0.2)),
         mu7=abs(rnorm(1,0.5,0.2)),mu8=abs(rnorm(1,0.5,0.2)),
         ks1=abs(rnorm(1,0.5,0.2)),ks2=abs(rnorm(1,0.5,0.2)),
         ks3=abs(rnorm(1,0.5,0.2)),ks4=abs(rnorm(1,0.5,0.2)),
         ks5=abs(rnorm(1,0.5,0.2)),ks6=abs(rnorm(1,0.5,0.2)),
         ks7=abs(rnorm(1,0.5,0.2)),ks8=abs(rnorm(1,0.5,0.2)),
         k1=abs(rnorm(1,0.5,0.2)),k2=abs(rnorm(1,0.5,0.2)),
         k3=abs(rnorm(1,0.5,0.2)),yc1=abs(rnorm(1,0.5,0.2)),
         yc2=abs(rnorm(1,0.5,0.2)),yc3=abs(rnorm(1,0.5,0.2)),
         yc4=abs(rnorm(1,0.5,0.2)),yc5=abs(rnorm(1,0.5,0.2)),
         yc6=abs(rnorm(1,0.5,0.2)),yc7=abs(rnorm(1,0.5,0.2)),
         yc8=abs(rnorm(1,0.5,0.2)),yc9=abs(rnorm(1,0.5,0.2)),
         yc10=abs(rnorm(1,0.5,0.2)),yc11=abs(rnorm(1,0.5,0.2)),
         yc12=abs(rnorm(1,0.5,0.2)),yc13=abs(rnorm(1,0.5,0.2)),
         yc14=abs(rnorm(1,0.5,0.2)),yc15=abs(rnorm(1,0.5,0.2)),
         yc16=abs(rnorm(1,0.5,0.2)),yc17=abs(rnorm(1,0.5,0.2)),
         yc18=abs(rnorm(1,0.5,0.2)),yc19=abs(rnorm(1,0.5,0.2)),
         yc20=abs(rnorm(1,0.5,0.2)),yc21=abs(rnorm(1,0.5,0.2)),
         ev1=abs(rnorm(1,0.5,0.2)),ev2=abs(rnorm(1,0.5,0.2)),
         ev3=abs(rnorm(1,0.5,0.2)),sigma=abs(rnorm(1,0.5,0.2))
    )
  }
}

#Stan call:

fit = stan(paste0("mm",args[2],".stan"),
           data=c("x0", "t0", "ts", "x", 'T','scl'),
           control=list(adapt_delta=as.numeric(args[3]),
                        stepsize=as.numeric(args[4]),
                        max_treedepth=as.numeric(args[5])),
           warmup = round(as.numeric(args[7])*1/3,0),
           init = ini,
           refresh = 5,
           cores = min(4, parallel::detectCores()),
           chains=4, iter=as.numeric(args[7]), seed=as.numeric(args[6]))

save(fit, file=paste0(lbl,"mm",args[2],".Rsave"))# saving output

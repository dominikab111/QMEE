#assignment 8: Bayesian models
# Setup ------------------------------------------------------------------------

#install required packages
#install.packages("rjags")
#install.packages("R2jags")
#install.packages("coda")
#install.packages("broom.mixed")
#install.packages("arm")
#install.packages("emdbooks")
#install.packages("dotwhisker")
#install.packages("effects")

#load your packages
#setwd('/Users/DOMO/Desktop')
library(vegan)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(lme4)
library(glmmTMB)
library(emmeans)
library(RColorBrewer)
library(rjags)
library(R2jags)
library(coda)
library(broom.mixed)
library(emdbook)
library(dotwhisker)
library(lattice)
library(arm)
library(effects)
namedList <- lme4:::namedList  ## utility

#import csv file
dat <- read.csv("calibration_data.csv")
head(dat)
str(dat)

## create a list

set.seed(20)
N <- nrow(dat)
y <- dat$protein
x <- dat$A750

#bayesian model 
jags1 <- jags(model.file='calibration_data.bug',
                 parameters=c("mx", "int"),
               data = list(N=nrow(dat), y=y, x=x),
                                n.chains = 3,
                                inits=NULL)

jags1 

dat1 <- as.mcmc.bugs(jags1$BUGSoutput)  ## extract the "BUGS output" component
plot(dat1)             ## large-format graph
## plot(mm)                ## trace + density plots, same as above
xyplot(dat1,layout=c(2,3))  ## prettier trace plot
densityplot(dat1,layout=c(2,3)) ## prettier density plot
print(dwplot(jags1) + geom_vline(xintercept=0, lty=2))              ## estimate + credible interval plot

#fit a frequentist linear model
summary(lm(A750~protein, data=dat)) #p value = 1.899e-11

##assignment discussion
# chose weak priors (0.0001) because there is a lot of variation in wavelength absorption 
# to protein concentration in the literature, and it has 3 decimal places (i.e. 0.123)
# citation: Siu Sylvia Lee, Britt Glaunsinger, Fiamma Mantovani, Lawrence Banks, Ronald T. Javier, 
# Multi-PDZ Domain Protein MUPP1 Is a Cellular Target for both Adenovirus E4-ORF1 and High-Risk Papillomavirus Type 18 E6 Oncoproteins, 
# Journal of Virology, 10.1128/JVI.74.20.9680-9693.2000, 74, 20, (9680-9693), (2000).

#the frequentist fit linear model states under the assumption that the true value of the y-intercept 
#is zero in the presence of x, random sampling of the same number of pairs
#would result in a least squares best fit line with a y-int of 1.9e-11.

#the bayesian model & inference plots
#rhat values are close to 1 and n.eff values are 3000 meaning our chain numbers are doing well
#the deviance's rhat value is 1.002 and n.eff value is 3000 meaning it has a good fit for an individual set of samples
#trace plots: looks like a caterpillar/noise look similar to one another therefore it is good
#density: normal distribution, density is sensible and looks fine
#dot whisker plot: intercept is below 0 and is not positive, mx is well above 10 and is very clearly positive. 
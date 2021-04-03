#assignment 8: Bayesian models

# Setup ------------------------------------------------------------------------

## BMB: thanks for commenting all of this out! If you want you can
## vectorize it, e.g.
##  pkgs <- c("R2jags","coda","broom.mixed","arm","emdbook")
##  install.packages(pkgs)
## (note, "emdbook" is misspelled below, and I don't think you need "rjags")
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
## redundant with tidyverse ...
## library(dplyr)
## library(ggplot2)
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

## BMB: do we need all of these packages for this exercise? Try to
## use only necessary packages ...

namedList <- lme4:::namedList  ## utility

#import csv file
dat <- read.csv("calibration_data.csv")
head(dat)
str(dat)

## create a list

set.seed(20)
N <- nrow(dat)
## BMB: I think these are backward! (Your lm predicts A750 from protein,
##  which is I think what you want? In any case they should be consistent.)
## And, you don't have to do this as a separate step (see below)
y <- dat$protein
x <- dat$A750

#bayesian model 
jags1 <- jags(model.file='calibration_data.bug',
              parameters=c("mx", "int"),
              data = list(N=nrow(dat), y=dat$A750, x=dat$protein),
                                n.chains = 3,
                                inits=NULL)

jags1 

dat1 <- as.mcmc.bugs(jags1$BUGSoutput)  ## extract the "BUGS output" component
dat1[] <- lapply(dat1, function(x) mcmc(x[-(1:10),]))  ## drop first few rows - not sure why we need this?
plot(dat1)             ## large-format graph
## plot(mm)                ## trace + density plots, same as above
xyplot(dat1,layout=c(2,3))  ## prettier trace plot
densityplot(dat1,layout=c(2,3)) ## prettier density plot
print(dwplot(jags1) + geom_vline(xintercept=0, lty=2))              ## estimate + credible interval plot



#fit a frequentist linear model
lm1 <- lm(A750~protein, data=dat)
summary(lm1) #p value = 1.899e-11

tt2 <- tidy(lm1,conf.int=TRUE)
tt1 <- tidy(jags1,conf.int=TRUE) %>% mutate(term=tt2$term)

print(dwplot(
    bind_rows(jags=tt1, lm=tt2, .id="model")) +
    geom_vline(xintercept=0, lty=2))              ## estimate + credible interval plot

## BMB: I'm not quite sure why JAGS CI is wider -- should be almost identical in this
## case?

##assignment discussion
# chose weak priors (0.0001) because there is a lot of variation in wavelength absorption
# to protein concentration in the literature, and it has 3 decimal places (i.e. 0.123)

# citation: Siu Sylvia Lee, Britt Glaunsinger, Fiamma Mantovani, Lawrence Banks, Ronald T. Javier, 
# Multi-PDZ Domain Protein MUPP1 Is a Cellular Target for both Adenovirus E4-ORF1 and High-Risk Papillomavirus Type 18 E6 Oncoproteins, 
# Journal of Virology, 10.1128/JVI.74.20.9680-9693.2000, 74, 20, (9680-9693), (2000).

## BMB: that seems reasonable. Not sure what three decimal places has to do with it though

#the frequentist fit linear model states under the assumption that the true value of the y-intercept 
#is zero in the presence of x, random sampling of the same number of pairs
#would result in a least squares best fit line with a y-int of 1.9e-11.

## BMB ??? that doesn't make any sense to me.
## under the null hypothesis (slope==0) the y-intercept would be equal to the mean
## of the data.
## you're quoting the p-value here ???

#the bayesian model & inference plots
#rhat values are close to 1 and n.eff values are 3000 meaning our chain numbers are doing well
#the deviance's rhat value is 1.002 and n.eff value is 3000 meaning it has a good fit for an individual set of samples
#trace plots: looks like a caterpillar/noise look similar to one another therefore it is good
#density: normal distribution, density is sensible and looks fine
#dot whisker plot: intercept is below 0 and is not positive, mx is well above 10 and is very clearly positive. 

## BMB: grade 2/3

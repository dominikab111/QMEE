#Assignment 6 - linear models

# Setup ------------------------------------------------------------------------

#load your packages
library(vegan)
library(tidyverse)
## BMB: tidyverse includes both dplyr and ggplot2, avoid redundant code ...
## library(dplyr)
## library(ggplot2)
library(phyloseq)
library(microbiome)


#if(!requireNamespace("BiocManager", quietly = TRUE)){
#install.packages("BiocManager")
#} 
#BiocManager::install("phyloseq")
#install.packages(c("devtools", "RcppEigen", "RcppParallel", "Rtsne", "ggforce", "units"))
## install.packages("BiocManager"); BiocManager::install("microbiome")

#setwd('/Users/DOMO/Documents/McMaster_University/Surette_lab/Weston_analysis')
#import csv files ----
asvtab = read.csv("asvtab.csv") #asv frequency table
taxtab = read.csv("taxatab.csv") #taxonomy file
mapfile= read.csv("updated_samplemap.csv") #sample information file

#Format tables ----

#convert the rownames of your sample file to DBIDs. 
#These are the fastq IDs so they are essential to use here

rownames(mapfile) = mapfile$WPNum
rownames(mapfile)
colnames(mapfile) 

#convert the asv sequences to the rownames for both the asvtab and taxtab
#this is important because we need to have matching rownames to convert these to 
#phyloseq objects. 

rownames(asvtab) = asvtab$ASV
rownames(taxtab) = taxtab$X
head(taxtab)

#checks if the rownames in samdat match the colnames of asvtab. prints a logical response
#next line will print out the sample IDs that do not match

## BMB: it's good to put stopifnot() around these: it will automatically
## stop your script if the conditions you want to be true aren't ...
all(rownames(mapfile) %in% colnames(asvtab))
length(rownames(mapfile)[!(rownames(mapfile) %in% colnames(asvtab))])==0

## BMB: didn't you print out the rownames of the mapfile just a few lines ago?
rownames(mapfile)
colnames(asvtab) 

#checks if the colnames in asvtab match the rownames of the sample file. logical response
#next line will print out the sample IDs that do not match

all(colnames(asvtab) %in% rownames(mapfile))
colnames(asvtab)[!(colnames(asvtab) %in% rownames(mapfile))]
## BMB: should this be empty?

#Convert to matrices -----------------------------------------------------------
#To use this data in a phyloseq object we must convert our dataframes 
#to matrices. The asvtab must be an integer, and the taxtab must be character. 

#converting data frame to integer matrix, use dim(asvtab)
dim(asvtab)
asv <- asvtab[,2:274] 
asvtable <- data.matrix(asv,rownames.force = NA)
class(asvtable)
str(asvtable) #make sure asvtab is integers and characters, checks the structure of the data

#converting dataframe into character matrix
dim(taxtab) #to determine how to truncate
tax <- taxtab[,-1] #remove first column

taxtable <- as.matrix(tax,rownames.force = NA)
class(taxtable) #make sure it is a matrix at this point
stopifnot(inherits(taxtable,"matrix"))
str(taxtable) #make sure output is all characters

#Making the phyloseq object ----------------------------------------------------

dat = phyloseq(otu_table(asvtable, taxa_are_rows = TRUE),
               tax_table(taxtable),
               sample_data(mapfile))
dat

#check sample sums
summary(sample_sums(dat))
sample_sums(dat)


sample_data(dat)[sample_sums(dat) < 2500,]

#filter out mitochondria bacteria from host
dat = subset_taxa(dat, Kingdom=="Bacteria", Family!="Mitochondria")

samdat_clean = data.frame(sample_data(dat)) #make a dataframe

#transform asv counts into relative abundance data (i.e. calculate relative abundance) 

dat_rel = transform_sample_counts(dat, function(x) x/sum(x)) 

#calculate alpha diversity
alpha_div <- microbiome::alpha(dat_rel, index = "diversity_shannon")
permdat <- merge(mapfile, alpha_div, by="row.names") #merging the alpha diversity calculations to the sample mapfile

## linear regression model: predicting diversity based on age
lmage <- lm(diversity_shannon~Age, data= permdat)
par(mfrow=c(2,2)) ##views al 4 plots at once
summary(lmage) ##p value = 0.1372 
plot(lmage, id.n=5) #approx 5 outliers passed Cook's distance line

## linear regression model: predicting diversity based on sample type
lmst <- lm(diversity_shannon~Sample.Type, data= permdat)
par(mfrow=c(2,2))
summary(lmst) ## p value = 0.01502
plot(lmst, id.n=5)

## BMB: any inferential plots?



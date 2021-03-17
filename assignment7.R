#Assignment 7 - generalized linear models 

# Setup ------------------------------------------------------------------------

#load your packages
library(vegan)
library(tidyverse)
library(phyloseq)
library(microbiome)
library(dplyr)
#if(!requireNamespace("BiocManager", quietly = TRUE)){
#install.packages("BiocManager")
#} 
#BiocManager::install("phyloseq")
#install.packages(c("devtools", "RcppEigen", "RcppParallel", "Rtsne", "ggforce", "units"))
library(ggplot2)
## install.packages("BiocManager"); BiocManager::install("microbiome")

setwd('/Users/DOMO/Documents/McMaster_University/Surette_lab/Weston_analysis')
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

all(rownames(mapfile) %in% colnames(asvtab))
rownames(mapfile)[!(rownames(mapfile) %in% colnames(asvtab))]

rownames(mapfile)
colnames(asvtab) 

#checks if the colnames in asvtab match the rownames of the sample file. logical response
#next line will print out the sample IDs that do not match

all(colnames(asvtab) %in% rownames(mapfile))
colnames(asvtab)[!(colnames(asvtab) %in% rownames(mapfile))]

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

#remove samples with less than 2500 reads
dat = prune_samples(sample_sums(dat)>2500, dat)

#filter out mitochondria bacteria from host
dat = subset_taxa(dat, Kingdom=="Bacteria", Family!="Mitochondria")

#rarefy the data
dat_rare = rarefy_even_depth(dat) 

#transform asv counts into relative abundance data (i.e. calculate relative abundance) 
dat_rel = transform_sample_counts(dat, function(x) x/sum(x)) 

-------------------------------------------------------------------------------
  ##Assignment 6
-------------------------------------------------------------------------------
  
#calculate alpha diversity
alpha_div <- microbiome::alpha(dat_rel, index = "diversity_shannon")
permdat <- merge(mapfile, alpha_div, by="row.names") #merging the alpha diversity calculations to the sample mapfile

## using GLM to predict diversity based on age
glm_age <- glm(diversity_shannon~Age, family=Gamma(link="log"), data=permdat)
par(mfrow=c(2,2)) #allows viewing of all 4 plots at once
summary(glm_age) #p value = 0.217
plot(glm_age)

#GLM predicting diversity based on sample type
glm_st <- glm(diversity_shannon~Sample.Type, family=Gamma(link="log"), data=permdat)
par(mfrow=c(2,2))
summary(glm_st)
plot(glm_st)

#comparing to LM
lm_age1 <- lm(log(diversity_shannon)~Age, data=permdat)
par(mfrow=c(2,2)) 
summary(lm_age1) # p value = 0.104
plot(lm_age1)

## GLM predicting diversity based on sample type
lm_st <- lm(log(diversity_shannon)~Sample.Type, data=permdat)
par(mfrow=c(2,2))
summary(lm_st) ## p value = 0.00257
plot(lm_st, id.n=5)

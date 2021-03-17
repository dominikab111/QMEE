#Assignment 6 - generalized linear models 

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


## BMB: please DON'T use setwd() (if it has to be in here, comment it out)
## setwd('/Users/DOMO/Documents/McMaster_University/Surette_lab/Weston_analysis')
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

#remove samples with less than 2500 reads
dat = prune_samples(sample_sums(dat)>2500, dat)

#filter out mitochondria bacteria from host
dat = subset_taxa(dat, Kingdom=="Bacteria", Family!="Mitochondria")

#rarefy the data
dat_rare = rarefy_even_depth(dat)
## BMB: are you sure want to do this? This is a contentious step

#transform asv counts into relative abundance data (i.e. calculate relative abundance) 
dat_rel = transform_sample_counts(dat, function(x) x/sum(x)) 

-------------------------------------------------------------------------------
  ##Assignment 6
-------------------------------------------------------------------------------
  
#calculate alpha diversity
alpha_div <- microbiome::alpha(dat_rel, index = "diversity_shannon")
permdat <- merge(mapfile, alpha_div, by="row.names") #merging the alpha diversity calculations to the sample mapfile

## BMB: this is assignment 6, not 7 ... any particular reason
## you're fitting GLMs? Are you combining these two assignments?

## using GLM to predict diversity based on age
glm_age <- glm(diversity_shannon~Age, family=Gamma(link="log"), data=permdat)
par(mfrow=c(2,2)) #allows viewing of all 4 plots at once
summary(glm_age) #p value = 0.217
plot(glm_age)
## BMB: any conclusions about this plot? (Looks pretty good to me;
## slight trend toward lower variability with increased Shannon diversity
## (also, seems like a pretty narrow range of diversities: 0.78-0.88 ?

#GLM predicting diversity based on sample type
glm_st <- glm(diversity_shannon~Sample.Type, family=Gamma(link="log"), data=permdat)
par(mfrow=c(2,2))
summary(glm_st)
plot(glm_st)
## BMB: pretty boring (that's not a bad thing), as is typical of models with only one two-level categorical predictor

#comparing to LM
lm_age1 <- lm(log(diversity_shannon)~Age, data=permdat)
par(mfrow=c(2,2)) 
summary(lm_age1) # p value = 0.104
plot(lm_age1)
## BMB: the Q-Q plot at least looks worse .. 

## GLM predicting diversity based on sample type
lm_st <- lm(log(diversity_shannon)~Sample.Type, data=permdat)
par(mfrow=c(2,2))
summary(lm_st) ## p value = 0.00257
plot(lm_st, id.n=5)

## BMB: any inferential plots?


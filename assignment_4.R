# Assignment 4 - computing permutation tests

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


sample_data(dat)[sample_sums(dat) < 2500,]

#filter out mitochondria bacteria from host
dat = subset_taxa(dat, Kingdom=="Bacteria", Family!="Mitochondria")

samdat_clean = data.frame(sample_data(dat)) #make a dataframe

#transform asv counts into relative abundance data (i.e. calculate relative abundance) 

dat_rel = transform_sample_counts(dat, function(x) x/sum(x)) 

#subset data into nasal and oral samples
dat.nasal = subset_samples(dat_rel, Sample.Type=='Nasal')
dat.sal = subset_samples(dat_rel, Sample.Type=='Saliva')

# testing compositional differences between sample types

bray_d <- phyloseq::distance(dat_rel, method = "bray") #calculate bray curtis dissimilarity

vegan::adonis(bray_d ~ phyloseq::sample_data(dat_rel)$Sample.Type)
dispr <- vegan::betadisper(bray_d, phyloseq::sample_data(dat_rel)$Sample.Type)
dispr
plot(dispr, main = "Ordination Centroids and Dispersion Labelled: Bray Curtis Dissimilarity Sample Types")

permutest(dispr) 

# testing compositional differences between sexes (m/f)

vegan::adonis(bray_d ~ phyloseq::sample_data(dat_rel)$Sex)
dispr_s <- vegan::betadisper(bray_d, phyloseq::sample_data(dat_rel)$Sex)
dispr_s
plot(dispr_s, main = "Ordination Centroids and Dispersion Labelled: Bray Curtis Dissimilarity Sexes")
permutest(dispr_s)

#calculating alpha diversity

alpha_div <- microbiome::alpha(dat_rel, index = "diversity_shannon")
print(alpha_div)

as.data.frame(mapfile)

#brute force permutation test
set.seed(101) ##for reproducibility
nsim <- 9999
res <- numeric(nsim) ##set aside space for results
for (i in 1:nsim) {
  ## standard approach: scramble response value
  perm <- sample(nrow(mapfile))
  bdat <- transform(mapfile,Age=Age[perm])
  ## compute & store difference in means; store the value
  res[i] <- mean(bdat$mapfile[bdat$Time=="Baseline"])-
    mean(bdat$mapfile[bdat$Time=="FollowUp"])
}
obs <- mean(mapfile$Age[mapfile$Time=="Baseline"])-
  mean(mapfile$Age[mapfile$Time=="FollowUp"])
## append the observed value to the list of results
res <- c(res,obs)

hist(res,col="grey",las=1, main="")
abline(v=obs,col="red")

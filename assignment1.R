#analysis for Weston-Boris 2016-2019 Data
# Goals - to perform a substantive calculation on my data set.

# Setup -------------------------------------

#set your working directory for your computer
## JD: No. Don't put this in your script.
## It turns out I don't have a directory called Users, so I can't make a directory called '/Users/DOMO/Documents/McMaster_University/Surette_lab/Weston_analysis' even if I wanted to.
## setwd('/Users/DOMO/Documents/McMaster_University/Surette_lab/Weston_analysis')

## Be selective about which libraries you call; if you are trying to be reproducible and careful

## BMB: I would repeat JD's comment here; do you really need all of this stuff?
library(vegan)
library(tidyverse)
library(phyloseq)

## BMB: I can run
if (FALSE) {
library(ggplot2)
library(lme4)
library(glmmTMB)
library(emmeans)
library(boot)
}

## BMB: may be worth noting that this is a Bioconductor package
## install.packages("BiocManager"); BiocManager::install("microbiome")
## (that depends on who you are collaborating with and whether
##  they are familiar with Bioconductor machinery)
library(microbiome)
library(knitr)

## JD: I have no access to these files
#import csv files ----
seqtab = read.csv("asvtab.csv") #asv frequency table
taxatab = read.csv("taxatab.csv") #taxonomy file
mapfile= read.csv("updated_samplemap.csv") #sample information file

#rename variables 
samdat_clean <- mapfile
taxtab <- taxatab
asvtab <- seqtab

#Format tables ----

#check to see if your tables have the same orientation
head(asvtab)
head(taxtab)
dim(asvtab)
dim(taxtab)

# formatting your asv and taxa tables

#convert the rownames of your sample file to DBIDs. 
#These are the fastq IDs so they are essential to use here

rownames(samdat_clean) = samdat_clean$WPNum
rownames(samdat_clean)
colnames(samdat_clean) 

#convert the asv sequences to the rownames for both the asvtab and taxtab
#this is important because we need to have matching rownames to convert these to 
#phyloseq objects. 

rownames(asvtab) = asvtab$ASV
rownames(taxtab) = taxtab$X
head(taxtab)

#head is a way to look at top 6 rows of table, allows you to make sure the data looks the way you think it should
#next we need to check that the sample IDs present in the sample file 
#match the fastq_ID's

#checks if the rownames in samdat match the colnames of asvtab. prints a logical response
#next line will print out the sample IDs that do not match

## BMB: you could add stopifnot() here (this will automatically stop
##  your script if some condition is not as you expected)
all(rownames(samdat_clean) %in% colnames(asvtab))
rownames(samdat_clean)[!(rownames(samdat_clean) %in% colnames(asvtab))]

rownames(samdat_clean)
colnames(asvtab)

#checks if the colnames in asvtab match the rownames of the sample file. logical response
#next line will print out the sample IDs that do not match

all(colnames(asvtab) %in% rownames(samdat_clean))
colnames(asvtab)[!(colnames(asvtab) %in% rownames(samdat_clean))]


#Convert to matrices -----
#To use this data in a phyloseq object we must convert our dataframes 
#to matrices. The asvtab must be an integer, and the taxtab must be character. 

#converting data frame to integer matrix, use dim(asvtab)
dim(asvtab)
asv <- asvtab[,2:274] #numbers change based on your dataset. use dim(asvtab) to determine
asvtable <- data.matrix(asv,rownames.force = NA)
class(asvtable)
str(asvtable) #make sure asvtab is integers and characters, checks the structure of the data

#converting dataframe into character matrix
dim(taxtab) #to determine how to truncate
tax <- taxtab[,2:7] #numbers are determined using the above function  
taxtable <- as.matrix(tax,rownames.force = NA)
class(taxtable) #make sure it is a matrix at this point 
str(taxtable) #make sure output is all characters

#Making the phyloseq object ----

dat = phyloseq(otu_table(asvtable, taxa_are_rows = TRUE),
               tax_table(taxtable),
               sample_data(samdat_clean))
dat #make sure the numbers match the number of samples you think you should have

#check sample sums
summary(sample_sums(dat))
sample_sums(dat)

#host bacteria is removed after this step
dat = subset_taxa(dat, Kingdom=="Bacteria", Family!="Mitochondria")

#rarefy the data by removing less than 2500 reads
dat = prune_samples(sample_sums(dat)>2500, dat)

#check to see how many reads were removed
dat #in this case it was 17

samdat_clean = data.frame(sample_data(dat)) #make a dataframe

## BMB: you probably don't need phyloseq machinery if all you want
## to do is filter ...

#transform asv counts into relative abundance data (i.e. calculate relative abundance) 

dat_rel = transform_sample_counts(dat, function(x) x/sum(x)) 

# Calculating Alpha diversity -----------------------------------------------------------

alpha_diversity1 <- microbiome::alpha(dat_rel, index = "diversity_shannon")
alpha_diversity2 <- microbiome::alpha(dat_rel, index = "diversity_inverse_simpson")
alpha_diversityrel <- microbiome::alpha(dat_rel, index = "dominance_relative")
alpha_diversityabs <- microbiome::alpha(dat_rel, index = "dominance_absolute")

## BMB: this line didn't work
## print(alpha_diversitycalc)
print(alpha_diversity2)
print(alpha_diversityrel)
print(alpha_diversityabs)

# Calculating Beta- Diversity ----------------------------------------------------------

ns_rel_otu <- data.frame(phyloseq::otu_table(dat_rel))
ns_rel_otu <- t(ns_rel_otu)
dist.mat <- vegan::vegdist(ns_rel_otu, method = "bray")
## BMB: it probably doesn't make sense to print this giant object out ...
##  what are you going to see?
##  have you considered heatmap(as.matrix(dist.mat)) ... ?
print(dist.mat)

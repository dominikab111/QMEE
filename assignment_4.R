#Assignment 4 - computing permutation tests

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

source("cleandata.R")
dat_rel <- readRDS("dat_rel.rds")

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

#brute force permutation test
set.seed(101)
nsim <- 9999
res <- numeric(nsim)
for (i in 1:nsim) {
  perm <- sample(nrow(alpha_div))
  bdat <- transform(alpha_div,diversity_shannon=diversity_shannon[perm])
  res[i] <- mean(bdat$alpha_div[bdat$diversity_shannon])}

# Assignment 3 - producing beta diveristy and relative abundance plots

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

## plotting a histogram to explore age ranges in the data 
hist(mapfile$Age,
     main= "Histogram of Participants Age",
     xlab="Age",
     border= "navy blue",
     col="maroon",
     breaks=10,
     prob=TRUE)

lines(density(mapfile$Age))

## Plotting alpha diversity 

#Plot richness for nasal samples
plot_richness(dat.nasal,measures=c("Shannon", "Simpson"))  +
  theme(panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

#Plot richness for saliva samples
plot_richness(dat.sal, measures = c("Shannon", "Simpson")) +
  theme(panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

#different way of plotting alpha diversity (arguably more efficient)

types <- unique(sample_data(dat_rel)$Sample.Type)
for (tt in types) {
  s <- subset_samples(dat, Sample.Type==tt)
  print(plot_richness(s,measures = c("Shannon", "Simpson")))
}

#the alpha diversity richness plots will reproduce if you run the line individually, otherwise a warning message is produced.
#I spoke to a lab member and the singletons have been removed with DADA2 processing but 
# they recommended dismissing the warning message and running the code one line at a time.

# Plotting beta diversity

ordu <- ordinate(dat.nasal, "PCoA", "bray", weighted = TRUE)
plot_ordination(dat.nasal, ordu, color="Age", shape="Time",
                title = "Bray Curtis Dissimilartiy of Nasal samples against Time")
plot_ordination(dat.nasal, ordu, color="Age", shape="Sex",
                title = "Bray Curtis Dissimilarity of Nasal samples against Sex")

# Plotting relative abundance

#to get the top 20 genera in the samples 
top20 <- names(sort(taxa_sums(dat_rel), decreasing=TRUE))[1:20]
dat.top20 <- transform_sample_counts(dat_rel, function(OTU) OTU/sum(OTU))
dat.top20 <- prune_taxa(top20, dat.top20)

# plotting the taxa bar plots
phyloseq::plot_bar(dat.top20, fill = "Genus") +
  geom_bar(aes(color = Genus, fill = Genus), stat = "identity", position = "stack", width = 1) +
  labs(x = "", y = "Relative Abundance\n") +
  theme(panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

phyloseq::plot_bar(dat.top20, fill = "Genus") +
  geom_bar(aes(color = Genus, fill = Genus), stat = "identity", position = "stack", width = 1) +
  labs(x = "", y = "Relative Abundance\n") +
  facet_wrap(~ Sample.Type, scales = "free")
theme(panel.background = element_blank(),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank())
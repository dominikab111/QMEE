# Assignment 3 objective- producing beta diveristy and relative abundance plots

# Setup ------------------------------------------------------------------------

#setwd('~/Documents/McMaster_University/QMME/qmee)


#load your packages
library(vegan)
library(tidyverse)
library(phyloseq)
library(microbiome)
library(dplyr)
library(ggplot2)
#if(!requireNamespace("BiocManager", quietly = TRUE)){
#install.packages("BiocManager")
#} 
#BiocManager::install("phyloseq")
#install.packages(c("devtools", "RcppEigen", "RcppParallel", "Rtsne", "ggforce", "units"))
## install.packages("BiocManager"); BiocManager::install("microbiome")

samdat_clean <- readRDS("samdat_clean.RData")

#transform asv counts into relative abundance data (i.e. calculate relative abundance) 

## dat_rel = transform_sample_counts(datp, function(x) x/sum(x))
dat_rel <- readRDS("dat_rel.rds")

#subset data into nasal and oral samples
dat.nasal = subset_samples(dat_rel, Sample.Type=='Nasal')
dat.sal = subset_samples(dat_rel, Sample.Type=='Saliva')

## plotting a histogram to explore age ranges in the data 
hist(mapfile$Age,
     main= "Histogram of Participants Age",
     xlab="Age",
     border= "navy blue",
     col="maroon",
     breaks=10,
     prob=TRUE)

lines(density(mapfile$Age)) # add lines to observe density

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

#these plots will reproduce if you run the line individually, otherwise a warning message is produce.
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
  facet_wrap(~ Sample.Type, scales = "free") +
  theme(panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

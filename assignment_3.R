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
  
#subset data into nasal and oral samples
dat.nasal = subset_samples(dat_rel, Sample.Type=='Nasal')
dat.sal = subset_samples(dat_rel, Sample.Type=='Saliva')
  
#-------------------------------------------------------------------------------
## Assignment 3
# -------------------------------------------------------------------------------
  
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
  
  #different way of plotting alpha diversity (potentially more efficient)
  
  types <- unique(sample_data(dat_rel)$Sample.Type)
  for (tt in types) {
    s <- subset_samples(dat, Sample.Type==tt)
    print(plot_richness(s,measures = c("Shannon", "Simpson")))
  }
  
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

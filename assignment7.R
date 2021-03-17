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

#only working with nasal samples
dat_nas = subset_samples(dat_rare, Sample.Type=='Nasal')

#transform asv counts into relative abundance data (i.e. calculate relative abundance) 

dat_rel = transform_sample_counts(dat_nas, function(x) x/sum(x)) 

#-------------------------------------------------------------------------------
  ##Assignment 7
#-------------------------------------------------------------------------------
  
#calculate alpha diversity
alpha_div <- microbiome::alpha(dat_rel, index = "diversity_shannon")
permdat <- merge(mapfile, alpha_div, by="row.names") #merging the alpha diversity calculations to the sample mapfile

# hypothesis: alpha diversity will be negatively associated with age
glm_age <- glm(diversity_shannon~Age, family=Gamma(link="log"), data=permdat)
par(mfrow=c(2,2)) #allows viewing of all 4 plots at once
plot(glm_age)
summary(glm_age) #p value = 0.217
#plot discussion
#residuals vs fitted: the red line is relatively flat, therefore we can infer there are no
#non-linear relationships present
#normal q-q: the data points remain diagonally straight therefore there is little data that is not 
#normally distributed
#scale-location: the red line seems relatively straight which indicates there is little heteroscedasticity within this model
#residuals vs leverage: I cannot see the Cook's contours/distance lines so this is not a concern - I am assuming everything is okay here 
#and there are no data points that should be treated as outliers.

#comparing to linear model
lm_age <- lm(log(diversity_shannon)~Age, data=permdat)
par(mfrow=c(2,2)) 
plot(lm_age)
summary(lm_age) # p value = 0.104
#plot discussion
#residuals vs fitted: the red line is almost perfectly flat, indicating there are no non-linear relationships
# in the model with alpha diversity and age
#normal q-q: the data points form an upside down 'U' shape and do not entirely follow the diagonal line. There may
# be some data that are not normally distributed, but are mostly okay.
#scale - location: the red line is sloped downwards which could indicate some heteroscedasticity (change in variance) within this dataset/model,
#and may be concerning
#residuals vs leverage: I can just vaguely see the Cook's contour/distance line as '0.5' in the bottom right hand corner - but no points 
#lie past this line so again, I will not treat any data points as outliers. 

#Rather than an inferential model, I decided to produce a linear regression model showing the 
#relationship between diversity & age in samples collected from the nasal cavity region

lm_age_div_nas <- ggplot(permdat, aes(x=Age, y=diversity_shannon, ylab = "Alpha diversity")) +
  geom_point(shape=1) + #use hollow circles
  geom_smooth(method=lm) #add a linear regression line
lm_age_div_nas
#alpha diversity seems to slightly decrease with age in nasal swab samples

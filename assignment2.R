#analysis for Weston-Boris 2016-2019 Data
# Goals - to perform a substantive calculation & produce one or two plots

# Setup ------------------------------------------------------------------------

#load your packages
library(vegan)
library(tidyverse)
library(phyloseq)
library(microbiome) ## install.packages("BiocManager"); BiocManager::install("microbiome")

## library(knitr) ## BMB: what is this for??

#import csv files ----
seqtab = read.csv("asvtab.csv") #asv frequency table
taxatab = read.csv("taxatab.csv") #taxonomy file
mapfile= read.csv("updated_samplemap.csv") #sample information file

#rename variables
## BMB: why read these in with one name and immediately rename them?
## why not assign them to the name/symbol you want to use directly?
samdat_clean <- mapfile
taxtab <- taxatab
asvtab <- seqtab

#Format tables ----

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

#checks if the rownames in samdat match the colnames of asvtab. prints a logical response
#next line will print out the sample IDs that do not match

## BMB: stopifnot() is good to use here: it will stop your
## script if the condition is FALSE
all(rownames(samdat_clean) %in% colnames(asvtab))
rownames(samdat_clean)[!(rownames(samdat_clean) %in% colnames(asvtab))]

## BMB: could do something like this:
rs <- rownames(samdat_clean)
bad <- which(!rs %in% colnames(asvtab))
if (length(bad)>0) {
    stop("non-matching rownames: ",
         paste(rs[bad], collapse=", "))
}
## BMB: and similarly for the next condition
## setdiff() could also be useful (look up ?setdiff and experiment)

rownames(samdat_clean)
colnames(asvtab) 
#stopifnot() #automatically stops script if condition is not as expected

#checks if the colnames in asvtab match the rownames of the sample file. logical response
#next line will print out the sample IDs that do not match

all(colnames(asvtab) %in% rownames(samdat_clean))
colnames(asvtab)[!(colnames(asvtab) %in% rownames(samdat_clean))]


#Convert to matrices -----------------------------------------------------------
#To use this data in a phyloseq object we must convert our dataframes 
#to matrices. The asvtab must be an integer, and the taxtab must be character. 

#converting data frame to integer matrix, use dim(asvtab)
dim(asvtab)
asv <- asvtab[,2:274]  ## BMB: avoid hard-coded numbers
## why 2:274?  would janitor::remove_empty() be useful?
asvtable <- data.matrix(asv,rownames.force = NA)
class(asvtable)
str(asvtable) #make sure asvtab is integers and characters, checks the structure of the data
## BMB: again, it would be good to test and AUTOMATICALLY
## stop if something looks wrong

#converting dataframe into character matrix
dim(taxtab) #to determine how to truncate
tax <- taxtab[,2:7] #numbers are determined using the above function
## BMB how are these determined?  do you always want to leave out the first
## column? tax <- taxtab[,-1] is probably best
## (or taxtab[,2:ncol(taxtab)] but using -1 is shorter/clearer
taxtable <- as.matrix(tax,rownames.force = NA)
class(taxtable) #make sure it is a matrix at this point
## BMB: stopifnot(inherits(taxtable,"matrix"))
## although since you _just_ explicitly converted it to a matrix
## it might be redundant to check that ...
str(taxtable) #make sure output is all characters

#Making the phyloseq object ----------------------------------------------------

dat = phyloseq(otu_table(asvtable, taxa_are_rows = TRUE),
               tax_table(taxtable),
               sample_data(samdat_clean))
dat

#check sample sums
summary(sample_sums(dat))
sample_sums(dat)

#filter out mitochondria bacteria from host
dat = subset_taxa(dat, Kingdom=="Bacteria", Family!="Mitochondria")
## BMB: I suspect but am not sure that this doesn't work the way you
## think it does?  dplyr::filter() allows multiple criteria (to be
## combined with AND, i.e. all criteria must be TRUE to retain an
## observation), but base::subset() doesn't.  I would suggest
## Kingdom=="Bacteria" & Family!="Mitochondria".  (I would also
## suggest making a very small, simple example and experimenting with
## it to see if I'm right or not)

#rarefy the data by removing less than 2500 reads
## BMB: changed this to 'datp' (pruned data) so I would be able
## to go back and work with unpruned data
datp = prune_samples(sample_sums(dat)>2500, dat)

## http://evomics.org/wp-content/uploads/2016/01/phyloseq-Lab-01-Answers.html#taxa-total-counts-histogram
tt <- tibble(tax_table(dat),
           TotalCounts = taxa_sums(dat),
           OTU = taxa_names(dat))
hist(log10(1+tt$TotalCounts))
tt %>% filter(TotalCounts<20) %>% count(TotalCounts)
tt %>% filter(TotalCounts<50) %>% pull(TotalCounts) %>% table() %>% plot()
## why do we have *no* OTUs with exactly 1 count ...? and we have OTUs
## with 0 counts?

#check to see how many reads were removed
dat #in this case it was 17
## BMB: there ought to be a programmatic way to retrieve this,
## e.g.
ncol(get_sample(datp))

samdat_clean = data.frame(sample_data(datp)) #make a dataframe

#transform asv counts into relative abundance data (i.e. calculate relative abundance) 

dat_rel = transform_sample_counts(datp, function(x) x/sum(x)) 

#subset data into nasal and oral samples
dat.nasal = subset_samples(dat_rel, Sample.Type=='Nasal')
dat.sal = subset_samples(dat_rel, Sample.Type=='Saliva')

## BMB: can you figure out a way to do these computations without
## repeating lots of code?  for() loop or tidyverse?  It doesn't
## look like phyloseq objects work naturally with group_by():
## https://github.com/joey711/phyloseq/issues/1105

types <- unique(sample_data(dat_rel)$Sample.Type)
for (tt in types) {
    s <- subset_samples(dat, Sample.Type==tt)
    print(plot_richness(s,measures = c("Shannon", "Simpson")))
}

# Calculating Alpha diversity -----------------------------------------------------------

alpha_divnasal <- microbiome::alpha(dat.nasal, index = "diversity_shannon")
alpha_divsal <- microbiome::alpha(dat.sal, index = "diversity_shannon")

## BMB: in general, printing out large objects doesn't make much sense
print(alpha_divnasal)
print(alpha_divsal)

# Plotting alpha diversity

#plot richness with shannon and simpson measurements
plot_richness(dat.nasal, measures = c("Shannon", "Simpson"))
plot_richness(dat.sal, measures = c("Shannon", "Simpson"))

samp_vec <- c("Nasal","Saliva")


#I received a warning message when running my plot_richness commands... I have pruned my samples previously
#and tried to prune my taxa using dat_prune <- prune_taxa((dat_nas_rel) > 0, dat_rel) but it did not work
# any suggestions?

#to make up for this I produced a heat map with bray-curtis dissimilarity as well 

# Plotting Beta-Diversity ----------------------------------------------------------

ns_rel_otu <- data.frame(phyloseq::otu_table(dat_rel))
ns_rel_otu <- t(ns_rel_otu)
dist.mat <- vegan::vegdist(ns_rel_otu, method = "bray")
heatmap(as.matrix(dist.mat))

# Plotting relative abundance
phyloseq::plot_bar(dat_rel, fill = "Genus") +
  geom_bar(aes(color = Genus, fill = Genus), stat = "identity", position = "stack") +
  labs(x = "", y = "Relative Abundance\n") +
  facet_wrap(~ Sample.Type, scales = "free") +
  theme(panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

#histogram of age
hist(mapfile$Age,
     main= "Histogram of Participants Age",
     xlab="Age",
     border= "navy blue",
     col="maroon",
     breaks=10,
     prob=TRUE)

lines(density(mapfile$Age))


#Plotting PCoA relative abundance 

relabunbray.ord <- ordinate(rel_abun_lessHOST_prune, method = "PCoA", distance = "bray")
relabunbray.plot <- plot_ordination(rel_abun_lessHOST_prune, relabunbray.ord,
                                    color = "illrun",
                                    axes = c(2,3),
                                    title = "Bray-Curtis relative abundance")
relabunbray.plot+geom_point(size=3)+scale_color_manual(values=colours)
p2.plot <- plot_ordination(rel_abun_lessHOST_prune, relabunbray.ord,
                           color = "illrun",
                           axes = c(1, 2),
                           title = "Bray-Curtis relative abundance")
p2.plot+geom_point(size=3)+scale_color_manual(values=c("#e6194B", "#3cb44b", 
                                                       "#4363d8", "#f58231", "#911eb4",
                                                       "#42d4f4", "#f032e6", "#bfef45",
                                                       "#fabebe", "#469990", "#e6beff",
                                                       "#9A6324", "#fffac8", "#800000",
                                                       "#aaffc3", "#808000", "#ffd8b1",
                                                       "#000075", "#a9a9a9", "#ffffff",
                                                       "#ffe119","#000000"))
rarefiedbray.ord <- ordinate(dat_r, method = "PCoA", distance = "bray")
p3.plot <- plot_ordination(dat_r, rarefiedbray.ord,
                           color = "illrun",
                           axes = c(1,2),
                           title = "Bray-Curtis rarefied 10000")
p3.plot+geom_point(size=3)+scale_color_manual(values=c("#e6194B", "#3cb44b", "#ffe119",
                                                       "#4363d8", "#f58231", "#911eb4",
                                                       "#42d4f4", "#f032e6", "#bfef45",
                                                       "#fabebe", "#469990", "#e6beff",
                                                       "#9A6324", "#fffac8", "#800000",
                                                       "#aaffc3", "#808000", "#ffd8b1",
                                                       "#000075", "#a9a9a9", "#ffffff",
                                                       "#000000"))

## Taxa bar plots --------------------------------------
top20 <- names(sort(taxa_sums(dat_lessExp2_lessHOST), decreasing=TRUE))[1:20]
dat.top20 <- transform_sample_counts(dat_lessExp2_lessHOST, function(OTU) OTU/sum(OTU))
dat.top20 <- prune_taxa(top20, dat.top20)

plot_bar(dat.top20, fill='Genus')


# Plotting alpha diversity -----------------------------

#alpha diversity
plot_richness(dat_rel, x="illrun", measures=c("Shannon", "Simpson"))
#alpha = estimate_richness(dat_lessHOST, measures=c("Shannon", "Simpson"))
#write.csv(alpha, file="alpha_Shannon_Simpson_lessHOST.csv")
alpha1 = plot_richness(dat_r, x="illrun", measures=c("Shannon", "Simpson"), color = "illrun")
alpha1 + scale_color_manual(values = colours2) + geom_point(size=3)


ggplot(tbl, aes(Treatment, Count)) +
  geom_boxplot() +
    theme_bw()+facet_wrap("Group", scales="free")+ggtitle("GFF2 Experiment")

saveRDS(dat_rel, file="dat_rel.rds")

  theme_bw()+facet_wrap("Group", scales="free")+ggtitle("GFF2 Experiment")

## BMB: I don't know what's going on with your missing singletons.
## Not clear whether you should be normalizing counts before doing richness
## plots?

## mark: 2
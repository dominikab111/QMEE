# qmee
# Assignment 1
The data set I have for this weeks assignment includes a taxanomy table, asv table and sample map. It is a data set from a previous cross-sectional study in the Bowdish lab, conducted during 2016-2018 to observe the compositional changes in the anterior nares and oral cavity of the upper respiratory tract (URT) pre- and post-respiratory infections. The sample map includes the ID numbers, age, sex, sample type and time point of donation in 44 individuals. Samples were provided from nasal and saliva swabs, and were collected during Baseline (September to November), Infection (or 'Event', October to February) and Follow Up (March to August). The data was sequenced with 16S rRNA sequencing in the v3 and v4 region of both the nasal and saliva samples. A few biological questions I have for this data set include, "Does the URT microbiome of older individuals re-stabilize after infection?", "How do compositional changes in the URT microbiome contribute to respiratory infections?" and "How does the URT microbiome stability fluctuate pre- and post- infection within young and older individuals? "

# Assignment 2
I would like to investigate age-related changes of the URT microbiome within healthy adults. I would investigate these changes with stability over different time points using alpha and beta diversity measurements. I would also investigate compositional changes within the URT microbiome of healthy adults using taxonomic bar plots. I would divide these investigations into the following components:

1. Alpha diversity measurements between Baseline and Follow Up samples of healthy adults.
2. Beta diversity measurements between 18-39 (young), 40-64(mid aged), and 65+ (older) adults at Baseline.
3. Beta diversity measurements between healthy adults at nasal and oral sample sites.
4. Taxonomic bar plots to investigate compositional changes at Baseline and between Baseline and Follow Up samples.
5. Taxonomic bar plots between Baseline and Infection (or 'Event') samples to assess changes in microbial dominance during respiratory infections.


#Assignment 3

For this assignment I wanted to explore the age categories within my sample set, the diversity of the nasal samples and the composition of the upper respiratory tract (URT) microbiome. 
I decided to explore the age categories with a histogram and plotted a line density to understand the different age distributions within my dataset, and how I should categorized the ages (ex. 18-39, 40-64, 65+) for future plots. I chose maroon and red colours because they were visually appealing.

I chose to calculate and plot Bray-Curtis Dissimilarity to understand the microbial structure and stability within the nasal swab samples. I also chose to plot the metrics to visualize if there is  clustering or an association between various timepoints that the samples were collected (Baseline, Event, FollowUp) or sex (Male, Female). Furthermore, I plotted these categories against age to determine if there was an age association with the different timepoints and sex. I chose a coloured gradient to display age on the plot, and various shapes for time and sex to differentiate between variables.

Lastly, I plotted the relative abundance of both the nasal and oral samples to observe the compositional differences between these two sample sites. I chose to visualize the top 20 genera within the samples to narrow down the most common microbial organisms. From these plots, I am able to observe and assess the changes in microbial dominance within both the nasal and oral sample sites.

Sidenote: I have also tried to reproduce alpha diversity plots, to observe diversity changes within both the nasal and oral swab samples. When I run the code, there is a warning message that states there are no singletons within my data. I spoke to a lab member and she stated this was a result of the DADA2 processing. We are aiming to see how we can fix this, but if the lines are run individually or seperately from the rest of the code, the plots will still be produced.
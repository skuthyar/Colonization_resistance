# Colonization_resistance

This repository includes the files needed to analyze the effect of pathogenic and commensal invaders on different gut microbiome communities. 

The qPCR folder includes the following: a file of invader abundances in media types, the R code used to test if there were differences, a file of invader abundances in different microbial communities post invasion, and a R markdown file that was used for microbial abundance analyses. 

The 16S control folder includes the following: the metadata file, the ASV counts file, the ASV taxonomy file (both of which were generated with dada2), and the scripts for testing if initial microbiome composition and diversity predicted invader susceptibility. The 16S challenged folder uses the same metadata, ASV counts, and ASV taxonomy files as the 16S control folder but includes scripts that analyze if and how microbial communities changed after invasion. 

The machine learning folder includes the following: the phyloseq object and a R markdown file that was used for the random forest model. 

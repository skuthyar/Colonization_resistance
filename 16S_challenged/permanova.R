# E.coli
E.coli <- prune_samples(sample_data(CR_16S_clean)$Challenger == "E.coli", CR_16S_clean)
E.coli_df <- sample_data(E.coli)
to.remove <-c("AL1_repeat", "AL2", "CA3", "CA7", "CA8", "CalPoly11_rep", "CalPoly6_rep", "NC_12_rep", "NC1_10_rep", "NC1_2_rep", "NC1_6", "NC2_1", "NC2_16", "NC2_5", "NC2_6_rep", "SM3", "SM4", "SM5_rep", "SSR1_repeat", "SSR3_rep", "UCDavis_13_rep", "UCDavis11", "VT1_rep", "VT3", "VT6_rep")
E.coli_clean <- subset_samples(E.coli, !(Individual %in% to.remove))

Domcat <- sample_data(E.coli_clean)$Domcat
Treatment <- sample_data(E.coli_clean)$Treatment
Individual <- sample_data(E.coli_clean)$Individual
Individual <- as.factor(Individual)
Population <- sample_data(E.coli_clean)$Population
Challenger_abundance <- sample_data(E.coli_clean)$Challenger_abundance_cells
Challenger_abundance <- as.numeric(Challenger_abundance)

perm_Ecoli <- adonis2(distance(E.coli_clean, method="bray") ~ Treatment + Domcat + Challenger_abundance + Population, strata=Individual, na.action = na.omit)
perm_Ecoli

# S.enterica
Sal <- prune_samples(sample_data(CR_16S_clean)$Challenger == "S.enterica", CR_16S_clean)
Domcat <- sample_data(Sal)$Domcat
Treatment <- sample_data(Sal)$Treatment
Individual <- sample_data(Sal)$Individual
Individual <- as.character(Individual)
Population <- sample_data(Sal)$Population
Challenger_abundance <- sample_data(Sal)$Challenger_abundance_cells
Challenger_abundance <- as.numeric(Challenger_abundance)

perm_Sal_int <- adonis2(distance(Sal, method="bray") ~ Treatment*Domcat + Challenger_abundance + Population, strata=Individual)
perm_Sal_int

# L.reuteri
Lacto <- prune_samples(sample_data(CR_16S_clean)$Challenger == "L.reuteri", CR_16S_clean)
Domcat <- sample_data(Lacto)$Domcat
Treatment <- sample_data(Lacto)$Treatment
Individual <- sample_data(Lacto)$Individual
Individual <- as.character(Individual)
Population <- sample_data(Lacto)$Population
Challenger_abundance <- sample_data(Lacto)$Challenger_abundance_cells
Challenger_abundance <- as.numeric(Challenger_abundance)
perm_Lacto_int <- adonis2(distance(Lacto, method="bray") ~ Treatment*Domcat + Challenger_abundance + Population, strata=Individual)
perm_Lacto_int

# remove ASVs belonging to invader
# E.coli: ASV9
# S.enterica: ASV8 
# L.reuteri: ASV2

#Change in alpha diversity after invasion
diversity_postinvasion <- read.csv("Shannon_diversity_postinvasion.csv")
diversity_postinvasion$Change_Shannon <- as.numeric(diversity_postinvasion$Change_Shannon)
diversity_postinvasion_clean <- subset(diversity_postinvasion, diversity_postinvasion$Change_Shannon != "N/A")
diversity_postinvasion_clean$Challenger <- factor(diversity_postinvasion_clean$Challenger, levels=c("E.coli", "S.enterica","L.reuteri"))

Shannon_postinvasion <- ggplot(diversity_postinvasion_clean, aes(Challenger, Change_Shannon, fill=Domcat)) + geom_boxplot(outlier.shape=NA) + 
  theme(axis.title.x=element_blank()) + ylab("Change in Shannon diversity post invasion") + 
  theme_bw() + theme(axis.title.x = element_blank()) + geom_point(aes(color=Domcat), position=pd, alpha=0.4, size=3) + geom_point(shape = 1,position=pd, size=3,colour = "black") + theme(text = element_text(size = 25)) + theme(legend.title=element_blank()) + geom_hline(yintercept=0, linetype="dashed")

# E.coli control vs. challenge without E.coli ASV
Ecoli_div <- subset(diversity_postinvasion, diversity_postinvasion$Challenger == "E.coli")
# remove treatment = NA
Ecoli_div_clean <- subset(Ecoli_div, Ecoli_div$Treatment != "NA")
wilcox.test(Shannon ~ Treatment, Ecoli_div_clean, paired=TRUE) 

# S.enterica control vs. challenged without S.enterica ASV
Sal_div <- subset(diversity_postinvasion, diversity_postinvasion$Challenger == "S.enterica")
# remove frw_CA_12_patho_S_S83
Sal_div_clean <- subset(Sal_div, Sal_div$X != "frw_CA_12_control_S_S127")
wilcox.test(Shannon ~ Treatment, Sal_div_clean, paired=TRUE) 

# L.reuteri control vs. challenged without L.reuteri ASV
Lacto_div <- subset(diversity_postinvasion, diversity_postinvasion$Challenger == "L.reuteri")
# remove res_UCDavis_11_control_L_S30
Lacto_div_clean <- subset(Lacto_div, Lacto_div$X != "res_UCDavis_11_control_L_S30")
wilcox.test(Shannon ~ Treatment, Lacto_div_clean, paired=TRUE)

# Cloud plots for alternative visualization
df1_clean$Challenger.Treatment = paste(df1_clean$Challenger,"-",df1_clean$Treatment)
multi.group <- 
  df1_clean %>%
  dabest(Challenger.Treatment, Shannon, 
         idx = list(c("E.coli - control", "E.coli - challenged", "E.coli - abx_challenged"), 
                    c("S.enterica - control", "S.enterica - challenged", "S.enterica - abx_challenged"),
                    c("L.reuteri - control", "L.reuteri - challenged", "L.reuteri - abx_challenged")),
         paired = FALSE
  )

multi.group.mean_diff <- multi.group %>% mean_diff() 

plot(multi.group.mean_diff, color.column = Challenger.Treatment, axes.title.fontsize = 10)

df2_clean$Challenger.Treatment = paste(df2_clean$Challenger,"-",df2_clean$Treatment)
multi.group <- 
  df2_clean %>%
  dabest(Challenger.Treatment, Richness, 
         idx = list(c("E.coli - control", "E.coli - challenged", "E.coli - abx_challenged"), 
                    c("S.enterica - control", "S.enterica - challenged", "S.enterica - abx_challenged")),
         paired = FALSE
  )

multi.group.mean_diff <- multi.group %>% mean_diff() 

plot(multi.group.mean_diff, color.column = Challenger.Treatment, axes.title.fontsize = 10)

# domestication groups within each challenger
# S.enterica
df1_Sal <- subset(df1_clean, df1_clean$Challenger == "S.enterica")
df1_Sal$Domcat.Treatment = paste(df1_Sal$Domcat,"-",df1_Sal$Treatment)
multi.group <- 
  df1_Sal %>%
  dabest(Domcat.Treatment, Shannon, 
         idx = list(c("free-ranging domestic - control", "free-ranging domestic - challenged", "free-ranging domestic - abx_challenged"), 
                    c("wild - control", "wild - challenged", "wild - abx_challenged"),
                    c("research - control","research - challenged","research - abx_challenged"),
                    c("industrial - control","industrial - challenged", "industrial - abx_challenged")),
         paired = FALSE
  )

multi.group.mean_diff <- multi.group %>% mean_diff() 

plot(multi.group.mean_diff, color.column = Domcat.Treatment, axes.title.fontsize = 10)

# L.reuteri
df1_Lacto <- subset(df1_clean, df1_clean$Challenger == "L.reuteri")
df1_Lacto$Domcat.Treatment = paste(df1_Lacto$Domcat,"-",df1_Lacto$Treatment)
multi.group <- 
  df1_Lacto %>%
  dabest(Domcat.Treatment, Shannon, 
         idx = list(c("free-ranging domestic - control", "free-ranging domestic - challenged", "free-ranging domestic - abx_challenged"), 
                    c("wild - control", "wild - challenged", "wild - abx_challenged"),
                    c("research - control","research - challenged","research - abx_challenged"),
                    c("industrial - control","industrial - challenged", "industrial - abx_challenged")),
         paired = FALSE
  )

multi.group.mean_diff <- multi.group %>% mean_diff() 

plot(multi.group.mean_diff, color.column = Domcat.Treatment, axes.title.fontsize = 10)

# Estimates of alpha diversity
CR_16S_nofecal <- prune_samples(sample_data(CR_16S_clean)$Treatment != "fecal", CR_16S_clean)
CR_16S_nofecal_noabx <- prune_samples(sample_data(CR_16S_nofecal)$Treatment != "abx_challenged", CR_16S_nofecal)
richness_CR <- CR_16S_nofecal_noabx %>% breakaway
plot(richness_CR, physeq=CR_16S_nofecal_noabx, color="Treatment", shape = "Domcat")
richness_CR_df <- summary(richness_CR) %>% as_tibble

meta <- CR_16S_nofecal_noabx %>%
  sample_data %>%
  as_tibble %>%
  mutate("sample_names" = CR_16S_nofecal_noabx %>% sample_names )

combined_richness <- meta %>%
  left_join(summary(richness_CR),
            by = "sample_names")
#divnet
fam <- CR_16S_nofecal_noabx %>%
  tax_glom(taxrank="Family")
dv <- DivNet::divnet(fam, X=NULL) 

combined_shannon <- meta %>%
  dplyr::left_join(dv$shannon %>% summary,
                   by = "sample_names")
combined_shannon

bt_day_fixed_id_random <- betta_random(formula = estimate ~ Treatment + Domcat + Challenger_abundance_cells | Individual, 
                                       ses = error,  data = combined_shannon)
bt_day_fixed_id_random$table

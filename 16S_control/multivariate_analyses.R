control <- prune_samples(sample_data(CR_16S_clean)$Treatment == "control", CR_16S_clean)
control_fam <- tax_glom(control, taxrank="Family")
sample_data(control_fam)$Shannon <- estimate_richness(control_fam, measure=c("Shannon"))$Shannon
sample_data(control_fam)$Richness <- estimate_richness(control_fam, measure=c("Observed"))$Observed
control_rel  = transform_sample_counts(control_fam, function(x) x / sum(x) )

control_rel_sub <- subset_taxa(control_rel, Family=="Enterobacteriaceae" | Family=="Lactobacillaceae")
control_rel_sub_df <- psmelt(control_rel_sub)

# E.coli
control_Ecoli_df <- subset(control_rel_sub_df, control_rel_sub_df$Challenger == "E.coli")
control_Ecoli_df <- subset(control_Ecoli_df, control_Ecoli_df$Family != "Lactobacillaceae")
control_Ecoli_df$Challenged_ResAbund <- as.numeric(control_Ecoli_df$Challenged_ResAbund)
control_Ecoli_df_clean <- subset(control_Ecoli_df, control_Ecoli_df$Challenged_ResAbund != "N/A")
mean(control_Ecoli_df_clean$Challenged_ResAbund)
var(control_Ecoli_df_clean$Challenged_ResAbund)
control_Ecoli_df_clean$Challenged_ResAbund <- round(control_Ecoli_df_clean$Challenged_ResAbund)

# Fit a Negative Binomial mixed effects model
nb_model <- glmmTMB(Challenged_ResAbund ~ Abundance + Shannon + 
                    (1|Population), 
                    data = control_Ecoli_df_clean, 
                    family = nbinom2)  # 'nbinom2' is for Negative Binomial distribution

summary(nb_model)
confint(nb_model)

# S.enterica
control_Sal_df <- subset(control_rel_sub_df, control_rel_sub_df$Challenger == "S.enterica")
control_Sal_df <- subset(control_Sal_df, control_Sal_df$Family != "Lactobacillaceae")
control_Sal_df <- subset(control_Sal_df, control_Sal_df$Sample != "frd_VT_7_control_S_S107")
control_Sal_df$Challenged_ResAbund <- as.numeric(control_Sal_df$Challenged_ResAbund)
control_Sal_df_clean <- subset(control_Sal_df, control_Sal_df$Challenged_ResAbund != "N/A")
mean(control_Sal_df_clean$Challenged_ResAbund)
var(control_Sal_df_clean$Challenged_ResAbund)
control_Sal_df_clean$Challenged_ResAbund <- round(control_Sal_df_clean$Challenged_ResAbund)

# check for multicollinearity: 
vif(lm(Challenged_ResAbund ~ Abundance + Domcat + Shannon, data = control_Sal_df_clean))

# fit a negative binomial mixed effects model
nb_model2 <- glmmTMB(Challenged_ResAbund ~ Abundance + Domcat + Shannon + 
                     (1|Population), 
                     data = control_Sal_df_clean, 
                     family = nbinom2)

summary(nb_model2)
confint(nb_model2)

# no wild
control_Sal_df_nowild <- subset(control_Sal_df, control_Sal_df$Domcat != "wild")
nb_model2_nowild <- glmmTMB(Challenged_ResAbund ~ Abundance + Domcat + Shannon + 
                     (1|Population), 
                     data = control_Sal_df_nowild, 
                     family = nbinom2)

summary(nb_model2_nowild)
confint(nb_model2_nowild)

# L.reuteri 
control_Lacto_df <- subset(control_rel_sub_df, control_rel_sub_df$Challenger == "L.reuteri")
control_Lacto_df <- subset(control_Lacto_df, control_Lacto_df$Family == "Lactobacillaceae")
control_Lacto_df$Challenged_ResAbund <- as.numeric(control_Lacto_df$Challenged_ResAbund)
control_Lacto_df_clean <- subset(control_Lacto_df, control_Lacto_df$Challenged_ResAbund != "N/A")
mean(control_Lacto_df_clean$Challenged_ResAbund)
var(control_Lacto_df_clean$Challenged_ResAbund)
control_Lacto_df_clean$Challenged_ResAbund <- round(control_Lacto_df_clean$Challenged_ResAbund)

# check for multicollinearity: 
vif(lm(Challenged_ResAbund ~ Abundance + Domcat + Shannon, data = control_Lacto_df_clean))

# fit a negative binomial mixed effects model
nb_model3 <- glmmTMB(Challenged_ResAbund ~ Abundance + Domcat + Shannon + 
                     (1|Population), 
                     data = control_Lacto_df_clean, 
                     family = nbinom2)

summary(nb_model3)
confint(nb_model3)

# no wild
control_Lacto_df_nowild <- subset(control_Lacto_df, control_Lacto_df$Domcat != "wild")
nb_model3_nowild <- glmmTMB(Challenged_ResAbund ~ Abundance + Shannon + 
                     (1|Population), 
                     data = control_Lacto_df_nowild, 
                     family = nbinom2)

summary(nb_model3_nowild)
confint(nb_model3_nowild)

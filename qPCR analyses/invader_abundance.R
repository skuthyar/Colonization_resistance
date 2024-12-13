full_qpcr <- read_excel("all_metadata_samesample.xlsx")
full_qpcr_clean <- subset(full_qpcr, full_qpcr$Domcat != "N/A")
full_qpcr_clean <- subset(full_qpcr_clean, full_qpcr_clean$Treatment != "fecal")
full_qpcr_clean <- subset(full_qpcr_clean, full_qpcr_clean$Treatment != "ABX")
full_qpcr_clean$Treatment <- factor(full_qpcr_clean$Treatment, levels=c("control","challenged","abx_challenged"))
full_qpcr_clean$Domcat <- factor(full_qpcr_clean$Domcat, levels=c("wild","free-ranging domestic","research","industrial"))
full_qpcr_clean$Challenger <- factor(full_qpcr_clean$Challenger, levels=c("E.coli","S.enterica","L.reuteri"))
full_qpcr_clean$Challenger_abundance <- as.numeric(full_qpcr_clean$Challenger_abundance)
full_qpcr_clean$Challenger_abundance_cells <- as.numeric(full_qpcr_clean$Challenger_abundance_cells)

# Create Figure 2A
all <- ggplot(full_qpcr_clean, aes(Treatment, Challenger_abundance_cells, fill=Treatment)) + geom_boxplot(outlier.shape=NA) + geom_point(aes(color=Treatment), position=pd, alpha=0.4, size=3) + geom_point(shape = 1,position=pd, size=3,colour = "black") + facet_wrap(~Challenger) + ylab("Abundance of invader (cells)") + theme_bw() + 
  scale_y_log10() + geom_hline(yintercept=1000000, linetype="dashed") + theme(legend.title = element_blank()) + theme(axis.title.x = element_blank()) + theme(legend.position="none") + 
  theme(text = element_text(size = 25)) + scale_x_discrete(labels=c("Control","Challenged","ABX+Challenged")) + scale_fill_brewer(palette="Paired") + scale_color_brewer(palette="Paired")
all

# Statistics for control vs. challenged 
full_qpcr_noabx <- subset(full_qpcr_clean, full_qpcr_clean$Treatment != "abx_challenged")
# E.coli
Ecoli_noabx <- subset(full_qpcr_noabx, full_qpcr_noabx$Challenger == "E.coli")
to.remove <- c("frw_AL_1_patho_repeat_S194", "frw_CA_3_control", "frw_CA_8_control", "res_CalPoly_6_rep_control", "ind_NC1_10_rep_control", "ind_NC2_6_rep_control", "frd_SSR_1_control_repeat", "frd_SSR_3_rep_control", "frd_VT_3_control", "frd_VT_6_rep_control")
Ecoli_clean <- Ecoli_noabx[!(Ecoli_noabx$Sampleid %in% to.remove), ]
wilcox.test(Challenger_abundance_cells ~ Treatment, Ecoli_clean, paired=TRUE, p.adjust.method="bonferroni", alternative = "two.sided") 

# S.enterica
Sal_noabx <- subset(full_qpcr_noabx, full_qpcr_noabx$Challenger == "S.enterica")
# remove frd_VT_7_control_S
Sal_noabx_clean <- subset(Sal_noabx, Sal_noabx$Sampleid != "frd_VT_7_control_S")
wilcox.test(Challenger_abundance_cells ~ Treatment, Sal_noabx_clean, paired=TRUE, p.adjust.method="bonferroni", alternative = "two.sided")

# L.reuteri
Lacto_noabx <- subset(full_qpcr_noabx, full_qpcr_noabx$Challenger == "L.reuteri")
# remove res_UCDavis_11_control_L
Lacto_noabx_clean <- subset(Lacto_noabx, Lacto_noabx$Sampleid != "res_UCDavis_11_control_L")
wilcox.test(Challenger_abundance_cells ~ Treatment, Lacto_noabx_clean, paired=TRUE, p.adjust.method="bonferroni", alternative = "two.sided") 

# Statistics for challenged vs. ABX+challenged
full_qpcr_nocontrol <- subset(full_qpcr_clean, full_qpcr_clean$Treatment != "control")

# E.coli
Ecoli_nocontrol <- subset(full_qpcr_nocontrol, full_qpcr_nocontrol$Challenger == "E.coli")
to.remove <- c("frw_CA_3_ABX_patho", "frw_CA_8_ABX_patho", "res_CalPoly_6_rep_ABX_patho", "ind_NC1_10_rep_ABX_patho", "ind_NC2_6_rep_ABX_patho", "frd_SSR_3_rep_ABX_patho", "frd_VT_3_ABX_patho", "frd_VT_6_rep_ABX_patho")
Ecoli_clean <- Ecoli_nocontrol[!(Ecoli_nocontrol$Sampleid %in% to.remove), ]
wilcox.test(Challenger_abundance_cells ~ Treatment, Ecoli_clean, paired=TRUE, p.adjust.method="bonferroni") 
# S.enterica
Sal_nocontrol <- subset(full_qpcr_nocontrol, full_qpcr_nocontrol$Challenger == "S.enterica")
# remove frd_VT_7_ABX_patho_S
Sal_clean <- subset(Sal_nocontrol, Sal_nocontrol$Sampleid != "frd_VT_7_ABX_patho_S")
wilcox.test(Challenger_abundance_cells ~ Treatment, Sal_clean, paired=TRUE, p.adjust.method="bonferroni", alternative = "two.sided")

# L.reuteri
Lacto_nocontrol <- subset(full_qpcr_nocontrol, full_qpcr_nocontrol$Challenger == "L.reuteri")
# remove res_UCDavis_11_ABX_patho_L
Lacto_clean <- subset(Lacto_nocontrol, Lacto_nocontrol$Sampleid != "res_UCDavis_11_ABX_patho_L")
wilcox.test(Challenger_abundance_cells ~ Treatment, Lacto_clean, paired=TRUE, p.adjust.method="bonferroni", alternative = "two.sided") 

full_qpcr_clean_noecoli <- subset(full_qpcr_clean, full_qpcr_clean$Challenger != "E.coli")
# Create Figure 2B
challenged_genetics <- ggplot(full_qpcr_clean_noecoli, aes(Treatment, Challenger_abundance_cells, fill=Genetics)) + geom_boxplot(outlier.shape=NA) + geom_point(aes(color=Genetics), position=pd, alpha=0.4, size=3) + geom_point(shape = 1,position=pd, size=3,colour = "black")  + theme_bw() + geom_hline(yintercept=1000000, linetype="dashed") +
  theme(legend.title = element_blank()) + theme(axis.title.x = element_blank()) + ylab("Abundance of invader (cells)") + theme(legend.position="bottom") + theme(text = element_text(size = 25)) + scale_fill_manual(values=c("wild"="#36651e","domestic"="#CD7F32")) + scale_color_manual(values=c("wild"="#36651e","domestic"="#CD7F32")) + facet_wrap(~Challenger) + scale_y_log10() + scale_x_discrete(labels=c("Control","Challenged","ABX+Challenged"))
challenged_genetics

# S.enterica - control domestic vs. challenged domestic
full_qpcr_noabx_noecoli <- subset(full_qpcr_clean_noecoli, full_qpcr_clean_noecoli$Treatment != "abx_challenged")
Sal <- subset(full_qpcr_noabx_noecoli, full_qpcr_noabx_noecoli$Challenger == "S.enterica")
Sal_dom <- subset(Sal, Sal$Genetics == "domestic")
# remove frd_VT_7_control_S
Sal_dom_clean <- subset(Sal_dom, Sal_dom$Sampleid != "frd_VT_7_control_S")
wilcox.test(Challenger_abundance_cells ~ Treatment, Sal_dom_clean, paired=TRUE, p.adjust.method="bonferroni") 

# S.enterica - control wild vs. challenged wild
Sal_wild <- subset(Sal, Sal$Genetics == "wild")
wilcox.test(Challenger_abundance_cells ~ Treatment, Sal_wild, paired=TRUE, p.adjust.method="bonferroni") 

# S.enterica - challenged, wild vs. domestic
Sal_challenged <- subset(Sal, Sal$Treatment == "challenged")
wilcox.test(Challenger_abundance_cells ~ Genetics, Sal_challenged, p.adjust.method="bonferroni")

# L.reuteri - control domestic vs. challenged domestic
Lacto <- subset(full_qpcr_noabx_noecoli, full_qpcr_noabx_noecoli$Challenger == "L.reuteri")
Lacto_dom <- subset(Lacto, Lacto$Genetics == "domestic")
# remove res_UCDavis_11_control_L
Lacto_dom_clean <- subset(Lacto_dom, Lacto_dom$Sampleid != "res_UCDavis_11_control_L")
wilcox.test(Challenger_abundance_cells ~ Treatment, Lacto_dom_clean, paired=TRUE, p.adjust.method="bonferroni")

# L.reuteri - control wild vs. challenged wild
Lacto_wild <- subset(Lacto, Lacto$Genetics == "wild")
wilcox.test(Challenger_abundance_cells ~ Treatment, Lacto_wild, paired=TRUE, p.adjust.method="bonferroni")

# L.reuteri - control domestic vs. wild
Lacto_control <- subset(Lacto, Lacto$Treatment == "control")
wilcox.test(Challenger_abundance_cells ~ Genetics, Lacto_control) 

# L.reuteri - challenged wild vs. domestic
Lacto_challenged <- subset(Lacto, Lacto$Treatment == "challenged")
wilcox.test(Challenger_abundance_cells ~ Genetics, Lacto_challenged, p.adjust.method="bonferroni")

full_qpcr_clean_noecoli_nowild <- subset(full_qpcr_clean_noecoli, full_qpcr_clean_noecoli$Domcat != "wild")
# Create Figure 2C
treatment_domcat <- ggplot(full_qpcr_clean_noecoli_nowild, aes(Treatment, Challenger_abundance_cells, fill=Domcat)) + geom_boxplot(outlier.shape=NA) + geom_point(aes(color=Domcat), position=pd, alpha=0.4, size=3) + geom_point(shape = 1,position=pd, size=3,colour = "black") + facet_wrap(~Challenger)
treatment_domcat <- treatment_domcat + theme_bw() + geom_hline(yintercept=1000000, linetype="dashed")
treatment_domcat <- treatment_domcat + theme(legend.title = element_blank()) + theme(axis.title.x = element_blank()) + theme(axis.title.y = element_blank()) + theme(legend.position="bottom") + scale_y_log10()
treatment_domcat <- treatment_domcat + theme(text = element_text(size = 25)) + scale_x_discrete(labels=c("Control","Challenged", "ABX+Challenged")) + scale_fill_manual(values = c("industrial"="#963a3a","research"="#e8aa23","free-ranging domestic"="#eb9696")) + scale_color_manual(values = c("industrial"="#963a3a","research"="#e8aa23","free-ranging domestic"="#eb9696"))
treatment_domcat 

# Statistical analyses for control vs. challenged
full_qpcr_clean_noecoli_nowild_noabx <- subset(full_qpcr_clean_noecoli_nowild, full_qpcr_clean_noecoli_nowild$Treatment != "abx_challenged")
# S.enterica - frd
Sal <- subset(full_qpcr_clean_noecoli_nowild_noabx, full_qpcr_clean_noecoli_nowild_noabx$Challenger == "S.enterica")
Sal_frd <- subset(Sal, Sal$Domcat == "free-ranging domestic")
Sal_frd_clean <- subset(Sal_frd, Sal_frd$Sampleid != "frd_VT_7_control_S")
wilcox.test(Challenger_abundance_cells ~ Treatment, Sal_frd_clean, paired=TRUE, p.adjust.method="bonferroni")

# S.enterica - research
Sal_res <- subset(Sal, Sal$Domcat == "research")
wilcox.test(Challenger_abundance_cells~ Treatment, Sal_res, paired=TRUE, p.adjust.method="bonferroni") 

# S.enterica - industrial
Sal_ind <- subset(Sal, Sal$Domcat == "industrial")
wilcox.test(Challenger_abundance_cells ~ Treatment, Sal_ind, paired=TRUE, p.adjust.method="bonferroni") 

# L.reuteri - frd
Lacto <- subset(full_qpcr_clean_noecoli_nowild_noabx, full_qpcr_clean_noecoli_nowild_noabx$Challenger == "L.reuteri")
Lacto_frd <- subset(Lacto, Lacto$Domcat == "free-ranging domestic")
wilcox.test(Challenger_abundance_cells ~ Treatment, Lacto_frd, paired=TRUE, p.adjust.method="bonferroni") 

# L.reuteri - research
Lacto_res <- subset(Lacto, Lacto$Domcat == "research")
# remove res_UCDavis_11_ABX_patho_L
Lacto_res_clean <- subset(Lacto_res, Lacto_res$Sampleid != "res_UCDavis_11_control_L")
wilcox.test(Challenger_abundance_cells ~ Treatment, Lacto_res_clean, paired=TRUE, p.adjust.method="bonferroni") 

# L.reuteri - industrial
Lacto_ind <- subset(Lacto, Lacto$Domcat == "industrial")
wilcox.test(Challenger_abundance_cells ~ Treatment, Lacto_ind, paired=TRUE, p.adjust.method="bonferroni") 


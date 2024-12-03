# Shannon diversity 
sample_data(CR_16S_clean)$Shannon <- estimate_richness(CR_16S_clean, measure=c("Shannon"))$Shannon
df1 <- data.frame(sample_data(CR_16S_clean))
df1_clean <- subset(df1, df1$Treatment != "fecal")
df1_clean <- subset(df1_clean, df1_clean$Challenger !="N/A")

df1_control <- subset(df1_clean, df1_clean$Treatment == "control")
df1_control <- subset(df1_control, df1_control$Challenged_ResAbund != "N/A")
df1_control$Challenged_ResAbund <- as.numeric(df1_control$Challenged_ResAbund)
df1_control$Domcat <- factor(df1_control$Domcat, levels=c("wild","free-ranging domestic", "research","industrial"))

# E.coli
df1_control_Ecoli <- subset(df1_control, df1_control$Challenger == "E.coli")
df1_control_Ecoli$Challenged_ResAbund <- as.numeric(df1_control_Ecoli$Challenged_ResAbund)
df1_control_Ecoli_clean <- subset(df1_control_Ecoli, df1_control_Ecoli$Challenged_ResAbund != 0)
Shannon_inv_abundance_Ecoli <- ggplot(df1_control_Ecoli, aes(Shannon, log(Challenged_ResAbund))) + geom_point(size=3) + theme_bw() + stat_poly_line() + theme(axis.title.x = element_blank()) + 
  ylab("E.coli abundance (challenged communities (log))") + theme(text = element_text(size = 20)) 
Shannon_inv_abundance_Ecoli 
cor.test(df1_control_Ecoli$Challenged_ResAbund, df1_control_Ecoli$Shannon,
                                           method = "spearman", exact=FALSE, na.action=na.exclude) 

# S.enterica
df1_control_Sal <- subset(df1_control, df1_control$Challenger == "S.enterica")
df1_control_Sal$Challenged_ResAbund <- as.numeric(df1_control_Sal$Challenged_ResAbund)
Shannon_inv_abundance_Sal <- ggplot(df1_control_Sal, aes(Shannon, Challenged_ResAbund)) + geom_point(size=4) + theme_bw() + scale_y_log10() + theme(axis.title.x = element_blank()) + 
  ylab("S.enterica abundance (challenged)") + theme(text = element_text(size = 20))
Shannon_inv_abundance_Sal 
cor.test(df1_control_Sal$Challenged_ResAbund, df1_control_Sal$Shannon,
                                           method = "spearman", exact=FALSE, na.action=na.exclude) 

# by domestication context 
Shannon_inv_abundance_Sal_domcat <- ggplot(df1_control_Sal, aes(Shannon, Challenged_ResAbund, color=Domcat)) + geom_point(size=4) + scale_y_log10() + theme_bw() + stat_poly_line(data=subset(df1_control_Sal, Domcat=="industrial"), aes(Shannon,Challenged_ResAbund,color=factor(Domcat)),se=FALSE) + 
 theme(text = element_text(size = 20)) + ylab("S.enterica abundance (challenged)") + theme(axis.title.x=element_blank()) + theme(legend.title=element_blank()) + scale_color_manual(values=c("wild"="#36651e", "industrial"="#963a3a","research"="#e8aa23","free-ranging domestic"="#eb9696"))
Shannon_inv_abundance_Sal_domcat 
# Spearman for frd
frd <- subset(df1_control_Sal, df1_control_Sal$Domcat == "free-ranging domestic")
cor.test(frd$Challenged_ResAbund, frd$Shannon,
                                         method = "spearman", p.adjust.method="BH") 
# Spearman for ind
ind <- subset(df1_control_Sal, df1_control_Sal$Domcat == "industrial")
cor.test(ind$Challenged_ResAbund, ind$Shannon,
                                 method = "spearman", p.adjust.method="BH") 
# Spearman for res
res <- subset(df1_control_Sal, df1_control_Sal$Domcat == "research")
cor.test(res$Challenged_ResAbund, res$Shannon,
                                 method = "spearman", p.adjust.method="BH") 
# Spearman for wild
wild <- subset(df1_control_Sal, df1_control_Sal$Domcat == "wild")
cor.test(wild$Challenged_ResAbund, wild$Shannon,
                                 method = "spearman", p.adjust.method="BH") 

# L.reuteri
df1_control_Lacto <- subset(df1_control, df1_control$Challenger == "L.reuteri")
df1_control_Lacto$Challenged_ResAbund <- as.numeric(df1_control_Lacto$Challenged_ResAbund)
Shannon_inv_abundance_Lacto <- ggplot(df1_control_Lacto, aes(Shannon, Challenged_ResAbund)) + geom_point(size=4) + theme_bw() + scale_y_log10() +
  ylab("L.reuteri abundance (challenged)") + theme(axis.title.x=element_blank()) + theme(text = element_text(size = 20)) 
Shannon_inv_abundance_Lacto 
cor.test(df1_control_Lacto$Challenged_ResAbund, df1_control_Lacto$Shannon,
                                         method = "spearman",  p.adjust.method="BH", exact=FALSE) 

# by domestication context
Shannon_inv_abundance_Lacto_domcat <- ggplot(df1_control_Lacto, aes(Shannon, Challenged_ResAbund, color=Domcat)) + geom_point(size=4) + scale_y_log10() + theme_bw() + stat_poly_line(data=subset(df1_control_Lacto, Domcat=="free-ranging domestic" | Domcat=="industrial"), aes(Shannon,Challenged_ResAbund,color=factor(Domcat)),se=FALSE) + 
  theme(text = element_text(size = 20)) + ylab("L.reuteri abundance (challenged)") + theme(axis.title.x=element_blank()) + theme(legend.title=element_blank()) + scale_color_manual(values=c("wild"="#36651e", "industrial"="#963a3a","research"="#e8aa23","free-ranging domestic"="#eb9696"))
Shannon_inv_abundance_Lacto_domcat 

# Spearman for frd
frd <- subset(df1_control_Lacto, df1_control_Lacto$Domcat == "free-ranging domestic")
cor.test(frd$Challenged_ResAbund, frd$Shannon,
         method = "spearman", exact=FALSE, p.adjust.method="BH") 

# Spearman for ind
ind <- subset(df1_control_Lacto, df1_control_Lacto$Domcat == "industrial")
cor.test(ind$Challenged_ResAbund, ind$Shannon,
         method = "spearman", exact=FALSE, p.adjust.method="BH") 

# Spearman for res
res <- subset(df1_control_Lacto, df1_control_Lacto$Domcat == "research")
cor.test(res$Challenged_ResAbund, res$Shannon,
         method = "spearman", exact=FALSE, p.adjust.method="BH") 

# Spearman for wild
wild <- subset(df1_control_Lacto, df1_control_Lacto$Domcat == "wild")
cor.test(wild$Challenged_ResAbund, wild$Shannon,
         method = "spearman", exact=FALSE, na.action=na.exclude, p.adjust.method="BH") 

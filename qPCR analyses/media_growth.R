# Read in file
media_growth <- read_excel("media_microbe_qpcr.xlsx")

# Make variables ready for plotting
media_growth$Bacterial_cells_adj <- as.numeric(media_growth$Bacterial_cells_adj)
media_growth$Media <- factor(media_growth$Media, levels=c("PBS","regular_media","abx_media"))
media_growth$Microbe <- factor(media_growth$Microbe, levels= c("E.coli","S.enterica","L.reuteri"))
media_growth_noPBS <- subset(media_growth, media_growth$Media != "PBS")

# Create Figure 1B
my_comp_media <- list(c("regular_media","abx_media"))
media_growth_noPBS_plot <- ggplot(media_growth_noPBS, aes(Media, Bacterial_cells_adj, fill=Media)) + geom_boxplot(outlier.shape=NA) + geom_point(position = position_jitterdodge(jitter.width = 0.1), size=4, alpha=0.6) + theme_bw() + scale_y_log10() + ylab("Abundance of invader (cells)") + theme(axis.title.x=element_blank()) + theme(text = element_text(size = 20)) 
media_growth_noPBS_plot <- media_growth_noPBS_plot + scale_x_discrete(labels=c("Regular media", "Antibiotic media")) + theme(legend.position="none") 
media_growth_noPBS_plot <- media_growth_noPBS_plot  + scale_fill_manual(values=c("#dfc27d","#a6611a")) + stat_compare_means(comp=my_comp_media) + facet_wrap(~Microbe)
media_growth_noPBS_plot

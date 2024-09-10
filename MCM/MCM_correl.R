
pattern <- 'MCM_D'
res.dirs <- grep(pattern, list.dirs('.', T, F), value = T)

res.dir      <- 'MCM_D6331_Emu_EMUdb'

plot_all <- data.table(NULL) 

for (res.dir in res.dirs[]) {

  file     <- paste0(res.dir, '/Observed-Theoretical Correlations.tsv')
  message('Importing: ', file)
  plot_sum <- fread(file)
  plot_sum[,source := file]
  #plot_sum   <- rbind(correl.exclunclass, correl.keepunclass)
  print(plot_sum[,.N,sample])
  
  plot_all <- rbind(plot_sum, plot_all, fill=T)
}

plot_all[,rank := factor(rank, levels = ranks)]


plot_sum   <- plot_all
sampl_freq <- plot_sum[,.N,sample]
sampl_freq <- plot_sum[,.N,.(GS_version, platform, db, workflow, Vregion, rank, sample_name, R2)]

plot_means <- unique(plot_sum[,.(GS_version, DNA_isolation_method, platform, db, workflow, Vregion, rank, sample_name, R2, source)])
plot_means[,rank := factor(rank, levels = ranks)]

## plotting

ggr2_B <- ggplot(plot_means[platform == 'Illumina' & rank != 'phylum'], 
               aes(x = sample_name, y = R2, fill = DNA_isolation_method)) +
  geom_col(position = position_dodge2()) +
  guides(fill = guide_legend(title = "DNA extraction method", direction = "horizontal", nrow = 1)) +
  scale_fill_manual(values=as.character(pal.man[])) +  # c(10:16)
  #facet_wrap(~ sample, scales = "free") +
  coord_cartesian(ylim = c(0,1)) +
  theme_ipsum() + 
  theme(legend.position="bottom",
        axis.text.x = element_text(angle=-60),
        panel.spacing.x = unit(0.1, 'mm')) +
  facet_nested(cols=vars(GS_version, platform, db, workflow),
               rows=vars(rank),
               # ncol=5, 
               scales = 'free') +
  ggtitle('b')


#ggsave('FIGURES_REVISION/Fig_6A_REVISION.jpg', ggr2_A, height = 12, width = 24)



ggr2_A <- ggplot(plot_means[platform == 'ONT' & rank != 'phylum'], 
                 aes(x = sample_name, y = R2, fill = DNA_isolation_method)) +
  geom_col(position = position_dodge2()) +
  guides(fill = guide_legend(title = "DNA extraction method", direction = "horizontal", nrow = 1)) +
  scale_fill_manual(values=as.character(pal.man[])) +  # c(10:16)
  #facet_wrap(~ sample, scales = "free") +
  coord_cartesian(ylim = c(0,1)) +
  theme_ipsum() + 
  theme(legend.position="bottom",
        axis.text.x = element_text(angle=-60),
        panel.spacing.x = unit(0.1, 'mm')) +
  facet_nested(cols=vars(GS_version, platform, db, workflow),
               rows=vars(rank),
               # ncol=5, 
               scales = 'free') +
  ggtitle('a')

#ggsave('FIGURES_REVISION/Fig_6B_REVISION.jpg', ggr2_B, height = 12, width = 24)



ggr2_AB <- cowplot::plot_grid(ggr2_A, ggr2_B, ncol = 1)

#ggsave('FIGURES_REVISION/Fig_6A-B_REVISION.jpg', ggr2_AB, height = 20, width = 24)



ggr2_C  <- readRDS('Portik_etal_2022/Fig_6C_REVISION.rds') + ggtitle('c') + theme(legend.position="bottom")

ggr2_ABC <- cowplot::plot_grid(ggr2_AB, cowplot::plot_grid(ggr2_C, NULL, ncol = 1), ncol = 2, rel_widths = c(3,1))

#ggsave('FIGURES_REVISION/Figure 6A-B-C REVISION.jpg', ggr2_ABC, height = 24, width = 24)



ggr2_D  <- readRDS('CAMISIM/Fig6D.ggr2.rds')


ggr2_ABCD <- cowplot::plot_grid(ggr2_AB, cowplot::plot_grid(ggr2_C, ggr2_D, ncol = 1), ncol = 2, rel_widths = c(3,1))

ggsave('FIGURES_REVISION/Figure 6A-B-C-D REVISION.jpg', ggr2_ABCD, height = 24, width = 24)


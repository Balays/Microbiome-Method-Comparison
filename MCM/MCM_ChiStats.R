
pattern <- 'MCM_D'
res.dirs <- grep(pattern, list.dirs('.', T, F), value = T)

res.dir      <- 'MCM_D6331_Emu_EMUdb'

plot_all <- data.table(NULL) 

for (res.dir in res.dirs[]) {
  
  file     <- paste0(res.dir, '/Detection statistics.tsv')
  message('Importing: ', file)
  plot_sum <- fread(file)
  plot_sum[,source := file]
  #plot_sum   <- rbind(correl.exclunclass, correl.keepunclass)
  print(plot_sum[,.N,sample])
  
  plot_all <- rbind(plot_sum, plot_all, fill=T)
}

plot_all[,level := factor(level, levels = c('Genus', 'Species'))]


## plotting

plot_chi <- data.table(pivot_wider(
  plot_all[level == 'Species' & (statistic == 'Chi_stat' | statistic == 'Chisq.test.sign'), ],
  names_from = 'statistic', values_from = 'value'
))


ggChi_A <- ggplot(plot_chi[platform == 'Illumina',] , 
                 aes(x = sample_name, y = Chi_stat, fill = DNA_isolation_method)) +
  geom_col(position = position_dodge2()) +
  geom_text(data=plot_chi[platform == 'Illumina' & Chisq.test.sign == 1,],
            aes(x = sample_name, y = Chi_stat + 10, fill = DNA_isolation_method),
            label = '*', position = position_dodge2()) +
  guides(fill = guide_legend(title = "DNA extraction method", direction = "horizontal", nrow = 1)) +
  scale_fill_manual(values=as.character(pal.man[])) +  # c(10:16)
  scale_y_continuous(name='Chi-statistic')+
  #facet_wrap(~ sample, scales = "free") +
  #coord_cartesian(ylim = c(0,1)) +
  theme_ipsum() + 
  theme(legend.position="bottom",
        axis.text.x = element_text(angle=-60),
        axis.title.x = element_blank(),
        panel.spacing.x = unit(0.1, 'mm')) +
  facet_nested(cols=vars(GS_version, platform, db, workflow),
               rows=vars(Detection_threshold_percent),
               # ncol=5, 
               scales = 'free') +
  ggtitle('a')


ggChi_B <- ggplot(plot_chi[platform == 'ONT',] , 
                  aes(x = sample_name, y = Chi_stat, fill = DNA_isolation_method)) +
  geom_col(position = position_dodge2()) +
  geom_text(data=plot_chi[platform == 'ONT' & Chisq.test.sign == 1,],
            aes(x = sample_name, y = Chi_stat + 10, fill = DNA_isolation_method),
            label = '*', position = position_dodge2()) +
  guides(fill = guide_legend(title = "DNA extraction method", direction = "horizontal", nrow = 1)) +
  scale_fill_manual(values=as.character(pal.man[])) +  # c(10:16)
  scale_y_continuous(name='Chi-statistic')+
  #facet_wrap(~ sample, scales = "free") +
  #coord_cartesian(ylim = c(0,1)) +
  theme_ipsum() + 
  theme(legend.position="bottom",
        axis.text.x = element_text(angle=-60),
        axis.title.x = element_blank(),
        panel.spacing.x = unit(0.1, 'mm')) +
  facet_nested(cols=vars(GS_version, platform, db, workflow),
               rows=vars(Detection_threshold_percent),
               # ncol=5, 
               scales = 'free') +
  ggtitle('b')


ggChi_C <- read_rds('K:/data/Portik_etal_2022/SuppFig_S10_ggChi_C.rds')


ggChi <- cowplot::plot_grid(cowplot::plot_grid(ggChi_A, ggChi_C, ncol = 2, rel_widths = c(3,1)), 
                            ggChi_B, 
                            ncol = 1)

ggsave('FIGURES_REVISION/Supp Fig S10 REVISION.jpg', ggChi, height = 24, width = 24)

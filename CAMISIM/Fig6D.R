
library(data.table)
library(ggplot2)
library(cowplot)


file     <-  '/Observed-Theoretical Correlations.tsv'
plot_sum <- fread(file)

sampl_freq <- plot_sum[,.N,sample_name]
sampl_freq <- plot_sum[,.N,.( platform, db, workflow, Vregion, rank, sample_name, R2)]

plot_means <- plot_sum[
  !sample_name %in% c('sample_10', 'sample_15', 'sample_35', 'sample20')
  ,.(mean_r2= mean(R2),
    sd_r2  = sd(R2)) ,
  by=.( platform, db, workflow, Vregion, rank)]

## plotting

ggr2_D1 <- ggplot(plot_means[rank != 'phylum'], 
               aes(x = platform, y = mean_r2, fill = platform)) +
geom_col(position = position_dodge2()) +
geom_errorbar(aes(ymin = mean_r2 - sd_r2, ymax= mean_r2 + sd_r2), width = 0.2) +
guides(fill = guide_legend(title = "Platform", direction = "horizontal", nrow = 1)) +
scale_fill_manual(values=as.character(pal.man[])) +  # c(10:16)
#facet_wrap(~ sample, scales = "free") +
coord_cartesian(ylim = c(0,1)) +
theme_ipsum() + 
theme(legend.position="bottom",
      axis.text.x = element_text(angle=-60),
      panel.spacing.x = unit(0.1, 'mm')) +
facet_nested(cols=vars(platform, db, workflow),
             rows=vars(rank),
             # ncol=5, 
             scales = 'free') +
ggtitle('d')

plot_means <- plot_sum[
  !sample_name %in% c('sample_17', 'sample_37', 'sample_57', 'sample_40', 'sample_50')
  ,.(mean_r2= mean(R2),
    sd_r2  = sd(R2)) ,
  by=.( platform, db, workflow, Vregion, rank)]

## plotting

ggr2_D2 <- ggplot(plot_means[rank != 'phylum'], 
                 aes(x = platform, y = mean_r2, fill = platform)) +
  geom_col(position = position_dodge2()) +
  geom_errorbar(aes(ymin = mean_r2 - sd_r2, ymax= mean_r2 + sd_r2), width = 0.2) +
  guides(fill = guide_legend(title = "Platform", direction = "horizontal", nrow = 1)) +
  scale_fill_manual(values=as.character(pal.man[])) +  # c(10:16)
  #facet_wrap(~ sample, scales = "free") +
  coord_cartesian(ylim = c(0,1)) +
  theme_ipsum() + 
  theme(legend.position="bottom",
        axis.text.x = element_text(angle=-60),
        panel.spacing.x = unit(0.1, 'mm')) +
  facet_nested(cols=vars(platform, db, workflow),
               rows=vars(rank),
               # ncol=5, 
               scales = 'free') +
  ggtitle('d')

cowplot::plot_grid(ggr2_D1, ggr2_D2, ncol = 2)


saveRDS(ggr2_D2, 'Fig6D.ggr2.rds')

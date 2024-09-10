

res.dir      <- 'HiFi_minitax.vs.sourmash_KeepUnclass'
correl.keepunclass <- fread(paste0(res.dir, '/Observed-Theoretical Correlations.tsv'))
correl.keepunclass$Keep_Unclassified <- T


res.dir      <- 'HiFi_minitax.vs.sourmash_exclUnclass'
correl.exclunclass <- fread(paste0(res.dir, '/Observed-Theoretical Correlations.tsv'))
correl.exclunclass$Keep_Unclassified <- F


plot_sum   <- rbind(correl.exclunclass, correl.keepunclass)
plot_sum$rank <- factor(plot_sum$rank, levels = ranks)

plot_means <- plot_sum %>% 
  group_by(across(any_of(c(colnames(metadata), "rank", 'Keep_Unclassified')))) %>% 
  summarise(R2_mean=mean(R2), R2_sd=sd(R2))

ggr2 <- ggplot(plot_sum, aes(x = Keep_Unclassified, y = R2, fill = workflow)) +
  geom_col(position = position_dodge2()) +
  #geom_errorbar(aes(ymin=R2_mean - R2_sd, ymax=R2_mean+R2_sd)) +
  #geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  geom_text(data=plot_means,
                  aes(x = Keep_Unclassified, y = R2_mean, fill = workflow,
                      label = round(R2_mean, 3)), # paste0('R2=',
                  #position = position_dodge2(),
                  #nudge_x = 0.02, 
                  nudge_y = 0.06, 
                  size = 4) +
  guides(color = guide_legend(direction = "vertical", ncol = 1)) +
  scale_fill_manual(values=as.character(pal.man[c(7, 14)])) +  # c(10:16)
  #facet_wrap(~ sample, scales = "free") +
  #coord_cartesian(xlim = c(0,1), ylim = c(0,1)) +
  theme_ipsum() + 
  theme(legend.position="right",
        axis.text.x = element_text(angle=-60),
        panel.spacing.x = unit(0.1, 'mm')) +
  facet_nested(cols=vars(Keep_Unclassified, workflow),
               rows=vars(rank),
               # ncol=5, 
               scales = 'free')


ggsave('Fig_6C_REVISION.jpg', height = 12, width = 8)

saveRDS(ggr2, 'Fig_6C_REVISION.rds')

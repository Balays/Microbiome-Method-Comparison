
##
res.dir <- 'E:/data/Portik_etal_2022/HiFi_minitax.vs.sourmash'

plot_sum_Portik   <- read_delim(paste0(res.dir, '/Observed-Theoretical Correlations.tsv'))
plot_means_Portik <- read_delim(paste0(res.dir, '/Observed-Theoretical Correlations means.tsv'))


##
res.dir <- 'E:/data/CAMISIM/CAMISIM_MOUSEGUT/minitax_CAMISIM_MG_Illumina_WGS_all_NCBI_genomes'

plot_sum_CAMI   <- read_delim(paste0(res.dir, '/Observed-Theoretical Correlations.tsv'))
plot_means_CAMI <- read_delim(paste0(res.dir, '/Observed-Theoretical Correlations means.tsv'))


##
res.dir <- 'D:/Gemini/Microbacontrol/MCM_V1-V2_Zymo_3600_all_comp'

plot_sum_Ill   <- read_delim(paste0(res.dir, '/Observed-Theoretical Correlations.tsv'))
plot_means_Ill <- read_delim(paste0(res.dir, '/Observed-Theoretical Correlations means.tsv'))


##
res.dir <- 'E:/data/Gemini/Minion/MCM_V1-V9_Zymo_3600_all_comp'

plot_sum_ONT   <- read_delim(paste0(res.dir, '/Observed-Theoretical Correlations.tsv'))
plot_means_ONT <- read_delim(paste0(res.dir, '/Observed-Theoretical Correlations means.tsv'))


##
plot_sum   <- plyr::rbind.fill( plot_sum_ONT,   plot_sum_Ill)
plot_means <- plyr::rbind.fill( plot_means_ONT, plot_means_Ill)


setDT(plot_sum)
setDT(plot_means)


## filter
plot_sum   <- plot_sum[workflow != 'minitax_v3']
plot_means <- plot_sum[workflow != 'minitax_v3']
#plot_sum[workflow == 'minitax_v3', workflow:='minitax']
#plot_means[workflow == 'minitax_v3', workflow:='minitax']

plot_sum[platform == 'illumina', platform:='Illumina']
plot_means[platform == 'illumina', platform:='Illumina']


plot_sum$rank <- factor(plot_sum$rank, levels = ranks)

##
ggr2A <- ggplot(plot_sum[platform == 'Illumina'],
                aes(x = DNA_isolation_method, y = R2, fill = DNA_isolation_method)) +
  geom_col(position = position_dodge2()) +
  #geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  #geom_text_repel(data=plot_means,
  #                aes(y = R2_mean, label = paste0('mean r2=', round(R2_mean, 3))),
  #                position = position_dodge2(),
  #                #nudge_x = 0.02, nudge_y = -0.02,
  #                size = 3) +
  guides(color = guide_legend(direction = "vertical", ncol = 1)) +
  scale_fill_manual(values=as.character(pal.man)) +  # c(10:16)
  #facet_wrap(~ sample, scales = "free") +
  coord_cartesian(ylim = c(0,1)) +
  theme_ipsum() +
  theme(legend.position="none",
        #axis.text.x = element_text(angle=-60),
        panel.spacing.x = unit(0.1, 'mm')) +
  facet_nested(cols=vars(platform, Vregion, db, workflow),
               rows=vars(rank),
               # ncol=5,
               scales = 'free') +
  ggtitle('A')


ggr2B <- ggplot(plot_sum[platform == 'ONT'],
                aes(x = DNA_isolation_method, y = R2, fill = DNA_isolation_method)) +
  geom_col(position = position_dodge2()) +
  #geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  #geom_text_repel(data=plot_means,
  #                aes(y = R2_mean, label = paste0('mean r2=', round(R2_mean, 3))),
  #                position = position_dodge2(),
  #                #nudge_x = 0.02, nudge_y = -0.02,
  #                size = 3) +
  guides(color = guide_legend(direction = "vertical", ncol = 1)) +
  scale_fill_manual(values=as.character(pal.man)) +  # c(10:16)
  #facet_wrap(~ sample, scales = "free") +
  coord_cartesian(ylim = c(0,1)) +
  theme_ipsum() +
  theme(legend.position="none",
        #axis.text.x = element_text(angle=-60),
        panel.spacing.x = unit(0.1, 'mm')) +
  facet_nested(cols=vars(platform, Vregion, db, workflow),
               rows=vars(rank),
               # ncol=5,
               scales = 'free') +
  ggtitle('B')




##
plot_sum_CAMI <- plot_means_CAMI
colnames(plot_sum_CAMI)[5] <- 'R2'
plot_sum_CAMI$project <- 'CAMISIM'

plot_means2   <- plyr::rbind.fill(plot_means_Portik, plot_means_CAMI)
plot_sum2     <- plyr::rbind.fill(plot_sum_Portik,   plot_sum_CAMI)
setDT(plot_means2)
setDT(plot_sum2)

## filter
plot_sum2   <- plot_sum2[rank != 'family']
plot_means2 <- plot_sum2[rank != 'family']

plot_sum2$rank <- factor(plot_sum2$rank, levels = ranks)

ggr2C <- ggplot(plot_sum2[plot_sum2$project == 'Portik_etal_2022'],
                aes(x = workflow, y = R2, fill = workflow)) +
  geom_col(position = position_dodge2()) +
  #geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  #geom_text_repel(data=plot_means2,
  #                aes(y = R2_mean, label = paste0('mean r2=', round(R2_mean, 3))),
  #                position = position_dodge2(),
  #                #nudge_x = 0.02, nudge_y = -0.02,
  #                size = 3) +
  guides(color = guide_legend(direction = "vertical", ncol = 1)) +
  scale_fill_manual(values=as.character(pal.man)) +  # c(10:16)
  #facet_wrap(~ sample, scales = "free") +
  coord_cartesian(ylim = c(0,1)) +
  theme_ipsum() +
  theme(legend.position="none",
        #axis.text.x = element_text(angle=-60),
        panel.spacing.x = unit(0.1, 'mm')) +
  facet_nested(cols=vars(platform, db, workflow),
               rows=vars(rank),
               # ncol=5,
               scales = 'free') +
  ggtitle('C')


ggr2D <- ggplot(plot_sum2[plot_sum2$project == 'CAMISIM'],
                aes(x = workflow, y = R2, fill = workflow)) +
  geom_col(position = position_dodge2()) +
  #geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  #geom_text_repel(data=plot_means2,
  #                aes(y = R2_mean, label = paste0('mean r2=', round(R2_mean, 3))),
  #                position = position_dodge2(),
  #                #nudge_x = 0.02, nudge_y = -0.02,
  #                size = 3) +
  guides(color = guide_legend(direction = "vertical", ncol = 1)) +
  scale_fill_manual(values=as.character(pal.man)) +  # c(10:16)
  #facet_wrap(~ sample, scales = "free") +
  coord_cartesian(ylim = c(0,1)) +
  theme_ipsum() +
  theme(legend.position="none",
        #axis.text.x = element_text(angle=-60),
        panel.spacing.x = unit(0.1, 'mm')) +
  facet_nested(cols=vars(platform, db, workflow),
               rows=vars(rank),
               # ncol=5,
               scales = 'free') +
  ggtitle('D')



#cowplot::plot_grid(ggr2, cowplot::plot_grid(ggr2.2, NULL, rel_widths = c(3,4),nrow=1), ncol = 1, rel_heights = c(1, 0.8))


ggr2 <- cowplot::plot_grid(ggr2A,
                           ggr2C,
                           ggr2B,
                           ggr2D,
                           #cowplot::plot_grid(ggr2.2, NULL, rel_widths = c(3,4),nrow=1),
                           ncol = 2
                           , rel_widths = c(1, 0.5)
                           )

ggr2
ggsave('E:/data/Gemini/MCM_R2.means.jpg', ggr2, width=16, height = 18)






rm(list = ls())
load("MCM_D6300_dada2_EMUdb.RData")
d6300_dada2_emudb    <- ps.all

#rm(list = grep('d6300_dada2_emudb', ls(), invert = T, value = T))
load("MCM_D6300_dada2_SILVA.RData")
d6300_dada2_silva    <- ps.all

#rm(list = grep('d6300_dada2_silva', ls(), invert = T, value = T))
load("MCM_D6300_EMU_EMUdb.RData")
d6300_emu_emudb      <- ps.all

#rm(list = grep('d6300_emu_emudb', ls(), invert = T, value = T))
load("MCM_D6300_minitax_EMUdb.RData")
d6300_minitax_emu    <- ps.all

#rm(list = grep('d6300_minitax_emu', ls(), invert = T, value = T))
load("MCM_D6300_minitax_NCBI.RData")
d6300_minitax_ncbi   <- ps.all

#rm(list = grep('d6300_minitax_ncbi', ls(), invert = T, value = T))
load("MCM_D6331_Emu_EMUdb.RData")
d6331_emu_emudb      <- ps.all

#rm(list = grep('d6331_emu_emudb', ls(), invert = T, value = T))
load("MCM_D6331_minitax_EMUdb.RData")
d6331_minitax_emudb  <- ps.all

#rm(list = grep('d6331_minitax_emudb', ls(), invert = T, value = T))
load("MCM_D6331_minitax_NCBI.RData")
d6331_minitax_ncbi   <- ps.all

d.all_minitax_ncbi   <- merge_phyloseq(d6300_minitax_ncbi,
                                       d6331_minitax_ncbi)

d.all_minitax_emudb   <- merge_phyloseq(d6300_minitax_emu,
                                        d6331_minitax_emudb)

d.all_emu_emudb   <- merge_phyloseq(d6300_emu_emudb,
                                    d6331_emu_emudb)

d.all_dada2   <- merge_phyloseq(d6300_dada2_emudb, 
                                d6300_dada2_silva)


d.all_minitax <- merge_phyloseq(d.all_minitax_ncbi,
                                d.all_minitax_emudb)


ps.all <- merge_phyloseq(d.all_minitax,
                         d.all_emu_emudb,
                         d.all_dada2)

ps <- prune_samples(grep('Zymo_6*', sample_names(ps.all), value = T, invert = T), ps.all)

metadata <- data.frame(sample_data(ps))

#  Alpha diversity estimations


richness <- merge(
  estimate_richness(ps, split = TRUE, measures = c("Shannon", "Simpson")),
  metadata,
  by=0) %>%
  dplyr::select(c(sample, sample_name, DNA_isolation_method, platform, Vregion, GS_version, workflow, db, Shannon, Simpson))

richness$Vregion[richness$Vregion == 'V1-V9'] <- richness$Vregion[richness$Vregion == 'V1_V9']

write_tsv(richness, 'MCM_all_Adiversity_estimates.tsv')



### DNA_isolation_method vs Library

measure <- 'Shannon'
#psh     <- plot_richness(ps, x="Vregion", measures=measure, color="DNA_isolation_method")

gg.shannon <- ggplot(richness, #[psh$data$workflow == 'minitax', ],
                     aes(x=GS_version, y=Shannon)) +
  geom_violin(aes(fill=GS_version), position = 'dodge', alpha=0.6) +
  geom_jitter(aes(color=GS_version), position = 'dodge') +
  scale_y_continuous(name='Shannon index') +
  scale_fill_manual(values=as.character(pal.man)) +
  scale_color_manual(values=as.character(pal.man)) +
  theme_ipsum() +
  theme(axis.title.x  = element_blank(),
        axis.text.x   = element_text(angle=-45, hjust=0),
        panel.spacing = unit(2, 'mm'), 
        legend.position = 'none') +
  facet_nested(cols=vars(workflow, db, DNA_isolation_method, GS_version ),
               rows=vars(Vregion),
               scales='free_x')

gg.shannon
  
  
  
  
cowplot::plot_grid(
  ,
  
  ggplot(psh$data, #[psh$data$workflow == 'minitax', ],
         aes(x=DNA_isolation_method, y=value)) +
    geom_violin(aes(fill=DNA_isolation_method), position = 'dodge', alpha=0.6) +
    geom_jitter(aes(color=DNA_isolation_method)) +
    scale_fill_manual(values=as.character(pal.man)) +
    scale_color_manual(values=as.character(pal.man)) +
    scale_y_continuous(name='Shannon index') +
    theme_ipsum() + 
    theme(axis.title.x  = element_blank(),
          axis.text.x = element_text(angle=-45, hjust=0)) + 
    facet_nested(cols=vars(workflow), scales='free_x')
  ,
  nrow=1, rel_widths = c(3,2)
  
  )




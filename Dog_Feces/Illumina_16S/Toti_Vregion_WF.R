rm(list = ls(all=T))
gc()
library(BiocParallel)
library(ggplot2)
library(dada2)
library(tidyr)
library(dplyr)
library(Biostrings)

#### Metadata, settings and parameters ####

### working directory
setwd('/mnt/f/data/Gemini/Toti_Vregion_comp3')

### output of dada2
par.dir <- '/mnt/f/data/Gemini/Toti_Vregion_comp3/dada2'
dir.create(par.dir)

### metadata and parent fastq directory
indir   <- '/mnt/f/data/Gemini/Toti_Vregion'

source('/mnt/f/data/Gemini/make.metadata.Toti.R')
#metadata$sample <- metadata$sample_ID
rownames(metadata) <- metadata$sample #_ID
#write.table(metadata, "/mnt/f/data/Gemini/Toti_metadata.tsv", row.names = F, quote = F, sep='\t')

### filtering version
filt.version <- 'filt.v2'

### config
config <- read.delim( paste0('/mnt/f/data/Gemini/Toti_Vregion_comp3/','Toti_config.', filt.version, '.tsv'))
config$sequencing_ID   <- paste(sep='_', config$sequencing_date, config$host, config$Vregion)
metadata <- merge(metadata, config, by=c('sequencing_date', 'sequencing_ID', 'host', 'Vregion' ))

### Sequence table for combined bimera removing
seqtab.all    <- data.frame(NULL)

### Remove previous results?
rm.prev.results <- F
#### ####
## 

#### carry out dada2 filtering, denoising and mergeing #### 
seq_ids <- unique(metadata$sequencing_ID)
n       <- c(5:length(seq_ids)) #  7,
n
for (i in seq_ids[n]) {
  
  ## fastq dir
  fastq.dir <- paste0(indir, '/', i)
  
  ## fastq files
  stopifnot(all.equal(
    list.files(fastq.dir, pattern = '.fastq$', full.names = T, recursive = T),
    metadata$fastq[metadata$sequencing_ID == i ]
  ))  
  fastq   <- metadata$fastq[metadata$sequencing_ID == i ]
  
  print(nrow( metadata[metadata$sequencing_ID == i, ] ))
  
  ## results directory for each sequencing (delete if previously created!!!)
  setwd(par.dir)
  res.dir   <- paste0(par.dir, '/', i) 
  
  if(rm.prev.results) {file.remove(res.dir); unlink(res.dir, recursive = TRUE) }
  dir.create(res.dir) 
  setwd(res.dir)
  
  ## sequencing config
  config.seq <- config[config$sequencing_ID == i, ]
  
  ## running dada2
  message('running dada2 on: ', i, '...')
  source('/mnt/f/data/Gemini/Toti_Vregion_comp2/Toti_Vregion_dada2.linux.R')
  
  print(paste0(i, ' ;', filt.version, " done!"))
  
  ## Sequence tables
  seqtab      <- as.matrix(t(seqtab.df))
  seqtab.df   <- data.frame(sequence=rownames(seqtab.df), seqtab.df)
  seqtab.gt   <- gather(seqtab.df, sample, count, -1)
  seqtab.all  <- plyr::rbind.fill(seqtab.gt, seqtab.all)
  
  ## Track tables
  
  track.all  <- plyr::rbind.fill(track, track.all)
  gc()
}

setwd(par.dir)
#### ####
## 

#### Remove bimeras after merging of samples ####
## run this in linux
filt.version <- 'filt.v2'
seqtab.all   <- read.delim(paste0(filt.version, '_all_sequence.table.tsv'))
seqtab.merge <- spread(seqtab.all, sample, count, fill=0)
rownames(seqtab.merge) <- seqtab.merge$sequence

vregions <- c('v1_3', 'v1_2', 'v3_4')

vr <- vregions[1]
seqtab.vr         <- seqtab.merge[,grep(vr, colnames(seqtab.merge))]
seqtab.vr         <- seqtab.vr[rowSums(seqtab.vr, na.rm = T) > 0, ]
write.table(seqtab.vr, paste0(vr, '_', filt.version, '_all_sequence.table.tsv'), sep="\t", row.names=T, quote=F)

vr <- vregions[2]
seqtab.vr         <- seqtab.merge[,grep(vr, colnames(seqtab.merge))]
seqtab.vr         <- seqtab.vr[rowSums(seqtab.vr, na.rm = T) > 0, ]
write.table(seqtab.vr, paste0(vr, '_', filt.version, '_all_sequence.table.tsv'), sep="\t", row.names=T, quote=F)

vr <- vregions[3]
seqtab.vr         <- seqtab.merge[,grep(vr, colnames(seqtab.merge))]
seqtab.vr         <- seqtab.vr[rowSums(seqtab.vr, na.rm = T) > 0, ]
write.table(seqtab.vr, paste0(vr, '_', filt.version, '_all_sequence.table.tsv'), sep="\t", row.names=T, quote=F)


seqtab.all.nochim.gt <- data.frame(NULL) 

for (vr in vregions) {
  message('importing ', vr, '...')
  seqtab.vr.df            <- read.delim(paste0(vr, '_', filt.version, '_all_sequence.table.tsv'))
  seqtab.vr               <- as.matrix(t(seqtab.vr.df))
  seqtab.vr.nochim        <- removeBimeraDenovo(seqtab.vr, method="consensus", multithread=multithread, verbose=TRUE)
  seqtab.vr.nochim.df     <- as.data.frame(t(seqtab.vr.nochim))
  seqtab.vr.nochim.df     <- data.frame(sequence=row.names(seqtab.vr.nochim.df), seqtab.vr.nochim.df)
  seqtab.vr.nochim.gt     <- gather(seqtab.vr.nochim.df, sample, count, -1)
  seqtab.all.nochim.gt    <- plyr::rbind.fill(seqtab.vr.nochim.gt, seqtab.all.nochim.gt)
}

seqtab.merge.nochim.df  <- spread(seqtab.all.nochim.gt, sample, count, fill=0)
rownames(seqtab.merge.nochim.df) <- seqtab.merge.nochim.df$sequence
seqtab.merge.nochim.df  <- seqtab.merge.nochim.df[,-1]
seqtab.merge.nochim     <- as.matrix(t(seqtab.merge.nochim.df))

write.table(seqtab.merge.nochim.df, paste0(filt.version, '_all_seqtab.nochim.tsv'), sep="\t", row.names=T, quote=F)
#### ####
## 

#### Taxonomic annotation ####
multithread <- T
seqtab.merge.nochim.df  <- read.delim(paste0(filt.version, '_all_seqtab.nochim.tsv'))
seqtab.merge.nochim     <- as.matrix(t(seqtab.merge.nochim.df))

### Databases 
dbfile <- c(
  "/mnt/f/data/databases/silva_nr99_v138.1_train_set.fa.gz",
  "/mnt/f/data/databases/silva_species_assignment_v138.1.fa.gz",
  "/mnt/f/data/databases/silva_nr99_v138.1_wSpecies_train_set.fa.gz" )

## OR for EMU db
dbfile <- rep("/mnt/f/data/databases/EMU/species_taxid.dada2.formatted.fasta", 3)

message('Assigning taxonomy in one step...')
spec_pool_silva      <- assignTaxonomy(seqtab.merge.nochim, 
                                       dbfile[3], 
                                       multithread=multithread)

write.table(spec_pool_silva, paste0(filt.version, '_all_spec_pool_EMU.tsv'), sep="\t", row.names=T, quote=F)
#write.table(spec_pool_silva, paste0(filt.version, '_all_spec_pool_silva.tsv'), sep="\t", row.names=T, quote=F)

#### ####
## 






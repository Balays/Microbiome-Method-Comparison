rm(list = ls(all=T))
gc()
library(tidyr)
library(dplyr)
library(ggplot2)
library(dada2)
library(Biostrings)
library(phyloseq)
#library(misc)
#library(gff)

##### source settings
setwd('F:/data/Gemini/Toti_Vregion_comp3')
source('settings.R')
vregions <- c('v1_3', 'v1_2', 'v3_4')
db <- 'EMU' ## OR db <- 'SILVA'

seqtab.merge.nochim.df  <- read.delim(paste0(filt.version, '_all_seqtab.nochim.tsv'))
seqtab.merge.nochim     <- as.matrix(t(seqtab.merge.nochim.df))

spec_pool_silva         <- read.delim(paste0(filt.version, '_all_spec_pool_', db, '.tsv'))

### CHECK
stopifnot(rownames(spec_pool_silva) == colnames(seqtab.merge.nochim))

### Tax table
##  Taxonomic annotation reformat
tax.annot <- spec_pool_silva
tax.annot <- data.frame(tax.annot, stringsAsFactors = F)

## CAT genus and species
#tax.annot[!is.na(tax.annot$Species) , 7] <- paste0(tax.annot[!is.na(tax.annot$Species), "Genus"], " ", tax.annot[!is.na(tax.annot$Species), "Species"])

## fill in NA spec_gen 
#tax.annot$Spec_Gen <-last.non.na(tax.annot[,1:(t-1)])$x
#tax.annot$Spec_Gen[is.na(tax.annot$Species)]  <- paste0('Unknown ', tax.annot$Spec_Gen[is.na(tax.annot$Species)] )

## order tax table
#tax.annot <- tax.annot[order(tax.annot$Spec_Gen), ]
tax.annot <- tax.annot[order(tax.annot$Species), ]
tax.annot$ASV_id <- paste0('ASV_', 1:nrow(tax.annot))
##

### Sample data table
metafilt <- metadata
metafilt$fastq <- gsub('.*/', '', metadata$fastq)
metafilt$fastq <- paste0(filt.version, '_', metafilt$fastq)


### subsetting data
metafilt <- metafilt[is.element(metafilt$fastq, colnames(seqtab.merge.nochim.df)), ]
metafilt$sample <- gsub('_R1', '', metafilt$sample)
rownames(metafilt) <- metafilt$sample

xdfa <- data.frame(seqtab=colnames(seqtab.merge.nochim.df), x="seqtab")
xdfb <- data.frame(metaf =metafilt$fastq, x='metaf')
xdf  <- merge(xdfa, xdfb, all=T, by=1)
message('These samples were filtered out from "seqtab.merge.nochim" : \n', paste0(xdf[is.na(xdf$x.y),1], collapse = '\n'))
seqtab.merge.nochim.df <- seqtab.merge.nochim.df[,metafilt$fastq]

### replacing sample names
colnames(seqtab.merge.nochim.df) <- metafilt$sample
seqtab.merge.nochim    <- as.matrix(t(seqtab.merge.nochim.df))

### OTU Table (counts)
all.asv.df <- merge(tax.annot, seqtab.merge.nochim.df, by=0)
colnames(all.asv.df)[1] <- 'sequence'
seqtab.merge.nochim.df <- all.asv.df[,-c(1:(ncol(tax.annot)+1))]
rownames(seqtab.merge.nochim.df) <- all.asv.df$sequence

### check
stopifnot(all.asv.df$sequence == colnames(seqtab.merge.nochim))

### Builid phyloseq object
taxtab   <- tax_table(as.matrix(tax.annot))
otutab   <- otu_table(as.matrix(seqtab.merge.nochim.df), taxa_are_rows = T)
samp_dat <- sample_data(metafilt)

ps <- phyloseq(otutab,
               samp_dat, 
               taxtab
               # OR if there's a phylogenetic tree: 
               #, phy_tree(fitGTR$tree)
)

### Sequences
dna        <- DNAStringSet(getSequences(seqtab.merge.nochim))
names(dna) <- all.asv.df$ASV_id

### Read counts


gc()
image.name   <- paste0('Toti_Vregion_dada2_', filt.version, '_', db, '.RData')
save.image(image.name)

stop()

## Phylogenetic tree?
if (phytree) {
  source('phytree.R')
  
  ps <- phyloseq(otu_table(s.dada2, taxa_are_rows=F), ## d.n0.clr ## s.dada2  ## s.n0 ## or d.n0
                 sample_data(metafilt), 
                 tax_table(tax.annot)
                 # OR if there's a phylogenetic tree: 
                 , phy_tree(fitGTR$tree)
  )
} else {
  ps <- phyloseq(otu_table(s.dada2, taxa_are_rows=F), ## d.n0.clr ## s.dada2  ## s.n0 ## or d.n0
                 sample_data(metafilt), 
                 tax_table(tax.annot)
                 # OR if there's a phylogenetic tree: 
                 #, phy_tree(fitGTR$tree)
  )
}

## Prevalence filtering?
source('prev.filter.R')

save.image(image.name)


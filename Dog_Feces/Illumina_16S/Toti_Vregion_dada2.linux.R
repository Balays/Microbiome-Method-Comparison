
### This script should be run in linux!
#mypaths <- c("/mnt/d/R/site-library", .libPaths())
#.libPaths(mypaths)
getN <- function(x) sum(getUniques(x))

#### Settings #### 
multithread <- T

plot.qc <- T

dir.create(paste0(res.dir, '/', filt.version))

merge.all <- T
justConcatenate <- config.seq$justConcatenate

#### Importing fastq file #### 
Fs <- list.files(fastq.dir, pattern="_R1.fastq$", full.names=T)
Rs <- list.files(fastq.dir, pattern="_R2.fastq$", full.names=T)

## Remove empty files
Fs <- Fs[file.size(Fs) != 0 ]
Rs <- Rs[file.size(Rs) != 0 ]

Fnames <- basename(gsub("_R1.fastq","", Fs))
Rnames <- basename(gsub("_R2.fastq","", Rs))

if(!all.equal(Fnames,Rnames)) {stop('sample names in forward and reverse reads are different!!')}


## Plot QC
if (plot.qc) {
  qc.Fs <- plotQualityProfile(Fs); ggsave(paste0('./QC.frw.png'), qc.Fs, height = 16, width = 20)
  qc.Rs <- plotQualityProfile(Rs); ggsave(paste0('./QC.rev.png'), qc.Rs, height = 16, width = 20)
}


#### Filtering ####

## Trimming primers
Fprimer <- 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCA'  # readLines("./fprimer.txt")
Rprimer <- 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT'  # readLines("./rprimer.txt")

Flen <- config.seq$Flen ## 9 #nchar(Fprimer) ## 
Rlen <- config.seq$Rlen ## 9 #nchar(Rprimer) ## 

filtFs <- paste0(filt.version, "_", Fnames, "_R1.fastq")
filtRs <- paste0(filt.version, "_", Rnames, "_R2.fastq")

message("filtering reads...")
out <- filterAndTrim(Fs, filtFs, 
                     Rs, filtRs, 
                     matchIDs=T, 
                     trimLeft = c(Flen, Rlen), 
                     truncLen = config.seq$truncLen, #c(245, 245), 
                     maxEE = config.seq$maxEE, #c(2,2), 
                     maxN=0, truncQ = 2, rm.phix = T, 
                     compress = F, multithread = multithread)
#c(220,160), 'Inf', c(2,2)

## Plot QC filtered
if (plot.qc) {
  qc.filt.Fs <- plotQualityProfile(filtFs); 
  ggsave(paste0('QC.frw.', filt.version, '.png'), qc.filt.Fs, height = 16, width = 20)
  qc.filt.Rs <- plotQualityProfile(filtRs); 
  ggsave(paste0('QC.rev.', filt.version, '.png'), qc.filt.Rs, height = 16, width = 20)
}

#### Generating track-file ####
track <- as.data.frame(cbind(Fnames, out))
#write.table(track, paste0("track.", filt.version, ".tsv"), sep="\t", row.names=T, quote=F)

#### learning errors ####
message("learning errors...")
errF <- learnErrors(filtFs, multithread=multithread)
errR <- learnErrors(filtRs, multithread=multithread)

err.Fs <- plotErrors(errF, nominalQ=multithread); 
ggsave(paste0('Err.frw.', filt.version, '.png'), err.Fs, height = 16, width = 20)
err.Rs <- plotErrors(errR, nominalQ=multithread); 
ggsave(paste0('Err.rev.', filt.version, '.png'), err.Rs, height = 16, width = 20)

#### dereplicating ####
message("dereplicating... ")
derepFs <- derepFastq(filtFs)
derepRs <- derepFastq(filtRs)

#### running dada2 ####
message("running dada2... ")
pool_dadaFs <- dada(derepFs, err=errF, multithread=multithread, pool=T)
pool_dadaRs <- dada(derepRs, err=errR, multithread=multithread, pool=T)

#### merging pairs, removing bimeras and assigning taxonomy ####
message("merging pairs... ")
mergers <- mergePairs(pool_dadaFs, derepFs, 
                      pool_dadaRs, derepRs, 
                      verbose=TRUE, justConcatenate = justConcatenate)
seqtab    <- makeSequenceTable(mergers)
seqtab.df <- as.data.frame(t(seqtab))

message('writing sequence tables...')
write.table(seqtab.df,        paste0(filt.version, '_sequence.table.tsv'), sep="\t", row.names=T, quote=F)
message('writing merged sequences in .fasta format...')

seq        <- DNAStringSet(getSequences(seqtab))
names(seq) <- paste0('ASV_', 1:nrow(seqtab.df))
writeXStringSet(seq, paste0(filt.version, '_sequences.fasta'), append = FALSE)

track <- as.data.frame(cbind(Fnames, out, sapply(pool_dadaFs, getN), sapply(pool_dadaRs, getN), sapply(mergers, getN)
                             #, rowSums(seqtab.nochim)
                             ))
colnames(track) <- c("sample_name", "input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- Fnames

  
message("writing image...")
save.image(paste0("dada2.", filt.version, ".RData"))

rm(derepFs, derepRs, pool_dadaFs, pool_dadaRs, qc.RS, qc.Fs, mergers, err.Fs, err.Rs)
image.name <- paste0("ASVs.dada2.", filt.version, ".RData")
save.image(image.name)






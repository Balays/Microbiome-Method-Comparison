
dirs <- c(
  grep('MCM_D6331', list.dirs('.', full.names = T, recursive = F), value = T ),
  grep('MCM_D6300', list.dirs('.', full.names = T, recursive = F), value = T )
)

all.sig.freq <- data.table(NULL)

for(dir in dirs) {
  
  try({
    sig.freq     <- fread(paste0(dir, '/sig.frq.tsv'))
    all.sig.freq <- rbind(sig.freq, all.sig.freq)
  })
  
}

all.sig.freq

all.sig.freq.sp <- dcast(all.sig.freq, sample + workflow + db  ~ is.significant, value.var='N')


all.l2f.data <- data.table(NULL)

for(dir in dirs) {
  
  try({
    l2f.data     <- fread(paste0(dir, '/l2F_diff.dt.tsv'))
    all.l2f.data <- rbind(l2f.data, all.l2f.data)
  })
}

all.l2f.data[,.N,by=.(sample, db, workflow, is.significant)]

#all.l2f.data[,DNA_isolation_method ]

fwrite(all.l2f.data, 'all.l2f.data.tsv', sep='\t')




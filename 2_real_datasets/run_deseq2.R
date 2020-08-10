 arguments <- commandArgs(trailingOnly = TRUE)
if(length(arguments) < 1 | length(arguments) > 3){
  stop('Uses: Rscript.R run_deseq2.R FunctionalProfile groupFile outputFile')
}else if(length(arguments) == 3){
  print(getwd())
  met <- read.table(file = gzfile(arguments[1]),header = T,sep = '\t',row.names = 1,check.names = F)
  coldata <- read.table(file = arguments[2],header = T,sep = '\t',row.names = 1,check.names = F)
  c <- intersect(colnames(met),rownames(coldata))
  coldata <- as.data.frame(coldata[c,])
  colnames(coldata) <- c('group')
  library('DESeq2')
  met <- met[,c]
  met <- round(met)
  #deseq2
  dds <- DESeqDataSetFromMatrix(countData = met,
                                colData = coldata,
                                design= ~ group)
  gm_mean = function(x, na.rm=TRUE){
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
  geoMeans = apply(counts(dds), 1, gm_mean)
  dds = estimateSizeFactors(dds, geoMeans = geoMeans)
  dds = DESeq(dds, fitType="local")
  res <- results(dds)
  res = res[order(res$padj, na.last=NA), ]
  alpha = 0.05
  sigtab = res[(res$padj < alpha), ]
  write.table(x = as.data.frame(res),file = arguments[3],quote = F,sep = '\t',col.names = NA)
}

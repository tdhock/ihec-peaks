library(data.table)
library(PeakError)
library(namedCapture)
library(ggplot2)

pre <- "http://hubs.hpc.mcgill.ca/~thocking/PeakSegFPOP-labels/H3K36me3_TDH_immune/"
read <- function(d){
  f <- sprintf("peaks_matrix_%s.tsv.gz", d)
  if(!file.exists(f)){
    u <- paste0(pre, f)
    download.file(u, f)
  }
  dt <- fread(paste("zcat", f))
  mat <- as.matrix(dt[,-1,with=FALSE])
  rownames(mat) <- dt$peak.name
  mat
}
group.mat <- read("group")
lik.mat <- read("likelihood")

regions.RData.vec <- Sys.glob("~/projects/chip-seq-paper/chunks/H3K36me3_TDH_other/*/regions.RData")
regions.dt.list <- list()
for(regions.RData in regions.RData.vec){
  load(regions.RData)
  regions.dt.list[[regions.RData]] <- data.table(regions)
}
regions.dt <- do.call(rbind, regions.dt.list)
regions.dt[, sample.path := paste0(cell.type, "/", sample.id)]
setkey(regions.dt, sample.path)

pattern <- paste0(
  "(?<chrom>chr[^:]+)",
  ":",
  "(?<peakStart>[0-9]+)",
  "-",
  "(?<peakEnd>[0-9]+)")
peak.pos.dt <- data.table(str_match_named(rownames(lik.mat), pattern, list(
  peakStart=as.integer,
  peakEnd=as.integer)))

diff.dt.list <- list()
total.dt.list <- list()
for(sample.path in colnames(lik.mat)){
  sample.regions <- regions.dt[sample.path, nomatch=0L]
  if(nrow(sample.regions)){
    lik.dt <- data.table(peak.pos.dt, lik=lik.mat[, sample.path])
    setkey(lik.dt, chrom, peakStart, peakEnd)
    setkey(sample.regions, chrom, chromStart, chromEnd)
    over.dt <- foverlaps(sample.regions, lik.dt, nomatch=0L)
    peak.dt <- unique(over.dt[, list(
      chrom, chromStart=peakStart, chromEnd=peakEnd, lik)])[order(-lik)]
    error.dt.list <- list()
    for(row.i in 0:nrow(peak.dt)){
      if(row.i==0){
        peaks <- Peaks()
        thresh <- Inf
      }else{
        peaks <- peak.dt[1:row.i]
        thresh <- peak.dt$lik[row.i]
      }
      error.df <- PeakError(peaks, sample.regions)
      error.dt.list[[paste(row.i)]] <- with(error.df, data.table(
        thresh, ##include peak if lik(peak) >= thresh.
        tp=sum(tp),
        fp=sum(fp)
      ))
    }
    error.dt <- do.call(rbind, error.dt.list)
    diff.dt.list[[sample.path]] <- error.dt[, data.table(
      sample.path,
      thresh=thresh[-1], fp=diff(fp), tp=diff(tp))]
    total.dt.list[[sample.path]] <- with(error.df, data.table(
      sample.path,
      possible.tp=sum(possible.tp),
      possible.fp=sum(possible.fp)))
  }
}
diff.dt <- do.call(rbind, diff.dt.list)[order(-thresh)]
total.dt <- do.call(rbind, total.dt.list)

possible.fp <- sum(total.dt$possible.fp)
possible.tp <- sum(total.dt$possible.tp)
roc.dt <- diff.dt[, data.table(
  thresh=c(Inf, thresh, 0),
  fp=c(0, cumsum(fp), possible.fp),
  tp=c(0, cumsum(tp), possible.tp)
  )]
roc.dt[, FPR := fp/possible.fp]
roc.dt[, TPR := tp/possible.tp]
WeightedROC::WeightedAUC(roc.dt[order(thresh)])
ggplot()+
  geom_path(aes(FPR, TPR), data=roc.dt)+
  coord_equal(xlim=c(0,1),ylim=c(0,1))

if(!file.exists("peaks_summary.tsv")){
  download.file("http://hubs.hpc.mcgill.ca/~thocking/PeakSegFPOP-labels/H3K36me3_TDH_immune/peaks_summary.tsv", "peaks_summary.tsv")
}
summary.dt <- fread("peaks_summary.tsv")
summary.dt[, table(n.samples.up, n.groups.up)]

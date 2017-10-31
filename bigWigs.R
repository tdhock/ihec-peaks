
library(data.table)
bigWigs.dt <- fread("2017-10-bigWigs.txt", header=FALSE, sep="\n")
setnames(bigWigs.dt, "url")

library(namedCapture)
pattern <- paste0(
  "http://epigenomesportal.ca/tracks/",
  "(?<project>[^/]+)",
  "/",
  "(?<genome>[^/]*)",
  "/",
  "(?<id>[0-9]+)",
  "[.]",
  "(?:[^.]+)",
  "[.]",
  "(?<cellType>.+)",
  "[.]",
  "(?<experiment>[^.]+)",
  "[.]",
  "(?<dataType>[^.]+)",
  "[.]",
  "(?<extension>bigWig|bigBed)")
match.mat <- str_match_named(bigWigs.dt$url, pattern)
match.dt <- data.table(match.mat, row=1:nrow(bigWigs.dt))
not.matched <- bigWigs.dt[is.na(match.dt$project)]
stopifnot(nrow(not.matched)==0)

bigWigs.dt[grepl("Duke", url)]

exp.counts <- match.dt[, list(count=.N), by=list(experiment)][order(-count)]
match.dt[experiment=="H3K27ac"]
stopifnot(match.dt[, sum(dataType=="H3K27ac")]==0)
exp.counts[1:10]
dput(exp.counts[grepl("^H", experiment)][1:6, experiment])

##    experiment count
## 1:    H3K27ac  2892
## 2:    H3K4me3  2076
## 3:    H3K4me1  1969
## 4:   H3K27me3  1650
## 5:   H3K36me3  1396
## 6:    H3K9me3  1316

## What is the difference between signal, peaks_bw, and signal_unstranded?
match.dt[experiment=="H3K27ac" & extension=="bigWig", table(dataType)]

## seems like peaks_bw is not coverage.
(peaks_bw <- match.dt[experiment=="H3K27ac" & dataType=="peaks_bw"])
bigWig.i <- peaks_bw$row[1]
bigWig.url <- bigWigs.dt$url[bigWig.i]
system(paste("bigWigToBedGraph", bigWig.url, "/dev/stdout -chrom=chr1 -start=0 -end=1000000"))

## what about signal? seems OK.
(signal <- match.dt[experiment=="H3K27ac" & dataType=="signal"])
bigWig.i <- signal$row[1]
bigWig.url <- bigWigs.dt$url[bigWig.i]
system(paste("bigWigToBedGraph", bigWig.url, "/dev/stdout -chrom=chr1 -start=0 -end=100000"))

## what about signal_unstranded? seems OK.
(signal <- match.dt[experiment=="H3K27ac" & dataType=="signal_unstranded"])
bigWig.i <- signal$row[1]
bigWig.url <- bigWigs.dt$url[bigWig.i]
system(paste("bigWigToBedGraph", bigWig.url, "/dev/stdout -chrom=chr1 -start=0 -end=100000 | tee test.bedGraph"))

core.hist.vec <- c(
  "H3K27ac", "H3K4me3", "H3K4me1", "H3K27me3", "H3K36me3", "H3K9me3")
match.dt[experiment %in% core.hist.vec, table(experiment, dataType)]
core.hist.dt <- match.dt[experiment %in% core.hist.vec & grepl("signal", dataType)]
core.hist.dt[, table(experiment, dataType)]


library(RJSONIO)
library(data.table)
Biosample.vec <- c(
  "skeletal muscle myoblast"="skeletalMuscle",
  "thyroid gland"="thyroid",
  "transverse colon"="colon",
  #"CD4-positive helper T cell", ##Roadmap.
  "thoracic aorta"="aorta")
dir.create("json", showWarnings=FALSE)
getExpList <- function(acc){
  accession.json <- file.path("json", paste0(acc, ".json"))
  if(!file.exists(accession.json)){
    u <- paste0(
      "https://www.encodeproject.org/experiments/",
      acc,
      "/?format=json")
    download.file(u, accession.json)
  }
  fromJSON(accession.json)
}
for(experiment in c("H3K36me3", "H3K4me3" )){
  samples.tsv <- paste0("ENCODE_", experiment, ".tsv")
  samples.dt <- fread(samples.tsv)
  setkey(samples.dt, Biosample)
  some.samples <- samples.dt[names(Biosample.vec)]
  for(sample.i in 1:nrow(some.samples)){
    one.sample <- some.samples[sample.i]
    data.dir <- file.path("labels", paste0(experiment, "_TDH_ENCODE"))
    elist <- getExpList(one.sample$Accession)
    clist <- getExpList(elist$possible_controls[[1]]$accession)
    sample.group.list <- structure(
      list(elist, clist),
      names=paste0(
        Biosample.vec[one.sample$Biosample],
        c("", "_Input")))
    for(sample.group in names(sample.group.list)){
      jlist <- sample.group.list[[sample.group]]
      afile.list <- jlist$files[sapply(
        jlist$files, "[[", "output_type") == "alignments" &
          sapply(jlist$files, "[[", "file_type") == "bam"]
      cat(sprintf("%4d / %4d %s\n", sample.i, nrow(some.samples), sample.group))
      for(afile in afile.list){
        sample.id <- sub(".bam$", "", basename(afile$href))
        alignments.bam <- file.path(
          data.dir, "samples",
          sample.group, sample.id, "alignments.bam")
        already.downloaded <- if(file.exists(alignments.bam)){
          md5.dt <- fread(paste("md5sum", alignments.bam), header=FALSE)
          md5.dt$V1 == afile$md5sum
        }else{
          FALSE
        }
        if(!already.downloaded){
          dir.create(
            dirname(alignments.bam), showWarnings=FALSE, recursive=TRUE)
          u <- paste0("https://www.encodeproject.org", afile$href)
          download.file(u, alignments.bam, method="wget")
        }
      }
    }
  }
}

library(httr)
library(RJSONIO)
library(data.table)
meta.dt <- fread("ENCODE_H3K4me3.tsv")

Biosample.vec <- c(
  "skeletal muscle myoblast"="skeletalMuscle",
  "thyroid gland"="thyroid",
  "transverse colon"="colon",
  #"CD4-positive helper T cell", ##Roadmap.
  "thoracic aorta"="aorta")
setkey(meta.dt, Biosample)
some.samples <- meta.dt[names(Biosample.vec)]

for(sample.i in 1:nrow(some.samples)){
  one.sample <- some.samples[sample.i]
  accession.json <- paste0(one.sample$Accession, ".json")
  if(!file.exists(accession.json)){
    u <- paste0(
      "https://www.encodeproject.org/experiments/",
      one.sample$Accession,
      "/?format=json")
    download.file(u, accession.json)
  }
  jlist <- fromJSON(accession.json)
  afiles <- jlist$files[sapply(
    jlist$files, "[[", "output_type") == "alignments" &
      sapply(jlist$files, "[[", "file_type") == "bam"]
  sample.group <- Biosample.vec[one.sample$Biosample]
  bam.suffix.vec <- sapply(afiles, "[[", "href")
  cat(sprintf("%4d / %4d %s\n", sample.i, nrow(some.samples), sample.group))
  for(bam.suffix in bam.suffix.vec){
    sample.id <- sub(".bam$", "", basename(bam.suffix))
    alignments.bam <- file.path(
      "H3K4me3_TDH_ENCODE", "samples",
      sample.group, sample.id, "alignments.bam")
    if(!file.exists(alignments.bam)){
      dir.create(dirname(alignments.bam), showWarnings=FALSE, recursive=TRUE)
      u <- paste0("https://www.encodeproject.org", bam.suffix)
      download.file(u, alignments.bam, method="wget")
      request <- GET(u)
      stop_for_status(request)
      writeBin(content(request), alignments.bam)
    }
  }
}

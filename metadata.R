
download.file("http://epigenomesportal.ca/cgi-bin/api/getReleases.py", "releases.json")
library(RJSONIO)

releases.list <- fromJSON("releases.json")

makeRow <- function(L){
  L <- as.list(L)
  length.vec <- sapply(L, length)
  if(any(1 < length.vec)){
    print(L)
    stop("too big")
  }
  L[length.vec==0] <- NA
  do.call(data.table, L)
}
makeDT <- function(L){
  dt <- do.call(rbind, lapply(L, makeRow))
  if(is.character(names(L))){
    dt$name <- names(L)
  }
  dt
}
(releases.dt <- makeDT(releases.list))

for(release.i in 1:nrow(releases.dt)){
  release <- releases.dt[release.i]
  release.json <- file.path("releases", paste0(release$id, ".json"))
  if(!file.exists(release.json)){
    u <- paste0("http://epigenomesportal.ca/cgi-bin/api/getDataHub.py?data_release_id=", release$id)
    dir.create(dirname(release.json), showWarnings=FALSE)
    download.file(u, release.json)
  }
  types.csv <- file.path("releases", paste0(release$id, "-types.csv"))
  if(!file.exists(types.csv)){
    cat(sprintf("%4d / %4d id=%d\n", release.i, nrow(releases.dt), release$id))
    release.list <- fromJSON(release.json)
    ##samples.dt <- makeDT(release.list$samples, c("cell_type", "tissue_type"))
    type.dt <- do.call(rbind, lapply(release.list$datasets, function(L){
      out.list <- as.list(L$ihec_data_portal[c("cell_type", "cell_type_category")])
      out.list$url <- sapply(L$browser, function(l)l[[1]]$big_data_url)
      do.call(data.table, out.list)
    }))
    fwrite(type.dt, types.csv)
  }
}

  

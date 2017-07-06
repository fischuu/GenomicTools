getRSLocation <- function(rs, species="Gallus_gallus",url="jul2016.archive.ensembl.org"){
#  rm(list=ls())
#rs="rs14225498"
#species="Gallus_gallus"
#url="jul2016.archive.ensembl.org"
  
  
  ensemblReturn <- readLines(paste("http://",url,"/",species,"/Variation/Explore?db=core;v=",rs,sep=""))
  tmp <- ensemblReturn[grepl('>Location</div><div class=\"rhs\"><p>',ensemblReturn)]
  out <- strsplit(strsplit(tmp,'>Location</div><div class=\"rhs\"><p>')[[1]][2],"<span class=\"text_separator\"")[[1]][1]

  out <- strsplit(strsplit(out,"<b>")[[1]][[2]],"</b>")[[1]][1]
  out <- strsplit(out,":")[[1]]
  out <- list(Chromosome=out[1], BP=out[2])
  out
}
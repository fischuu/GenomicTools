# Functions based on the REST API from KEGG, description can be found here:
# http://www.kegg.jp/kegg/rest/keggapi.html

getKEGGOrganisms <- function(url="http://rest.kegg.jp/list/organism"){
  tmp <- readLines(url)
  tmp <- strsplit(tmp,"\t")
  out <- data.frame(KEGG=sapply(tmp,"[",1),
                    code=sapply(tmp,"[",2),
                    name=sapply(tmp,"[",3),
                    ontology=sapply(tmp,"[",4))
  out
}

getKEGGOrthologOverview <- function () 
{
    tmp <- readLines("http://rest.kegg.jp/list/KO/")
    tmp <- strsplit(tmp, "\t")
    out <- data.frame(ortholog = sapply(tmp, "[", 1), description = sapply(tmp, 
        "[", 2))
    out
}

getKEGGPathwayOverview <- function(code="hsa"){
  tmp <- readLines(paste("http://rest.kegg.jp/list/pathway/",code,sep=""))
  tmp <- strsplit(tmp,"\t")
  out <- data.frame(pathway=sapply(tmp,"[",1),
                    description=sapply(tmp,"[",2))
  out
}

getKEGGInformation <- function(x){
  tmp <- readLines(paste("http://rest.kegg.jp/get/",x,sep=""))
  output <- list()
  for(i in 1:length(tmp)){
  # This is the case that we find a new headline
    if(substr(tmp[i],1,1)!=" "){
    # We found a new FIELD entry, store the old one into the output list
      if(i>1){
        output[[eval(tmpName)]] <- tmpEntry
      }
    # Extract the inforamtion, first field the the headline, third field the entry and the rest the explanation.
        tmpLine <- unique(strsplit(tmp[i], " ")[[1]])
    # Create a temporary matrix that stores the information for that particular field
        tmpEntry <- data.frame(Entry = tmpLine[3],
                               Description = paste(tmpLine[4:length(tmpLine)], collapse=" "), stringsAsFactors = FALSE)
        tmpName <- tmpLine[1]
  # This is then the entry
    } else {
      # Extract the information, first field the the headline, third field the entry and the rest the explanation.
      tmpLine <- unique(strsplit(tmp[i], " ")[[1]])
      # Create a temporary matrix that stores the information for that particular field
      tmpEntry2 <- data.frame(Entry = tmpLine[2],
                              Description = paste(tmpLine[3:length(tmpLine)], collapse=" "), stringsAsFactors = FALSE)
      tmpEntry <- rbind(tmpEntry, tmpEntry2)
    }
  }
  output
}

getKEGGModule <- function(module){
  tmp <- readLines(paste("http://rest.kegg.jp/get/",module,sep=""))
  ENTRY <- tmp[which(grepl("ENTRY",tmp)==TRUE)]
  NAME <- tmp[which(grepl("NAME",tmp)==TRUE)]
  DEFINITION <- tmp[which(grepl("DEFINITION",tmp)==TRUE)]
  CLASS <- tmp[which(grepl("CLASS",tmp)==TRUE)]

  ORTHOLOGY <- which(grepl("ORTHOLOGY",tmp)==TRUE)
  leadingChr <- substr(tmp,0,1)
  # Now get the orthology rows
  if(length(ORTHOLOGY)>0){
    orthologyStart <- ORTHOLOGY
    checkThis <- ORTHOLOGY + 1
    while(leadingChr[checkThis]==" "){
      checkThis <- checkThis + 1
    }
    orthologyEnd <- checkThis - 1
    ORTHOLOGY.tmp <- tmp[orthologyStart:orthologyEnd]
    
    ORTHOLOGY.tmp[1] <- gsub("^ORTHOLOGY","         ",ORTHOLOGY.tmp[1])
    ORTHOLOGY.tmp <- trimws(ORTHOLOGY.tmp)
    
    ORTHOLOGY.tmp <- strsplit(ORTHOLOGY.tmp, " ")
    
    ORTHOLOGY.out <- data.frame(Orthology=sapply(ORTHOLOGY.tmp,"[",1),
                                Desc=trimws(sapply(sapply(ORTHOLOGY.tmp,"[",-1), paste, collapse=" "))
                                )
    
  } else {
    ORTHOLOGY.out <- "Not available"
  }
  out <- list(Entry = ENTRY,
              Name = NAME,
              Definition = DEFINITION,
              Class = CLASS,
              Orthology = ORTHOLOGY.out)
  out
}
  
getKEGGPathway <- function(pathway){
  tmp <- readLines(paste("http://rest.kegg.jp/get/",pathway,sep=""))
  ENTRY <- tmp[which(grepl("^ENTRY",tmp)==TRUE)]
  NAME <- tmp[which(grepl("^NAME",tmp)==TRUE)]
  DESCRIPTION <- tmp[which(grepl("^DESCRIPTION",tmp)==TRUE)]
  CLASS <- tmp[which(grepl("^CLASS",tmp)==TRUE)]
  PATHWAY_MAP <- tmp[which(grepl("^PATHWAY_MAP",tmp)==TRUE)]
  ORGANISM <- tmp[which(grepl("^ORGANISM",tmp)==TRUE)]
  GENE <- which(grepl("^GENE",tmp)==TRUE)
  COMPOUND <- which(grepl("^COMPOUND",tmp)==TRUE)
  MODULE <- which(grepl("^MODULE",tmp)==TRUE)
  ORTHOLOGY <- which(grepl("^ORTHOLOGY",tmp)==TRUE)
  leadingChr <- substr(tmp,0,1)
  # Now get the gene rows
   if(length(GENE)>0){
      geneStart <- GENE
      checkThis <- GENE + 1
      while(leadingChr[checkThis]==" "){
        checkThis <- checkThis + 1
      }
      geneEnd <- checkThis - 1
      GENE <- tmp[geneStart:geneEnd]
   # Now bring the data into some better form
      GENE[1] <- gsub("GENE", "    ",GENE[1])
      GENE <- strsplit(GENE,";")
      GENEtmp1 <- strsplit(trim.leading(sapply(GENE,"[",1)),"  ")
      GENEtmp2 <- trim.leading(sapply(GENE,"[",2))
      GENE <- data.frame(GeneID = sapply(GENEtmp1,"[",1),
                         GeneName = sapply(GENEtmp1,"[",2),
                         GeneDesc = GENEtmp2)
   } else {
      GENE <- data.frame(GeneID = "Not available",
                         GeneName = "Not available",
                         GeneDesc = "Not available")
   }     
  # Now get the compound rows
    if(length(COMPOUND)>0){
      compoundStart <- COMPOUND
      checkThis <- COMPOUND + 1
      while(leadingChr[checkThis]==" "){
        checkThis <- checkThis + 1
      }
      compoundEnd <- checkThis - 1
      COMPOUND <- tmp[compoundStart:compoundEnd]
    } else {
      COMPOUND <- "Not available"
    }
  # Now get the module rows
  if(length(MODULE)>0){
    moduleStart <- MODULE
    checkThis <- MODULE + 1
    while(leadingChr[checkThis]==" "){
      checkThis <- checkThis + 1
    }
    moduleEnd <- checkThis - 1
    MODULE.tmp <- tmp[moduleStart:moduleEnd]
    
    MODULE.tmp[1] <- gsub("MODULE","      ",MODULE.tmp[1])
    MODULE.tmp <- trimws(MODULE.tmp)
    
    MODULE.out <- data.frame(Module=substr(MODULE.tmp,1,6),
                            Desc=trimws(substr(MODULE.tmp, 7,nchar(MODULE.tmp))))
    
  } else {
    MODULE.out <- "Not available"
  }
  if(length(ORTHOLOGY)>0){
    orthologyStart <- ORTHOLOGY
    checkThis <- ORTHOLOGY + 1
    while(leadingChr[checkThis]==" "){
      checkThis <- checkThis + 1
    }
    orthologyEnd <- checkThis - 1
    ORTHOLOGY.tmp <- tmp[orthologyStart:orthologyEnd]
    
    ORTHOLOGY.tmp[1] <- gsub("^ORTHOLOGY","         ",ORTHOLOGY.tmp[1])
    ORTHOLOGY.tmp <- trimws(ORTHOLOGY.tmp)
    
    ORTHOLOGY.tmp <- strsplit(ORTHOLOGY.tmp, " ")
    
    ORTHOLOGY.out <- data.frame(Orthology=sapply(ORTHOLOGY.tmp,"[",1),
                                Desc=trimws(sapply(sapply(ORTHOLOGY.tmp,"[",-1), paste, collapse=" "))
    )
    
  } else {
    ORTHOLOGY.out <- "Not available"
  }
  out <- list(Entry = ENTRY,
              Name = NAME,
              Description = DESCRIPTION,
              Class = CLASS,
              Pathway_map= PATHWAY_MAP,
              Organism = ORGANISM,
              Module = MODULE.out,
              Orthology = ORTHOLOGY.out,
              Gene = GENE,
              Compound = COMPOUND)
  out
}

getKEGGPathwayImage <- function(pathway, folder=NULL){
  pathway <- gsub("path:","", pathway)
  if(is.null(folder)) folder <- getwd()
  filename <- paste(pathway,".png",sep="")
  dlURL <- paste("http://rest.kegg.jp/get/",pathway,"/image",sep="")
  download.file(url=dlURL, destfile= file.path(folder,filename), mode="wb")
}

getKEGGModuleOverview <- function () 
{
  tmp <- readLines("http://rest.kegg.jp/list/module/")
  tmp <- strsplit(tmp, "\t")
  out <- data.frame(module = gsub("md:", "", sapply(tmp, "[", 1)), description = sapply(tmp, "[", 2))
  out
}

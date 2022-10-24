library("pathview")

# DEobject: The output of edgeR

# For now inputVector takes a vector with fold changes of the genes. The gene names are expected to be in the names() of the vector

inputVector <- runif(0,2, n=10)
names(inputVector) <- c("CS", "ACLY", "ACO2", "ACO1", "IDH2", "IDH1", "IDH3A", "IDH3B", "IDH3G", "OGDHL")

plotKEGGPathway <- function(pathway="bta00190",
                            DEobject=NA,
                            FDR=0.05,
                            outfolder=NA,
                            verbose=TRUE){

  if(length(DEobject)==1) if(is.na(DEobject)) stop("Please provide a DE object from edgeR for visualisation")
  
  genesOI <- which(p.adjust(DEobject$table$PValue, method="BH") < FDR)
  
  genesToMark <- DEobject$table$logFC[genesOI]
  names(genesToMark) <- DEobject$genes$GeneName[genesOI]
  
# Set the correct output dir
  if(!is.na(outfolder)){
    prevwd <- getwd()
    setwd(outfolder)
  }

  
#Get the pathway details
  pathway.details <- getKEGGPathway(pathway)
  if(verbose) cat("Pathway name:", pathway.details$Name ,"\n")
  pathway.genes <- pathway.details$Gene

  # Now filter out the pathways we want to mark
  if(verbose) cat("Genes found in pathway  : ", nrow(pathway.genes),"\n")
  if(verbose) cat("Significant genes found : ", length(genesToMark),"\n")

  genesToMark <- genesToMark[is.element(names(genesToMark), pathway.genes$GeneName)]
  if(verbose) cat("Genes found to mark     : ", length(genesToMark),"\n")
  
  pv.out <- pathview(gene.data = log(genesToMark),
                     pathway.id = gsub('[a-zA-Z]+', '', pathway),
                     species = gsub('[0-9]+', '', pathway),
                     gene.idtype="SYMBOL",
                     kegg.native = TRUE)

  if(!is.na(outfolder)) setwd(prevwd)
  pv.out
}
  
plotKEGGPathway <- function(pw="ko00020",
                            orthologsToMark=keggortholog.SE.interesting,
                            annotation=annotOI.KO.SE,
                            KEGG.var="V9",
                            Gene.var="GeneID",
                            expression=prodigal_counts_filtered.SE.tmm,
                            groups=emmSE,
                            outdir="/work/users/ejo138/RumenPredict/RuminomicsHighLow/KEGG.SE"){
  # First get the pathway Orthologs
  #tmp <- getKEGGPathway(pw)$Gene
  #pwOrthologs <- gsub("]","",paste0("ko:",sapply(strsplit(sapply(strsplit(tmp$GeneDesc, "\\[KO:"),"[",2),"\\] "),"[",1)))
  tmpPW <- getKEGGPathway(pw)
  tmpOrthologs <- tmpPW$Orthology
  
  if(as.matrix(tmpOrthologs)[1,1]=="Not available"){
    tmp <- tmpPW$Module
    KEGGOrthologs <- c()
    if(as.matrix(tmp)[1,1]=="Not available"){
      KEGGOrthologs <- "Not available"
    } else {
      for(j in 1:nrow(tmp)){
        moduleOrthologs <- getKEGGModule(tmp[j,1])$Orthology[,1]
        KEGGOrthologs <- unique(c(KEGGOrthologs,unique(unlist(strsplit(moduleOrthologs,",")))))
        KEGGOrthologs <- unique(unlist(strsplit(KEGGOrthologs,"\\+")))
      }
    }
  } else {
    KEGGOrthologs <- unique(tmpPW$Orthology[,1])
  }
  
  tmpPW$KEGGorthologs <- KEGGOrthologs
  
  # Now filter out the pathways we want to mark
  orthologsToMark <- orthologsToMark[is.element(orthologsToMark, tmpPW$KEGGorthologs)]
  foldChanges <- rep(0, length(orthologsToMark))
  
  # Now get the expression data for the figure
  for(i in 1:length(orthologsToMark)){
    tmp <- get(KEGG.var, annotation)      
    orthologGenes <- get(Gene.var,annotation)[grep(orthologsToMark[i],tmp)]
    #        orthologGenes.expressions <- expression[is.element(rownames(expression), orthologGenes),]
    geneExpressions <- matrix(0, ncol=2, nrow=length(orthologGenes))
    for(j in 1:length(orthologGenes)){
      tmpExpr <- expression[rownames(expression)==orthologGenes[j],]
      if(sum(names(tmpExpr)==names(groups))!=length(tmpExpr))stop("Error: Mismatch with the columns between groups and expressions")
      group1 <- tmpExpr[is.element(names(tmpExpr), names(groups)[groups==1])]
      group2 <- tmpExpr[is.element(names(tmpExpr), names(groups)[groups==3])]
      geneExpressions[j,1] <- mean(group1)
      geneExpressions[j,2] <- mean(group2)
    }
    tmpExp <- apply(geneExpressions,2,sum, na.rm=TRUE)
    foldChanges[i] <- tmpExp[1]/ tmpExp[2]
  }
  names(foldChanges) <- gsub("ko:","",orthologsToMark)
  cat(foldChanges, "\n")
  pv.out <- pathview(gene.data = log(foldChanges),
                     pathway.id = gsub("map","",pw),
                     species = "ko",
                     kegg.native = TRUE)
}

plotThesePWs <- gsub("path:","",ISpathways)
for(i in 1:length(plotThesePWs)){
  setwd("/work/users/ejo138/RumenPredict/RuminomicsHighLow/KEGG.SE")
  plotKEGGPathway(pw=plotThesePWs[i],
                  orthologsToMark=keggortholog.SE.interesting,
                  annotation=annotOI.KO.SE,
                  KEGG.var="V9",
                  Gene.var="GeneID",
                  expression=prodigal_counts_filtered.SE.tmm,
                  groups=emmSE)
  
  setwd("/work/users/ejo138/RumenPredict/RuminomicsHighLow/KEGG.IT")
  plotKEGGPathway(pw=plotThesePWs[i],
                  orthologsToMark=keggortholog.IT.interesting,
                  annotation=annotOI.KO.IT,
                  KEGG.var="V9",
                  Gene.var="GeneID",
                  expression=prodigal_counts_filtered.IT.tmm,
                  groups=emmIT)
  
  setwd("/work/users/ejo138/RumenPredict/RuminomicsHighLow/KEGG.UK")
  plotKEGGPathway(pw=plotThesePWs[i],
                  orthologsToMark=keggortholog.UK.interesting,
                  annotation=annotOI.KO.UK,
                  KEGG.var="V9",
                  Gene.var="GeneID",
                  expression=prodigal_counts_filtered.UK.tmm,
                  groups=emmUK)
}
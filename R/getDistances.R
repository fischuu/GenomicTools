# library("GenomicTools")
# 
# data <- importVCF("/home/fischuu/ownCloud/Luke/Projects/Genotype/Whitefish/GSC.vcf")
# 
# #test <- readLines("/home/fischuu/ownCloud/Luke/Projects/Genotype/Raccoon/plink.mendel")
# 
# importPlinkMendel <- function(file){
#   
#   x <- readLines(file)
#   
#   header <- strsplit(x[1], "\t")[[1]]
# 
#   splitPlink <- function(x){
#     tmp <- strsplit(x, "\t")[[1]]
#     out <- c(tmp[1:6], paste(tmp[7:11],collapse=""), tmp[12:13])
#     out
#   }
#     
#   
#   out <- do.call(rbind,lapply(x[-1], splitPlink))
#   colnames(out) <- header
#   
#   out
# }
# 
# plinkMendel  <- importPlinkMendel(file="/home/fischuu/ownCloud/Luke/Projects/Genotype/Whitefish/plink.mendel")
# plinkMendel[,5] <- gsub("MOCKREFGENOME", "MockRefGenome", plinkMendel[,5])
# #plinkMendel[,5] <- gsub("-SNV", "", plinkMendel[,5])
# plinkMendel[,5] <- gsub(":", ".", plinkMendel[,5])
# 
# errorSNPs <- data$map[is.element(data$map$snp.names,plinkMendel[,5]),]
# 
# errorSNPs.info <- data$genotypesInfo[is.element(data$map$snp.names,plinkMendel[,5]),]
# 
# coverage <- c()
# for(i in 1:nrow(plinkMendel)){
#   nrow <- which(errorSNPs$snp.names==plinkMendel[i,5])
#   #ncol <- which(colnames(errorSNPs.info)==plinkMendel[i,3])
#   coverage[i] <- errorSNPs.info[nrow,get(plinkMendel[i,3])]
# }
# 
# errorCoverage <- apply(sapply(strsplit(coverage, ","), as.numeric),2,sum)
# 
# totalCoverage <- c()
# for(i in 1:nrow(data$genotypesInfo)){
# 
#   tmp <- data$genotypesInfo[i,]
#   totalCoverage[i] <- mean(apply(sapply(strsplit(as.vector(as.matrix(tmp)),","),as.numeric),2,sum),na.rm=TRUE)
# }
# 
# nper <- 10
# randomCoverage <- matrix(c, ncol=nper, nrow=length(coverage))
# for(i in 1:1){}
# 
# 
# 
# # vcf: VCF object
# # n: n-fold comparison (2 or 3)
# # In case n=2 we just calculate the pairwise distances (abs sum over the snp-wise differences)
# # In case n=3 we add then still another layer over it and assume sample i and j are mother/father, and take the average of those and then we calculate the distance of this
# # average to the third sample
# getDistances <- function(vcf, n=2, n.snps=100){
#   
#   numericMatrix <- apply(as.matrix(data$genotypes[,1:n.snps]),1,as.numeric)
#   numericMatrix[numericMatrix==3] <- NA
#   
#   if(n==2){
#     out <- matrix(0,nrow(vcf$genotypes), ncol=nrow(vcf$genotypes))
#     for(i in 1:(nrow(vcf$genotypes)-1)){
#         diffMatrix <- abs(t(t(numericMatrix) - numericMatrix[i,])
#         out[i,] <- apply(diffMatrix,1,sum, na.rm=TRUE)
#                         
#             for(j in (i+1):nrow(vcf$genotypes)){
#         diffVector <- apply(numericMatrix,1,diff, na.rm=TRUE)
#         out[i,j] <- sum(abs(diffVector), na.rm=TRUE)/(length(diffVector)-sum(is.na(diffVector)))
#         out[j,i] <- out[i,j]
#         diag(out) <- 0
#       }
#     cat("Sample",i,"finished\n")
#     }
#   } else if(n==3){
#     out <- array(0,c(rep(nrow(vcf$genotypes),3)))
#     for(i in 1:(nrow(vcf$genotypes)-1)){
#       for(j in (i+1):nrow(vcf$genotypes)){
#         numericMatrix <- apply(as.matrix(data$genotypes[c(i,j),1:n.snps]),1,as.numeric)
#         numericMatrix[numericMatrix==3] <- NA
#      # This is the average genotype vector of the assumed partens
#         meanVector <- apply(numericMatrix,1,mean, na.rm=TRUE)
#         for(k in 1:nrow(vcf$genotypes)){
#         # Here we combine the offspring genotype vecotr with the assumed mean vector of the parents
#           numericMatrix <- cbind(meanVector, as.numeric(vcf$genotypes[k,1:n.snps]))
#           numericMatrix[numericMatrix==3] <- NA
#           diffVector <- apply(numericMatrix,1,diff, na.rm=TRUE)
#           out[i,j,k] <- sum(abs(diffVector), na.rm=TRUE)/(length(diffVector)-sum(is.na(diffVector)))
#           out[j,i,k] <- out[i,j,k]
#         }
#       }
#       cat("Sample",i,"finished\n")
#     }
#   }
#   out 
# }
# 
# res <- getDistances(data, n=3)
# 
# estimatePedigree <- function(x, vcf){
#   out <- matrix("NA", nrow=dim(x)[3], ncol=3)
#   for(i in 1:dim(x)[3]){
#     tmp <- x[,,i]
#     tmp[i,] <- 1
#     tmp[,i] <- 1
#     diag(tmp) <- 1
#     parentsIndex <- arrayInd(which.min(tmp), dim(tmp))
#     out[i,1] <- rownames(vcf$genotypes)[i]
#     out[i,2] <- rownames(vcf$genotypes)[parentsIndex[1,1]]
#     out[i,3] <- rownames(vcf$genotypes)[parentsIndex[1,2]]
#   }
#   out
# }
# 
# estimatePedigree(res, data)
# 
# 
# # Do some indeep stats for the problematic SNPs
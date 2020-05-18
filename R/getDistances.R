library("GenomicTools")

data <- importVCF("/home/fischuu/ownCloud/Luke/Projects/Genotype/Raccoon/GSC.vcf")

# vcf: VCF object
# n: n-fold comparison (2 or 3)
# In case n=2 we just calculate the pairwise distances (abs sum over the snp-wise differences)
# In case n=3 we add then still another layer over it and assume sample i and j are mother/father, and take the average of those and then we calculate the distance of this
# average to the third sample
getDistances <- function(vcf, n=2, n.snps=100){
  
  if(n==2){
    out <- matrix(0,nrow(vcf$genotypes), ncol=nrow(vcf$genotypes))
    for(i in 1:(nrow(vcf$genotypes)-1)){
      for(j in (i+1):nrow(vcf$genotypes)){
        numericMatrix <- apply(as.matrix(data$genotypes[c(i,j),1:n.snps]),1,as.numeric)
        numericMatrix[numericMatrix==3] <- NA
        diffVector <- apply(numericMatrix,1,diff, na.rm=TRUE)
        out[i,j] <- sum(abs(diffVector), na.rm=TRUE)/(length(diffVector)-sum(is.na(diffVector)))
        out[j,i] <- out[i,j]
        diag(out) <- 0
      }
    cat("Sample",i,"finished\n")
    }
  } else if(n==3){
    out <- array(0,c(rep(nrow(vcf$genotypes),3)))
    for(i in 1:(nrow(vcf$genotypes)-1)){
      for(j in (i+1):nrow(vcf$genotypes)){
        numericMatrix <- apply(as.matrix(data$genotypes[c(i,j),1:n.snps]),1,as.numeric)
        numericMatrix[numericMatrix==3] <- NA
     # This is the average genotype vector of the assumed partens
        meanVector <- apply(numericMatrix,1,mean, na.rm=TRUE)
        for(k in 1:nrow(vcf$genotypes)){
        # Here we combine the offspring genotype vecotr with the assumed mean vector of the parents
          numericMatrix <- cbind(meanVector, as.numeric(vcf$genotypes[k,1:n.snps]))
          numericMatrix[numericMatrix==3] <- NA
          diffVector <- apply(numericMatrix,1,diff, na.rm=TRUE)
          out[i,j,k] <- sum(abs(diffVector), na.rm=TRUE)/(length(diffVector)-sum(is.na(diffVector)))
          out[j,i,k] <- out[i,j,k]
        }
      }
      cat("Sample",i,"finished\n")
    }
  }
  out 
}

res <- getDistances(data, n=3)

estimatePedigree <- function(x, vcf){
  out <- matrix("NA", nrow=dim(x)[3], ncol=3)
  for(i in 1:dim(x)[3]){
    tmp <- x[,,i]
    tmp[i,] <- 1
    tmp[,i] <- 1
    diag(tmp) <- 1
    parentsIndex <- arrayInd(which.min(tmp), dim(tmp))
    out[i,1] <- rownames(vcf$genotypes)[i]
    out[i,2] <- rownames(vcf$genotypes)[parentsIndex[1,1]]
    out[i,3] <- rownames(vcf$genotypes)[parentsIndex[1,2]]
  }
  out
}

estimatePedigree(res, data)


# Do some indeep stats for the problematic SNPs
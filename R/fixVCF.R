# This function takes a vcf object and fills the missing variants, based on the read coverage ratios. By default it uses 1/3 and 2/3 as threshold for discretisation 0/0, 0/1 and 1/1

fixVCF <- function(vcf){
  
  for(sampleRun in 1:nrow(vcf$genotypes)){
    missingLoci <- which(vcf$genotypes[sampleRun,]=="03")
    
    counts <- as.vector(as.matrix(vcf$genotypesInfo[,..sampleRun]))
    
    counts_split <- strsplit(counts, ",")
    
    counts_ref <- sapply(counts_split, "[", 1)
    counts_alt <- sapply(counts_split, "[", 2)
    
    # Explicitly replace "." with NA before conversion
    counts_ref[counts_ref == "."] <- NA
    counts_alt[counts_alt == "."] <- NA
    
    # Convert to numeric without warnings
    counts_ref <- as.numeric(counts_ref)
    counts_alt <- as.numeric(counts_alt)
   
    # Calculate the ratio
    counts_ratios <- ifelse((counts_ref + counts_alt) == 0, NA, counts_ref / (counts_ref + counts_alt))
    
    # Replace missing genotype "03" based on counts_ratios thresholds
    new_values <- as.vector(as.matrix(vcf$genotypes[sampleRun, ]))
    
    new_values[missingLoci] <- ifelse(
      is.na(counts_ratios[missingLoci]), "03", # Keep "03" if NA
      ifelse(counts_ratios[missingLoci] >= 0.66, "00",
             ifelse(counts_ratios[missingLoci] >= 0.33, "01", "02"))
    )
    
    # Assign back row-wise using data.table syntax
    vcf$genotypes[sampleRun, (colnames(vcf$genotypes)) := as.list(new_values)]
  }
  
  return(vcf) # Return the modified VCF object
  
}


#This script loads the output of the SNP quality profile plugin from TASSEL in order to calculate
#a quality score based on minimum MAF (minor allele frequency), minimum F (inbreeding coefficient),
#and excluding InDels
#Change values in minF, no.gap, and minMAF below to use different thresholds
#The input file is called "outputStats.txt"
#Written by Marco Pessoa-Filho
#

#load file and exclude NaN's if they are present
snp.qualME <- read.delim("outputStatsME.txt", header = TRUE)
filtered_snps <- (snp.qualME[snp.qualME$inbredF_DGE2 != "NaN",])

#Define thresholds
minF <- filtered_snps$inbredF_DGE2 > 0.8
no.gap <- filtered_snps$gapDepthProp == 0 
minMAF <- filtered_snps$minorAlleleFreqGE2 > 0.2

#SNPs that follow all three rules above will have a score of 3
QUALITYSCORE_ME <- c()
rownumber <- nrow(filtered_snps)
for(j in 1:rownumber){
  QUALITYSCORE_ME[j] <- minF[j] + no.gap[j]  + minMAF[j]
}

#add a column with quality scores to table
filtered_snps$QUALITYSCORE <- QUALITYSCORE_ME

#generate quality score table to be used as input in the production pipeline
myQsFile_ME <- c()
myQsFile_ME$CHROM <- filtered_snps[,1]
myQsFile_ME$POS <- filtered_snps[,2]
myQsFile_ME$QUALITYSCORE <- filtered_snps[,14]
myQsFile_ME <- as.data.frame(myQsFile_ME, stringsAsFactors = TRUE)


#the commands below fix chromosome names from 1 to 9 to match those used in the reference
myQsFile_ME$CHROM[myQsFile_ME$CHROM == 1] <- as.character("01")
myQsFile_ME$CHROM[myQsFile_ME$CHROM == 2] <- as.character("02")
myQsFile_ME$CHROM[myQsFile_ME$CHROM == 3] <- as.character("03")
myQsFile_ME$CHROM[myQsFile_ME$CHROM == 4] <- as.character("04")
myQsFile_ME$CHROM[myQsFile_ME$CHROM == 5] <- as.character("05")
myQsFile_ME$CHROM[myQsFile_ME$CHROM == 6] <- as.character("06")
myQsFile_ME$CHROM[myQsFile_ME$CHROM == 7] <- as.character("07")
myQsFile_ME$CHROM[myQsFile_ME$CHROM == 8] <- as.character("08")
myQsFile_ME$CHROM[myQsFile_ME$CHROM == 9] <- as.character("09")

#the file below can be used as input in the TASSEL GBSv2 production pipeline
write.table(myQsFile_ME, file="myQsFileME", sep="\t", quote = FALSE, row.names = FALSE)

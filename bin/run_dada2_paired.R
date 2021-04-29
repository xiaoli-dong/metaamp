#!/usr/bin/env Rscript
library("optparse")
library(dada2)
packageVersion("dada2")
library(ggplot2)

option_list <- list(
  make_option(c("-d", "--basedir"), type = "character", default = ".", help = "output directory"),
  make_option(c("-r", "--seqdir"), type = "character", default = ".", help = "directory name for the run"),
  make_option(c("-n", "--aname"), type = "character", default = "", help = "analysis name"),
  make_option(c("-p", "--pattern"), type = "character", default = "", help = "file name pattern"),	
  make_option(c("-i", "--id"), type = "character", default = "", help = "run id, it must be unique for each run"),
  make_option("--truncLen_F", type="integer", default=0),
  make_option("--truncLen_R", type="integer", default=0),
  make_option("--maxN", type="integer", default=0),	
  make_option("--maxEE_F", type="integer", default=1),
  make_option("--maxEE_R", type="integer", default=1),	
  make_option("--truncqF", type="integer", default=2),
  make_option("--truncqR", type="integer", default=2)			                   
)
opt = parse_args(OptionParser(option_list = option_list, usage = "This Script is a test for arguments!"))

print (opt);

setwd(opt$basedir)

####################### Change this for each Run #######################

path<- opt$seqdir
seqtabname = paste0(opt$id, ".", opt$aname, ".rds");

#Pattern will changes based upon your file name

fnFs <- sort(list.files(path, pattern=paste0("_R1.", opt$pattern), full.names = TRUE))
fnRs <- sort(list.files(path, pattern=paste0("_R2.", opt$pattern), full.names = TRUE))
fnFs
fnRs

#What the sample names are infered from file name
basename(fnFs)
sample.names <- sapply(strsplit(basename(fnFs), paste0("_R1.", opt$pattern)), `[`, 1)
sample.names

#Filtered file name, you can modify the string between "" to change it
filtFs <- file.path(path, paste0(sample.names, "_R1.filt.", opt$id, ".fastq.gz"))
filtRs <- file.path(path, paste0(sample.names, "_R2.filt.", opt$id, ".fastq.gz"))
filtFs
filtRs

#Double checking all the file pairs are there
if(length(fnFs) != length(fnRs)) stop("Forward and reverse files do not match.")

#Filtering... Play around with the trunclen and trunQ based upon your run. 
#maxEE is subjective to run too. Normally the larger second number allows for more error in reverse read
#out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(280,260), maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, compress=TRUE, verbose=TRUE, multithread=TRUE)
#out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(opt$truncLen_F, opt$truncLen_R), maxN=opt$maxN, maxEE=c(opt$maxEE_F,opt$maxEE_R), truncQ=c(opt$truncqF,opt$truncqR), rm.phix=TRUE, compress=TRUE, verbose=TRUE, multithread=TRUE)

if(opt$truncLen_F == 0 |  opt$truncLen_R == 0){
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncQ=c(opt$truncqF,opt$truncqR), maxN=opt$maxN, maxEE=c(opt$maxEE_F,opt$maxEE_R),  rm.phix=TRUE, compress=TRUE, verbose=TRUE, multithread=TRUE)
} else {
  out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(opt$truncLen_F, opt$truncLen_R), maxN=opt$maxN, maxEE=c(opt$maxEE_F,opt$maxEE_R), truncQ=c(opt$truncqF,opt$truncqR), rm.phix=TRUE, compress=TRUE, verbose=TRUE, multithread=TRUE)
}

#Check out what happened, how many reads lost. Based upon this, you can go back and play around with the parameters in previos step
out
write.table(gsub("^(./+)", "", c(unname(filtFs),unname(filtRs))), file=paste0(opt$id, ".", "dada2_filter_trim_file.list.txt"), row.names = FALSE, na='', col.names = FALSE,  quote=F)
write.table(out, file=paste0(opt$id, ".", "dada2_filter_trim.stats.csv"), quote=F, sep=",", row.names=T, col.names=c("filename,readin","readout"))	

names(filtFs) <- sample.names
names(filtRs) <- sample.names
#Always set your own seed, helps with replication of data processing
set.seed(100)
#Following 2 steps will take a LOT of time
#Genereally, nread=1e6 is plenty, but a large number of samples you can use 2e6
errF <- learnErrors(filtFs, nbases=2e10, multithread=TRUE, randomize = TRUE, MAX_CONSIST = 20)
errR <- learnErrors(filtRs, nbases=2e10, multithread=TRUE, randomize = TRUE, MAX_CONSIST = 20)


mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names
#Check what the error profiles look like and how much do they diverge
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

# A for loop for large dataset, so that not the enitre seq data is loaded in memory
# It processes samples one by one
head(sample.names)                                    
for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  derepF <- derepFastq(filtFs[[sam]])
  ddF <- dada(derepF, err=errF, multithread=TRUE)
  ddF$map

  derepR <- derepFastq(filtRs[[sam]])
  ddR <- dada(derepR, err=errR, multithread=TRUE)

  merger <- mergePairs(ddF, derepF, ddR, derepR)
  #str(merger)
  
  head(merger)
  mergers[[sam]] <- merger
}

#Removing huge intermediate files
#rm(derepF); rm(derepR)

#Make seq table..akin to OTU table
cat("merge\n")
seqtab <- makeSequenceTable(mergers)
#How many seq variants you have
dim(seqtab)
cat("save rds\n");
saveRDS(seqtab, file=seqtabname)



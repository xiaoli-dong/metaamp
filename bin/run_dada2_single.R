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
  make_option("--truncLen", type="integer", default=150),
  make_option("--maxN", type="integer", default=0),	
  make_option("--maxEE", type="integer", default=1),
  make_option("--truncq", type="integer", default=2)			
  		                   
)
opt = parse_args(OptionParser(option_list = option_list, usage = "This Script is a test for arguments!"))

setwd(opt$basedir)

####################### Change this for each Run #######################
path<- opt$seqdir
seqtabname = paste0(opt$id, ".", opt$aname, ".rds");	   

#Pattern will changes based upon your file name
fnFs <- sort(list.files(path, pattern=opt$pattern, full.names = TRUE))
fnFs

#What the sample names are infered from file name
basename(fnFs)

sample.names <- sapply(strsplit(basename(fnFs), paste0(".", opt$pattern)), `[`, 1)
sample.names

#Filtered file name, you can modify the string between "" to change it
filtFs <- file.path(path, paste0(sample.names, ".filt.", opt$id, ".fastq.gz"))
filtFs

#Filtering... Play around with the trunclen and trunQ based upon your run. 
#maxEE is subjective to run too. Normally the larger second number allows for more error in reverse read
#out <- filterAndTrim(fnFs, filtFs, truncLen=200, maxN=0, maxEE=1, truncQ=2, rm.phix=TRUE, compress=TRUE, verbose=TRUE, multithread=TRUE)
out <- filterAndTrim(fnFs, filtFs, truncLen=opt$truncLen, maxN=opt$maxN, maxEE=opt$maxEE, truncQ=opt$truncq, rm.phix=TRUE, compress=TRUE, verbose=TRUE, multithread=TRUE)
#Check out what happened, how many reads lost. Based upon this, you can go back and play around with the parameters in previos step
out



#write.table(gsub("^(./+)", "", unname(filtFs)), file="dada2_filter_trim_file.list.txt", row.names = FALSE, na='', col.names = FALSE,  quote=F)
#write.table(out, file="dada2_filter_trim.stats.csv", quote=F, sep=",", row.names=T, col.names=c("filename,readin","readout"))

write.table(gsub("^(./+)", "", unname(filtFs)), file=paste0(opt$id, ".", "dada2_filter_trim_file.list.txt"), row.names = FALSE, na='', col.names = FALSE,  quote=F)
write.table(out, file=paste0(opt$id, ".", "dada2_filter_trim.stats.csv"), quote=F, sep=",", row.names=T, col.names=c("filename,readin","readout"))	
names(filtFs) <- sample.names

#Always set your own seed, helps with replication of data processing
set.seed(100)
#Following 2 steps will take a LOT of time
#Genereally, nread=1e6 is plenty, but a large number of samples you can use 2e6
errF <- learnErrors(filtFs, nbases=2e10, multithread=TRUE, randomize = TRUE, MAX_CONSIST = 20)

#Check what the error profiles look like and how much do they diverge
plotErrors(errF, nominalQ=TRUE)

dds <- vector("list", length(sample.names))
names(dds) <- sample.names

# A for loop for large dataset, so that not the enitre seq data is loaded in memory
# It processes samples one by one
head(sample.names)                                    
for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  derepF <- derepFastq(filtFs[[sam]])
  dds[[sam]] <- dada(derepF, err=errF, multithread=TRUE)
 }

#Removing huge intermediate files
#rm(derepF); rm(derepR)

#Make seq table..akin to OTU table
seqtab <- makeSequenceTable(dds)
#How many seq variants you have
dim(seqtab)
saveRDS(seqtab, file=seqtabname)

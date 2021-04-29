#!/usr/bin/env Rscript
library("optparse")
library(dada2)
packageVersion("dada2")

option_list <- list(
make_option(c("-d", "--basedir"), type = "character", default = ".", help = "output directory"),
make_option(c("-a", "--aname"), type = "character", default = "", help = "analysis name"),
make_option(c("-r", "--runids_str"), type = "character", default = "", help = "runids seperated by , for example: run1,run2 ")
)
opt = parse_args(OptionParser(option_list = option_list, usage = "This Script is a test for arguments!"))

basedir <- opt$basedir
setwd(basedir)
seqtabname = paste0(opt$aname, ".rds");
#string vector
rds_vector <- c();

runids <- unlist(strsplit(opt$runids_str, ","))
run_count <- length(runids);

if(run_count == 1){
rds_file_name <- paste0(basedir, "/", opt$runids_str, ".", opt$aname, ".rds");
cat(rds_file_name);
seqtab <- readRDS(rds_file_name);
class(seqtab)
} else{
for (runid in runids) {
rds_file_name <- paste0(basedir, "/", runid, ".", opt$aname, ".rds");
cat(rds_file_name);
rds_vector <- c(rds_vector, rds_file_name);
}
rds_vector
class(rds_vector);
seqtab <- mergeSequenceTables(tables=rds_vector)
class(seqtab)
}

#What's the size distribution of seq variants
table(nchar(getSequences(seqtab)))
#Chimera search
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

#What changed/how many chimeras removed
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)
#Seq variant lenght distribution after chimera checking
table(nchar(getSequences(seqtab.nochim)))
saveRDS(seqtab, file=seqtabname)

library(tibble)
seqTable <- seqtab.nochim
cs <- data.frame(t(colSums(seqTable)))
# Save the sequence variant table as a csv file
seqnum <- paste0("asv_", seq(ncol(seqTable))) #gives unique ID (asv1, asv2,...) to each sequence variant
cs <- data.frame(t(colSums(seqTable)))
fasta_header <- paste0(seqnum, "|", as.vector(cs[1,]))

uniqueSeqs <- as.list(colnames(seqTable)) #creates a list of the sequences
#seqTable.transposed <- t(seqTable) #transposes the matirx (asvs in rows, samples in columns)
#rownames(seqTable.transposed) <- as.character(seqnum) #changes rownames to asv IDs
colnames(seqTable) <- as.character(seqnum)
test <- rownames_to_column(as.data.frame(seqTable), var="Group")
test <- add_column(test, numasv=ncol(seqTable), .before=2)
test <- add_column(test, label="asv", .before=1)
write.table(test, file=paste0(opt$aname, ".shared"), sep="\t", row.names = FALSE, quote=FALSE)

# Write a fasta file containing the seq variants
library(seqinr)
write.fasta(uniqueSeqs, fasta_header, paste0(opt$aname,".asv.fasta"))

#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Long;
use FindBin;
my ($fasta, $analysisname, $scutoff, $cpus);

$scutoff = 0.97;

&GetOptions("fasta=s" =>\$fasta,
	    "s=f" =>\$scutoff,
            "an=s" =>\$analysisname,
	    "t=i" =>\$cpus
    );

($fasta && $analysisname) or 
    die "\nusage: $0 \n".
    "-fasta <fasta file>\n".
    "-an <analysis name>\n".
    "-s <similarity cutoff for generating OTUs, default is 0.97>\n".
    "-t <number of threads>\n".
    "cluster the input fasta file into otus\n";

my $bin = "$FindBin::RealBin";
my $usearch = "$bin/../programs/usearch64";
my $mothur = "$bin/../programs/mothur/mothur";

# dereplicate seqs
#my $cmd = "$usearch -fastx_uniques $fasta -threads 8 -fastaout $analysisname.derep.fasta -sizeout;";
my $cmd = "$usearch -fastx_uniques $fasta -threads $cpus -fastaout $analysisname.derep.fasta -sizeout;";
print STDERR "$cmd\n";
system $cmd;

# cluster OTUs and discarding singletons
$cmd = "$usearch -cluster_otus $analysisname.derep.fasta  -threads $cpus -otus $analysisname.otus.fasta -minsize 2 -relabel OTU_;";
print STDERR "$cmd\n";
system $cmd;

#$cmd = "$usearch -usearch_global $fasta -threads 8 -db $analysisname.otus.fasta -strand plus -id 0.97 -uc otu.map.uc -otutabout otu_table.txt -biomout otu_table.json";
$cmd = "$usearch -usearch_global $fasta -threads $cpus -db $analysisname.otus.fasta -strand plus -id $scutoff -uc otu.map.uc -otutabout otu_table.txt -biomout otu_table.json";

print STDERR "$cmd\n";
system $cmd;




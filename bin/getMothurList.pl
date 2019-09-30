#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;

my ($input, $fasta, $dissimilarity, $list, $group, $rep);


&GetOptions("i=s" =>\$input,
	    "fasta=s" =>\$fasta,
	    "ds=f" =>\$dissimilarity,
            "l=s" =>\$list,
	    "g=s" =>\$group,
	    "rep=s" =>\$rep
	    
    );

($input && $fasta && $list && $group && $rep) || 
    die "\nusage: $0\n".
    "-i <otu.map.uc USEARCH cluster format (UC)>\n".
    "-fasta <fasta format rep file output from metaamp>\n".
    "-ds <dissimilrity cutoff for generating OTUs, default value is: 0.03>\n".
    "-l <mothur list file name to produce>\n".
    "-g <mothur group file to produce>\n" .
    "-rep <reformatted fasta rep file, added otus size in the fasta header>\n".
    "Parse USEARCH cluster format (UC) file to Generate mothur compatile List and group file\n";

open UC, "<$input" or die "Could not open $input to read, $!\n";

my %otu2seqs = ();
my %seq2otuid = ();
my %seq2group = ();

#H       187     400     99.5    +       0       0       400M    M03263_22_000000000-ATTNV_1_1101_9855_1283;barcodelabel=LLMO_Batch_4                    OTU_188




while (<UC>) {
    if(/^H/){
	my @line = split(/\t+/, $_);
	#my ($otuid) = $line[9] =~ /OTU_(\d+)/;
	my ($otuid) = $_ =~ /OTU_(\d+)/;
	my $seqid = $line[8];
	
	#my($seqid, $groupid) = $line[8] =~ /^(\S+?);\w+?=(\w+)/;
	my ($groupid) = $seqid =~ /^\S+?;\w+?=(\w+)/;
	
	$seq2otuid{$seqid} = $otuid;
	push(@{$otu2seqs{$otuid}}, $seqid);
	
	$seq2group{$seqid} = $groupid;
    }
}
my $seq2groupsize = scalar keys %seq2group;

#print STDERR "$seq2groupsize\t$count\n";

close UC;
open LIST, ">$list" or die "Could not open $list to write, $!\n";
open GROUP, ">$group" or die "Could not open $group to write, $!\n";


foreach my $seq (sort {$seq2group{$a} cmp $seq2group{$b}} keys %seq2group){
    print GROUP "$seq\t$seq2group{$seq}\n";
}

print LIST "label\tnumOtus\t";

my $total_otu_count = keys  %otu2seqs;
#print the head of the list file
my %otuid2newid = ();
my $count = 1;
foreach my $otuid (sort {$a <=> $b} keys %otu2seqs){
    $otuid2newid{$otuid} = $count;
    print LIST "Otu", $count, "\t";
    $count++;
    #print LIST "Otu", $otuid, "\t";
    
}
print LIST "\n";

#print the otu one by one

print LIST "$dissimilarity\t", $total_otu_count, "\t";

foreach my $otuid (sort {$a <=> $b} keys %otu2seqs){
    
    print LIST join(",",@{$otu2seqs{$otuid}}), "\t";  
}
print LIST "\n";
close(LIST);
close(GROUP);

open FASTA, $fasta or die "Could not open $fasta file to read, $!\n";
open REP , ">$rep.fasta" or die "Cound not open $rep.fasta file to write, $!\n";
while (<FASTA>){
    next if !/\S+/;
    chomp;
    if(/^(>\S+?)(\d+)/){
	if(exists $otu2seqs{$2}){
	    my $size = scalar @{$otu2seqs{$2}};
	    my $newid = $otuid2newid{$2};
	    print REP ">OTU\_$newid|$size\n";
	}
    }
    else{

	print REP "$_\n";
    }
    
}

close(FASTA);
close(REP);

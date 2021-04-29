#!/usr/bin/env perl

use strict;
use warnings;

@ARGV == 4 or die "Usage: $0 <fastq file> <prefix of the fasta outfile> <quality threshold> <relable>\n",
                  "Writes a fasta file with lower case letters where quality is below threshold\n";
open(FASTQ, $ARGV[0]) or die "Cannot open $ARGV[0] for reading: $!\n";
open(FASTA, ">$ARGV[1].fasta") or die "Cannot open $ARGV[1] for writing: $!\n";
open(QUAL, ">$ARGV[1].qual") or die "Cannot open $ARGV[1].qual for writing: $!\n";
my $label = $ARGV[3];

while(<FASTQ>){
    chomp;
    my $line1 = $_;
    my $line2 = <FASTQ>;
    chomp($line2);
    
    my $line3 = <FASTQ>;
    chomp($line3);
    
    my $line4 = <FASTQ>;
    chomp($line4);
    $line1 =~ s/\s+//g;
    $line1 =~ s/^\@/>/;
    print FASTA "$line1;barcodelabel=$label\n";
    print QUAL "$line1;barcodelabel=$label\n";
    my @seq = ();
    if($line2 !~ /\S/){print STDERR $line1, "\n";}
    @seq = split //, $line2;
    my $masked_seq = "";
    my @qual_array = ();
    my $i = 0;
    
    for my $q (split //, $line4){
	my $phred_score = ord($q) - 33; # decode FASTQ quality encoding
	push(@qual_array, $phred_score);
	
	if($phred_score < $ARGV[2]){ # low qual
	    $masked_seq .= lc($seq[$i++]);
	}
	else{
	    $masked_seq .= uc($seq[$i++]);
	}
    }
    
    print FASTA $masked_seq, "\n";
    print QUAL join(" ", @qual_array), "\n"; 
}

close(FASTQ);
close(FASTA);
close(QUAL);
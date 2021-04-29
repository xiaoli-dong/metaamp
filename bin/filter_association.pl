#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Long;

my ($assoFile, $coef, $sig);
$coef = 0.6;
$sig = 0.01;

&GetOptions("a=s" =>\$assoFile,
	    "c=f" =>\$coef,
	    "p=f" =>\$sig
    );
($assoFile) or
    die "\nusage: $0 \n".
    "-a <otu.association output file>\n".
    "-c <correlation coefficient cutoff, the value should be between 0 and 1>\n".
    "-p <p value cutoff>\n";


open(ASSO, $assoFile) or die "Could not open $assoFile to read, $!\n";

while(<ASSO>){
    chomp;
    next if /Significance/;
    my @line = split(/\t/, $_);
    my $otuA = shift @line;
    my $otuB = shift @line;
    my $ab_coef = shift @line;
    my $ab_sig = shift @line;
    
    if(abs($ab_coef) >= $coef && $ab_sig <= $sig){ 
	print "$_\n";
    }
    
}
close(ASSO);


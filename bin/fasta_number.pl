#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Long;

my ($fasta, $prefix, $out, $tax);

&GetOptions("fasta=s" =>\$fasta,
	    "tax=s" =>\$tax,
            "prefix=s" =>\$prefix,
	    "out=s" =>\$out
    );

($fasta && $prefix && $out && $tax) or 
    die "\nusage: $0 \n".
    "-fasta <fasta file>\n".
    "-tax <classified tax file>\n".
    "-prefix <rename fasta seqname to prefixN>\n".
    "-out <output fasta and tax file prefix>\n".
    "rename the input fasta sequence name  into prefixN format, N is the input sequence order number\n";
open FASTA, "$fasta" or die "Could not open $fasta, $!\n";
open FOUT, ">$out.fasta" or die "Could not open $out.fasta, $!\n";
my $count = 0;
my %seqid2otuid = ();

while (<FASTA>){
    if(/^>(\S+)/){
	$count++;
	print FOUT ">$prefix$count\n";
	$seqid2otuid{$1} = "$prefix$count";
    }
    else{
	print FOUT $_;
    }
}
close(FASTA);
close(FOUT);

open TAX, "$tax" or die "Could not open $tax, $!\n";
open TOUT, ">$out.tax" or die "Could not open $out.tax, $!\n";
while (<TAX>){
    if(/^(\S+?)\s+(\S.*)$/){
	print TOUT $seqid2otuid{$1}, "\t$2\n";;
    }
   
}
close(TAX);
close(TOUT);

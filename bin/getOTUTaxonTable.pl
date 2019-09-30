#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;

my ($input1, $input2);


&GetOptions("i1=s" =>\$input1,
	    "i2=s" =>\$input2
    );

($input1 && $input2) || 
    die "\nusage: $0\n".
    "-i1 <taxonomic annotation of each otu: *.taxonomy>\n".
    "-i2 <mothur shared file>\n".
    "generate otu table summary file with taxonomic inforation\n";

#OTU_1|7275 Archaea(100);Euryarchaeota(100);Methanomicrobia(100);Methanomicrobiales(100);Methanoregulaceae(100);Methanoregula(100);
open TAXON, "<$input1" or die "Could not open $input1 to read, $!\n";
my %otus = ();

while (<TAXON>) {
    chomp;
    my ($otuid, $count, $taxon) = $_ =~ /^\w+?_?(\d+?)\|(\d+?)\s+?(\w.*)$/;
    $otus{$otuid}->{count} = $count;
    $otus{$otuid}->{taxon} = $taxon;
}
close(TAXON);

#OTUID 27 28 29
#OTU_1 2834 2553 1888

open TABLE, "<$input2" or die "Could not open $input2 to read, $!\n";
my %table = ();
my @samples = ();
my %totals = ();
while (<TABLE>) {
    next if /^label/;
    chomp;
    my @line = split(/\s+/, $_);
    my $distance = shift @line;
    my $samplename = shift @line;
    my $numOtus = shift @line;
   
    my $len = @line;
    
    for (my $i = 0; $i < $len; $i++) {
	my $num = $line[$i];
	my $id = $i+1;
	$table{$id}->{$samplename} = $num;
	if(exists $totals{$samplename}){
	    $totals{$samplename} += $num;
	}
	else{
	    $totals{$samplename} = $num;
	}
    }
    push(@samples, $samplename);
}
close(TABLE);
print "OTUID\tTotal\t". join("\t", @samples). "\t";
foreach (@samples){
    print "$_(%)\t";
}
print "Taxonomy\n";

#foreach my $otuid (sort {$a <=> $b} keys %table){
foreach my $otuid (sort {$a <=> $b} keys %otus){
    print "OTU_$otuid\t" ;
    print $otus{$otuid}->{count} if(exists $otus{$otuid}->{count});
#total of each otu
    print "\t";
    
    foreach  my $sample (@samples){
	if(exists $table{$otuid}->{$sample}){
	    print $table{$otuid}->{$sample}, "\t";
	}
	else{
	    print "0\t";
	}
	
    }
    foreach  my $sample (@samples){
        #print sprintf("%.4f", 100 * $table{$otuid}->{$sample}/$total), "\t";
	if(exists $table{$otuid}->{$sample}){
	    print sprintf("%.4f", 100 * $table{$otuid}->{$sample}/$totals{$sample}), "\t";
	}
	else{
	    print "0.0000\t";
	}
	
    }
    print $otus{$otuid}->{taxon} if(exists $otus{$otuid}->{taxon});
    print "\n";
}


#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;

my ($input1, $input2, $label);
$label = "asv";

&GetOptions("i1=s" =>\$input1,
	    "i2=s" =>\$input2,
	    "l=s" =>\$label
    );

($input1 && $input2) || 
    die "\nusage: $0\n".
    "-i1 <taxonomic annotation of each otu: *.taxonomy>\n".
    "-i2 <mothur shared file>\n".
    "-l <clustering distance 0.03|asv, default is asv>\n".
    "generate otu table summary file with taxonomic inforation\n";

my $prefix = $label eq "asv" ? "asv" : "OTU";

#OTU_1|7275 Archaea(100);Euryarchaeota(100);Methanomicrobia(100);Methanomicrobiales(100);Methanoregulaceae(100);Methanoregula(100);
open TAXON, "<$input1" or die "Could not open $input1 to read, $!\n";
my %otus = ();

while (<TAXON>) {
    chomp;
    my ($otuid, $count, $taxon) = $_ =~ /^\w+?_?(\d+?)\|(\d+?)\s+?(\w.*)$/;
    $taxon =~ s/unclassified;|_unclassified//g;
    $otus{$otuid}->{count} = $count;
    $otus{$otuid}->{taxon} = $taxon;
}
close(TAXON);

#OTUID 27 28 29
#OTU_1 2834 2553 1888

open TABLE, "<$input2" or die "Could not open $input2 to read, $!\n";
my $h = <TABLE>;
chomp($h);
my @hl = split(/\s+/, $h);
#shift the first three columns: lable, group, numOTU
shift @hl;
shift @hl;
shift @hl;

my %table = ();
my @samples = ();
my %totals = ();

while (<TABLE>) {
    chomp;
    my @line = split(/\s+/, $_);
    my $distance = shift @line;
    my $samplename = shift @line;
    my $numOtus = shift @line;
   
    my $len = @line;
    for (my $i = 0; $i < $len; $i++) {
	my $num = $line[$i];
	my $otuid = $hl[$i];
	
	#my $id = $i+1;
	$table{$otuid}->{$samplename} = $num;
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
#print "$prefix", "ID\tTotal\t". join("\t", @samples). "\t";

print "$prefix", "ID\tTotal\t". join("\t", @samples). "\t";

print "Taxonomy\n";
#foreach my $otuid (sort {$a <=> $b} keys %table){
foreach my $otuid (sort {$a <=> $b} keys %otus){
    print "$prefix\_$otuid\t" ;
    print $otus{$otuid}->{count} if(exists $otus{$otuid}->{count});
#total of each otu
    print "\t";
    
    foreach  my $sample (@samples){
	if(exists $table{"$prefix\_$otuid"}->{$sample}){
	    print $table{"$prefix\_$otuid"}->{$sample}, "\t";
	}
	else{
	    print "0\t";
	}
	
    }
    
    print $otus{"$otuid"}->{taxon} if(exists $otus{$otuid}->{taxon});
    print "\n";
}


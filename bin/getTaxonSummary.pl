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
    "-i2 <otu table file: *.OTU-table.txt>\n".
    
    "generate otu table summary file with taxonomic inforation\n";

my $prefix = $label eq "asv" ? "asv" : "OTU";
#OTU_1|7275 Archaea(100);Euryarchaeota(100);Methanomicrobia(100);Methanomicrobiales(100);Methanoregulaceae(100);Methanoregula(100);
open TAXON, "<$input1" or die "Could not open $input1 to read, $!\n";
my %otus = ();
my %taxons = ();
while (<TAXON>) {
    chomp;
    my ($otuid, $count, $taxon) = $_ =~ /^\w+?_?(\d+?)\|(\d+?)\s+?(\w.*)$/;
    $otus{$otuid}->{count} = $count;
    $taxon =~ s/unclassified;|_unclassified//g;
    $taxon =~ s/\(\d+\)//g;
#$otus{$otuid}->{taxon} = $taxon;
    push(@{$taxons{$taxon}}, $otuid);

}
close(TAXON);

#OTUID 27 28 29
#OTU_1 2834 2553 1888
#shared file
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

while (<TABLE>) {
    chomp;
    my @line = split(/\s+/, $_);
    my $distance = shift @line;
    my $samplename = shift @line;
    my $numOtus = shift @line;

    my $len = @line;

    for (my $i = 0; $i < $len; $i++) {
	my $num = $line[$i];
	#my $id = $i+1;
	my $otuid = $hl[$i];
	
	$table{$otuid}->{$samplename} = $num;

    }
    push(@samples, $samplename);
}

close(TABLE);

#header of the output
print "#Taxonomy\tTotal\t". join("\t", @samples). "\n";

foreach my $taxon (sort keys %taxons){

    print "$taxon\t";
    my $total = 0;
    my %sampletaxon = ();

    foreach my $otuid (@{$taxons{$taxon}}){
	$total += $otus{$otuid}->{count};
   #print STDERR "otuid = $otuid, total=$total;";
	foreach my $sample (@samples){
	    
	    $sampletaxon{$sample} += $table{"$prefix\_$otuid"}->{$sample} if exists $table{"$prefix\_$otuid"}->{$sample};
	}

    }

    print "$total\t";

    foreach my $sample (@samples){
	print $sampletaxon{$sample} if exists $sampletaxon{$sample};
	print "\t";
    }

    print "\n";
}



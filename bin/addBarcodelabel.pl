#!/usr/bin/env perl
use strict;
use warnings;

@ARGV == 2 or die "Usage: $0 <fastq file><relable>\n",
                  "add ;barcodelabel=samplename to each sequence headers\n";
open(FASTQ, $ARGV[0]) or die "Cannot open $ARGV[0] for reading: $!\n";
my $label = $ARGV[1];
my $head = "";
while(<FASTQ>){
    my $l1 = $_;
    my $l2 = <FASTQ>;
    my $l3 = <FASTQ>;
    my $l4 = <FASTQ>;
    chomp($l1);
    ($head) = $l1 =~ /^(@\S+)/;
    #$head =~ s/\s+//g;
    print "$head;barcodelabel=$label\n", $l2, "+\n", $l4;
}

close(FASTQ);


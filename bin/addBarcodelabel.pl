#!/usr/bin/env perl
use strict;
use warnings;

@ARGV == 2 or die "Usage: $0 <fastq file><relable>\n",
                  "add ;barcodelabel=samplename to each sequence headers\n";
open(FASTQ, $ARGV[0]) or die "Cannot open $ARGV[0] for reading: $!\n";
my $label = $ARGV[1];
my $head = "";
#print STDERR "*************************$ARGV[0]\n";
my $count = 0;
while(<FASTQ>){
    
    my $l1 = $_;
    $count++;
    my $l2 = <FASTQ>;
    $count++;
    my $l3 = <FASTQ>;
    $count++;
    my $l4 = <FASTQ>;
    $count++;
    
    $l1 ||= "";
    $l2 ||= "";
    $l3 ||= "";
    $l4 ||= "";
    
    next if $l1 eq "" || $l2 eq "" || $l3 eq "" || $l4 eq "";
    
    chomp($l1);
    
    ($head) = $l1 =~ /^(@\S+)/;
    print "$head;barcodelabel=$label\n", $l2, "+\n", $l4;
}

close(FASTQ);


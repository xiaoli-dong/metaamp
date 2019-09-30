#!/usr/bin/env perl
use warnings;
use strict;

calc_qvalue($ARGV[0], 172);
sub calc_qvalue{
    
    my ($assoFile, $total_test) = @_;
    open(ASSO, $assoFile) or die "Could not open $assoFile to read, $!\n";
    my %matrix = ();
    
    while(<ASSO>){
	next if /Significance/;
	chomp;
	my @line = split(/\t/, $_);
	$matrix{$_} = $line[3];
	
    }
    
    my $prev_sig = -1;
    my $rank = 0;

    foreach my $key (sort {$matrix{$a} <=> $matrix{$b}} keys %matrix){
	#print STDERR $key, "\n";
	my $value = $matrix{$key};
	
	if($value > $prev_sig){
	    $prev_sig = $value;
	    $rank++;
	}
	#print STDERR "$rank\n";
	#$matrix{$key}{rank} = $rank;
	
	my $qvalue = sprintf("%.6f", $rank/$total_test);
	#if($qvalue <= $sig){
	print "$key\t$rank\t$qvalue\n";
	#}
    }
    
}

#!/usr/bin/env perl

use Getopt::Long;
use warnings;
use strict;


my ($shared,$minpercent, $minnumsamples, $method);

$minpercent = 0.1;
$minnumsamples = 0;
$method = "otu";
    
GetOptions(
    "s=s" => \$shared,
    "mp=f" => \$minpercent,
    "m=s" => \$method,
    "mns=i" =>\$minnumsamples
    );

($shared) ||
    die "usage: $0 OPTIONS\n".
    "where options are:\n".
    "-s <shared file>\n".
    "-m <otu|asv>\n".
    "-mp <the minimum relative abundance in percent of an OTU in the sample, default: 0.1, it means 0.1%>\n".
    "-mns <the minimum number of samples present in an OTU, default: 0, If the number of samples present falls below the minimum, the OTU is removed>\n".
    "filter mothur shared file\n";

open(SHARED, $shared) or die "Could not open $shared file to read, $!\n";

my %otus1 = ();
my %sample2reads = ();
my @samples = ();
my %otus2 = ();
while(<SHARED>){

    chomp;

    if(/^label\s+Group/){
	next;
    }
    else{
	my @line = split(/\s+/, $_);
	my $dist = shift @line;
	my $sample = shift @line;
	my $numOtus = shift @line;
	push(@samples, $sample);
	my $index = 0;
	push(@{$otus1{$dist}->{$sample}}, @line);
	my $sum = 0;
	map { $sum += $_ } @line;

	$sample2reads{$dist}->{$sample} = $sum;
	for(my $i =0; $i < $numOtus; $i++){
	    push(@{$otus2{$dist}->{$i}}, $line[$i]);
	}
	
    }
}
close(SHARED);

my %keep_ids = ();
for my $dist (keys %otus2){
    my $keep_ids = ();

   for my $otuid (keys %{$otus2{$dist}}){

	my @col = @{$otus2{$dist}->{$otuid}};
	my $count = 0;
	for (my $i = 0; $i < @samples; $i++){
	    my $total = $sample2reads{$dist}->{$samples[$i]};
	    if(100*$col[$i]/$total >= $minpercent){
		$count++;
		if($count >= $minnumsamples){
		    last;
		}
	    }
	    
	}
	if($count >= $minnumsamples){
	    $keep_ids{$dist}->{$otuid} = 1;
	}
	
   }
}

for my $dist (keys %otus1){
    my @ids = sort {$a <=>$b} keys %{$keep_ids{$dist}};

    my $numOtus = @ids;
    my $prefix = $method eq "otu" ? "otu" : "asv";
    print "label\tGroup\tnum$prefix", "s\t", join("\t", map{"$prefix\_" . ($_+1)} @ids), "\n";
    
    for my $sample (keys %{$otus1{$dist}}){
	print "$dist\t$sample\t$numOtus";
	
	my @row = @{$otus1{$dist}->{$sample}};
	
	for(my $i = 0; $i < @ids; $i++){
	    my $id = $ids[$i];
	    print "\t", $row[$id];
	}
	print "\n";
	
	
    }
}

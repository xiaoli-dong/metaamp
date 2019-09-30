#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Long;

my ($uc);

&GetOptions("uc=s" =>\$uc
    );

($uc) or 
    die "\nusage: $0 \n".
    "-uc <otu mapping file in uc format file>\n".
    "read in mapping file in uc format, split the reads among the samples and  generate the otu table file\n";

open UC, "$uc" or die "Could not open $uc, $!\n";

my %group2otusize = ();
my %otu2group2size = ();
while (<UC>) {
    if(/^H/){
	my @line = split(/\t/, $_);
	my ($otuid) = $line[9] =~ /OTU_(\d+)/;
	my($seqid, $groupid) = $line[8] =~ /^(\S+?);\w+?=(\S+)/;
	$group2otusize{$groupid}->{$otuid}++;
	$otu2group2size{$otuid}->{$groupid}++;
    }
}

close(UC);

my @groups = sort keys %group2otusize;
print "OTUID\t", join("\t", @groups), "\n";

foreach my $id (sort {$a <=>$b} keys %otu2group2size){
    print "OTU_$id";
    
    foreach my $group (@groups){
	
	if(exists  $otu2group2size{$id}->{$group}){
	    print "\t", $otu2group2size{$id}->{$group};
	}
	else{

	    print "\t0";
	}
    }
    print "\n";
}


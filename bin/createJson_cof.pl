#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
my ($csv, $sep, $prefix);
$sep = "\t";
$prefix = "";

&GetOptions("c=s" =>\$csv,
	    "s=s" =>\$sep,
	    "p=s" =>\$prefix
	    
    );
($csv) or
    die "\nusage: $0 \n".
    "-c <csv file>\n".
    "-p <array name prefix, default is \"\">\n".
    "-s <column separator, default is \"\t\">\n" .
    "-ex <column index excluded, the first col is 0: 0:1, default is \"\">\n";

die "csv file $csv does not exist, $!\n" if (! -e $csv);
#perl ../../../../bin/createJson_heatmap.pl -c mock.shared -s "\t" -ex "0:2"

open(CSV, $csv) or die "Could not open $csv to read, $!\n";

my @cols = ();
my $l = <CSV>;
chomp($l);
my @header = split($sep, $l);

#write data variable
print $prefix, "data = [\n";

my $row_count = 0;
my $str = "";
my %rows = ();

my  %data = ();

while(<CSV>){
    chomp;
    my @l = split($sep, $_);
    my $rid = $l[0];
    my $cid = $l[1];
    my $cof = $l[2];
    my $sig = $l[3];
    $rows{$rid}++;
    $rows{$cid}++;
    $data{$rid}->{$cid}->{cof} = $cof;
    $data{$cid}->{$rid}->{cof} = $cof;
    $data{$rid}->{$cid}->{sig} = $sig;
    $data{$cid}->{$rid}->{sig} = $sig;
}

foreach my $rid (sort keys %rows){
    
    foreach my $cid (sort keys %rows){

	if($rid eq $cid){
	    $str .= (' ' x 2) . "{";
	    $str .= "\"rowid\": \"$rid\", ";
	    $str .= "\"colid\": \"$cid\", ";
	    $str .= "\"value\": \"1\", ";
	    $str .= "\"sig\": 0";
	    $str .= "},\n"; 
	}
	else{
	    #my $cof = exists $data{$rid}->{$cid}->{cof} ? $data{$rid}->{$cid}->{cof}: 0;
	    #my $sig = exists $data{$rid}->{$cid}->{sig} ? $data{$rid}->{$cid}->{sig}: 999;
	    if(exists  $data{$rid}->{$cid}->{cof}){
		my $cof = $data{$rid}->{$cid}->{cof};
		my $sig = $data{$rid}->{$cid}->{sig};
		$str .= (' ' x 2) . "{";
		$str .= "\"rowid\": \"$rid\", ";
		$str .= "\"colid\": \"$cid\", ";
		$str .= "\"value\": \"$cof\", ";
		$str .= "\"sig\": $sig";
		$str .= "},\n";
	    }
	    else{
		
		my $cof = "0";
		my $sig = "999";
		$str .= (' ' x 2) . "{";
		$str .= "\"rowid\": \"$rid\", ";
		$str .= "\"colid\": \"$cid\", ";
		$str .= "\"value\": \"$cof\", ";
		$str .= "\"sig\": $sig";
		$str .= "},\n";
	    }
	}
	
    }
}



$str =~ s/,\n$//;
print "$str\n];\n";


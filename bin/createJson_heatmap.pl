#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
my ($csv, $sep, $prefix, $exclude);
$sep = ",";
$prefix = "";

&GetOptions("c=s" =>\$csv,
	    "s=s" =>\$sep,
	    "p=s" =>\$prefix,
	    "ex=s" =>\$exclude
    );
($csv) or
    die "\nusage: $0 \n".
    "-c <csv file>\n".
    "-p <array name prefix, default is \"\">\n".
    "-s <column separator, default is \",\">\n" .
    "-ex <column index excluded, the first col is 0: 0:1, default is \"\">\n";

die "csv file $csv does not exist, $!\n" if (! -e $csv);
#perl ../../../../bin/createJson_heatmap.pl -c mock.shared -s "\t" -ex "0:2"

my %excludes = ();
foreach (split(":", $exclude)){
    $excludes{$_} = 1;
}
open(CSV, $csv) or die "Could not open $csv to read, $!\n";

my @cols = ();
my $l = <CSV>;
chomp($l);
my @header = split($sep, $l);

foreach my $i (0 .. $#header){
    if(not exists $excludes{$i}){
	push(@cols, $header[$i]);
    }
}

shift(@cols);


#write data variable
print $prefix, "data = [\n";

my $row_count = 0;
my $str = "";
while(<CSV>){
    chomp;
    my @new_l = ();
    my @l = split($sep, $_);
    foreach my $i (0 .. $#l){
	if(not exists $excludes{$i}){
	    push(@new_l, $l[$i]);
	}
    }
    
    my $rowid = shift @new_l;
    
    $row_count++;
    my $col_count = 0;
    
    for (my $i=0; $i < @cols; $i++) {
	$col_count++;
	$str .= (' ' x 2) . "{";
	#$str = ",\n" . $str if  $col_count > 1;
	
	my $col = $cols[$i];
	my $value = $new_l[$i];
	$str .= " " if($i > 0);
	
	$str .= "\"colid\": \"$col\", ";
	$str .= "\"rowid\": \"$rowid\", ";
	$str .= "\"value\": \"$value\"";
	$str .= "},\n";
	
    }
    
}
$str =~ s/,\n$//;
print "$str\n];\n";


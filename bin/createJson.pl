#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
my ($csv, $sep, $prefix);
$sep = ",";
$prefix = "";

&GetOptions("c=s" =>\$csv,
	    "s=s" =>\$sep,
	    "p=s" =>\$prefix
    );
($csv) or
    die "\nusage: $0 \n".
    "-c <csv file>\n".
    "-p <array name prefix, default is \"\">\n".
    "-s <column separator, default is \",\">\n";

die "csv file $csv does not exist, $!\n" if (! -e $csv);

open(CSV, $csv) or die "Could not open $csv to read, $!\n";

my @header = ();
my $l = <CSV>;
chomp($l);
$l =~ s/^\#?//g;
@header = split($sep, $l);

#print table header
my $str = "";
print $prefix, "tablecol = [\n";
my $col_count = 0;

foreach my $col (@header){
    
    $col =~ s/0\.03-//g;
    $col_count++;
    my $str = (' ' x 2) . "{";
    $str = ",\n" . $str if  $col_count > 1;
    $str .= "\"title\": \"$col\",";
    $str .= " \"data\": \"$col\"";
    $str .= "}";
    print $str;
}

print "\n];\n";



#write data variable
print $prefix, "data = [\n";

my $row_count = 0;

while(<CSV>){
    next if !/\S/;
    chomp;
    my @l = split($sep, $_);
    $row_count++;
    my $str = (' ' x 2) . "{";
    $str = ",\n" . $str if  $row_count > 1;
    
    
    for (my $i=0; $i < @header; $i++) {
	my $h = $header[$i];
	my $value = $l[$i];
	$str .= " " if($i > 0);
	
	$str .= "\"$h\": \"$value\"";
	$str .= "," if $i < @header -1
    }
    $str .= "}";
    print $str;
}
print "\n];\n";


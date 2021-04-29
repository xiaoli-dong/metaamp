#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
my ($csv, $prefix);

$prefix = "";

&GetOptions("c=s" =>\$csv,
	    "p=s" =>\$prefix
    );
($csv) or
    die "\nusage: $0 \n".
    "-c <csv file>\n".
    "-p <array name prefix, default is \"\">\n";
    

die "csv file $csv does not exist, $!\n" if (! -e $csv);
#perl ../../../../bin/createJson_heatmap.pl -c mock.shared

open(CSV, $csv) or die "Could not open $csv to read, $!\n";

my @cols = ();
my $l = <CSV>;
chomp($l);

@cols = split(/\s+/, $l);
shift(@cols);
shift(@cols);

#write data variable
print $prefix, "data = [\n";

my $row_count = 0;
my $str = "";
while(<CSV>){
    chomp;
    my @l = split(/\s+/, $_);
    shift(@l);
    my $rowid = shift(@l);
    my $colid = shift(@l);
    
    #print STDERR "****************$rowid, $colid\n";
    
    
    $row_count++;
    my $col_count = 0;
    $str .= (' ' x 2) . "{";
    $str .= "\"colid\": \"$colid\", ";
    $str .= "\"rowid\": \"$rowid\", ";
    for (my $i=0; $i < @cols; $i++) {
	my $colhead = $cols[$i];
	my $value = $l[$i];
	if($i == @cols -1){
	    $str .= "\"$colhead\": \"$value\"";
	}
	else{
	    $str .= "\"$colhead\": \"$value\",";
	}
    }
    $str .= "},\n";
}
$str =~ s/,\n$//;
print "$str\n];\n";


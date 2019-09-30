#!/usr/bin/env perl
use warnings;
use strict;

$#ARGV == -1
    and die "Usage: $0 <OTU associate output><shared file>\n";

open(ASSO, $ARGV[0]) or die "Could not open $ARGV[0] to read, $!\n";


my %edges = ();
my %all_otus = ();
while(<ASSO>){
    next if /Significance/;
    chomp;
    my @line = split(/\t/, $_);
    $edges{$line[0]}->{$line[1]}->{weight} = $line[2];
    $edges{$line[0]}->{$line[1]}->{sig} = $line[3];
    if(not exists $all_otus{$line[0]}){
	$all_otus{$line[0]} = 1;
    }
    if(not exists $all_otus{$line[1]}){
	$all_otus{$line[1]} = 1;
    }
}
close(ASSO)


open(SHARE, $ARGV[1]) or die "Could not open $ARGV[1] to read, $!\n";
my %otus = ();
while(<SHARE>){
    next if /label/;
    chomp;
    my @line = split(/\s+/, $_);
    my $dist = shift @line;
    my $sample = shift @line;
    shift @line;
    my $count = 1;
    foreach my $otu (@line){
	$otus{"Otu$count"} = $out;
	$count++;
    }
}
close(SHARE)


print "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
print "<gexf xmlns=\"http://www.gexf.net/1.2draft\" version=\"1.2\">\n";
print "<graph mode=\"static\" defaultedgetype=\"directed\">\n";
print "<attributes class=\"edge\">\n";
print "<attribute id=\"0\" title=\"Significance\" type=\"float\"/>\n";
print "</attribute>\n";

#nodes section
pint "<nodes>\n";
foreach my $otu (sort keys %all_otus){
    my $size = $otus{$otu};
    print "<node id=\"$otu\">\n";
    print "<vis:size value=\"$size\" />\n";
    print "</node>\n";
}
print "</nodes>\n";

#edge section
pint "<edges>\n";
my $edge_count = 0;
foreach my $otuA (sort keys %edges){
    
    foreach my $otuB (sort keys %{$edges{$otuA}}){
	
	my $size = $edges{$otuA}->{$otuB}->{weight};
	my $sig = $edges{$otuA}->{$otuB}->{sig};
	print "<edge id=\"$edge_count\" source=\"$outA\" target=\"$otuB\">\n";
	print "<attvalues><attvalue for=\"0\" value=\"$sig\"></attvalues>\n";
	print "<vis:size value=\"$size\" />\n";
	print "</edge>\n";
	$edge_count++;
    }
}
print "</edges>\n";

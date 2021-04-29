#!/usr/bin/env perl
use Getopt::Long;
use strict;
use File::Slurp;


my ($input, $template, $distance);
$distance = "asv";

&GetOptions("i=s"       => \$input,
	    "t=s" =>\$template,
	    "d=s" =>\$distance
    );

($input && $template) ||
    die "usage: $0 OPTIONS

where options are: -i  <*.groups.summary.json file generated from mothur analysis> <-t template for generate rarefaction curve> <-d custering distance or dissimilarity: asv|0.03, default is asv>\n";

my $type = "asv";
$type = "OTU" if $distance ne "asv";

my $template_text = read_file($template);

$template_text =~ s/ebgmetaampjson/..\/$input/g;
$template_text  =~ s/ebgmetaamptype/$type/g;
print $template_text;

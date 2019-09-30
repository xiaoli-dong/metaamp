#!/usr/bin/env perl
use Getopt::Long;
use strict;
use File::Slurp;


my ($input, $template);
&GetOptions("i=s"       => \$input,
	    "t=s" =>\$template
    );

($input && $template) ||
    die "usage: $0 OPTIONS

where options are: -i  <*.groups.summary file generated from mothur analysis> <-t template for generate rarefaction curve>\n";

my $template_text = read_file($template);

$template_text =~ s/ebgmetaamptsv/..\/$input/;

print $template_text;

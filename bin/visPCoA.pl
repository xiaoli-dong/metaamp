#!/usr/bin/env perl
use Getopt::Long;
use strict;
use File::Slurp;


my ($axesfile, $loadingfile, $template);
&GetOptions("a=s"       => \$axesfile,
	    "l=s" => \$loadingfile,
	    "t=s" =>\$template
    );

($axesfile && $loadingfile && $template) ||
    die "usage: $0 OPTIONS

where options are: -a  <*.pcoa.axes> -l <*.pcoa.loadings> <-t template for generate rarefaction curve>\n";

my $template_text = read_file($template);

$template_text =~ s/ebgmetaamppcoaaxes/..\/$axesfile/;
$template_text =~ s/ebgmetaamppcoaloading/..\/$loadingfile/;
print $template_text;

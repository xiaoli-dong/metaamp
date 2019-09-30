#!/usr/bin/env perl
use Getopt::Long;
use strict;
use File::Slurp;


my ($axesfile, $stressfile, $template);
&GetOptions("s=s"       => \$stressfile,
	    "a=s" => \$axesfile,
	    "t=s" =>\$template
    );

($axesfile && $stressfile && $template) ||
    die "usage: $0 OPTIONS

where options are: -a  <*.nmds.stresss> -l <*.nmds.axes> <-t template for generate rarefaction curve>\n";

my $template_text = read_file($template);

$template_text =~ s/ebgmetaampstress/..\/$stressfile/;
$template_text =~ s/ebgmetaampaxes/..\/$axesfile/;
print $template_text;

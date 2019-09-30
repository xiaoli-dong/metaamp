#!/usr/bin/env perl

use warnings;
use strict;
use Getopt::Long;
use File::Path qw(make_path remove_tree);
use File::Copy;
use FindBin;
use File::Copy::Recursive qw(dircopy);

my $analysisName;

GetOptions("an=s" => \$analysisName);
$analysisName or
    die "usage: $0 OPTIONS where options are: -an analysisName\n";

#umask 022;
my $bin = "$FindBin::RealBin/";
my $js = "$bin/../js";
my $css = "$bin/../css";
my $images = "$bin/../images";

my $outDir = "results";
remove_tree $outDir if -d $outDir;
make_path($outDir);

my @subdirs = (
    "$outDir/OTU_and_taxonomy",
    "$outDir/alpha_and_beta_diversity",
    "$outDir/hypothesis_testing",
    "$outDir/QC",
    "$outDir/html"
    );

dircopy($js,"$outDir/js") or die "Copy $js to $outDir/js failed, $!\n";
dircopy($css,"$outDir/css") or die "Copy $css to $outDir/css failed, $!\n";
dircopy($images,"$outDir/images") or die "Copy $js to $outDir/images failed, $!\n";

make_path (@subdirs);

my $file;

########## OTUs##########

copy "$analysisName.OTUs.fasta", $subdirs[0];
for $file (
    (glob "*.OTU-table.taxonomy"),
    (glob "*.otu"),
    (glob "*.taxonomy.*summary"),
    (glob "$analysisName.OTUs.*.taxonomy")
    ){
    copy $file,$subdirs[0];
}
########## Alpha_diversity ##########
copy "$analysisName.summary",$subdirs[1];
for $file (
    #(glob "*.rabund"),
    (glob "*.groups.rarefaction"),
    (glob "*.groups.summary"),
    (glob "*.corr"),
    (glob "*.pcoa.*"),
    (glob "*.tre"),
    (glob "*.nmds.*"),
    (glob "*.list"),
    (glob "*.groups"),
    (glob "*.shared"),
    (glob "*relabund")

    ){
    copy $file, $subdirs[1];
}

for $file (
    (glob "*.parsimony"),
    (glob "*.psummary"),
    (glob "*.weighted"),
    (glob "*.trewsummary"),
    (glob "*.unweighted"),
    (glob "*.uwsummary"),
    (glob "*.amova"),
    (glob "*.homova")

    ){#print STDERR "$file\n";
    copy $file,$subdirs[2];
}

copy "$analysisName.qc.summary",$subdirs[3];

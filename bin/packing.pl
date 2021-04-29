#!/usr/bin/env perl

use warnings;
use strict;
use Getopt::Long;
use File::Path qw(make_path remove_tree);
use File::Copy;
use FindBin;
use File::Copy::Recursive qw(dircopy);

my ($analysisName, $label);
$label = "asv";

GetOptions("an=s" => \$analysisName,
	   "l=s" =>\$label
    );
$analysisName or
    die "usage: $0 OPTIONS where options are: -an analysisName -l <clustering distance, asv|0.03|0.05..., default is asv>\n";

my $bin = "$FindBin::RealBin/";
my $js = "$bin/../js";
my $css = "$bin/../css";
my $images = "$bin/../images";
my $template = "$bin/../template";

my $outdir = "results";
remove_tree $outdir if -d $outdir;
make_path($outdir);

my @subdirs = ();

$label eq "asv" ? push(@subdirs, "$outdir/asv_and_taxonomy") : push(@subdirs, "$outdir/otu_and_taxonomy");
push(@subdirs, ("$outdir/alpha_and_beta_diversity","$outdir/hypothesis_testing","$outdir/qc","$outdir/html"));
dircopy($js,"$outdir/js") or die "Copy $js to $outdir/js failed, $!\n";
dircopy($css,"$outdir/css") or die "Copy $css to $outdir/css failed, $!\n";
dircopy($images,"$outdir/images") or die "Copy $js to $outdir/images failed, $!\n";
make_path (@subdirs);

my $file;

########## OTUs/asv##########
$label eq "asv" ? copy "$analysisName.asv.fasta", $subdirs[0] : copy "$analysisName.otu.fasta", $subdirs[0];

foreach (glob "$analysisName*taxonomy*"){
    copy $_, $subdirs[0] or die "Copy file=$_ failed, $!\n";
}

########## QC ##########
#copy "$analysisName.qc.summary",$subdirs[3] or die "Copy file=$analysisName.qc.summary failed, $!\n";

foreach (glob "*.stats.*"){
    copy $_ ,$subdirs[3] or die "Copy file=$_ failed, $!\n";
}
copy "$template/taxontree.html",$subdirs[4] or die "Copy file=$template/taxontree.html failed, $!\n";;
copy "$template/heatmap.js",$subdirs[4] or die "Copy file=$template/heatmap.js failed, $!\n";;
copy "$template/heatmap_relabund.js",$subdirs[4] or die "Copy file=$template/heatmap_relabund.js failed, $!\n";;
copy "$template/heatmap_network.js",$subdirs[4] or die "Copy file=$template/heatmap_network.js failed, $!\n";;
copy "$template/pca.js",$subdirs[4] or die "Copy file=$template/pca.js failed, $!\n";
copy "$template/nmds.js",$subdirs[4] or die "Copy file=$template/nmds.js failed, $!\n";
copy "$template/alpha_diversity_indexes.js",$subdirs[4] or die "Copy file=$template/alpha_diversity_indexes.js failed, $!\n";
copy "$template/bubble_plot.js",$subdirs[4] or die "Copy file=$template/bubble_plot.js failed, $!\n";
copy "$template/cor_bubble.js",$subdirs[4] or die "Copy file=$template/cor_bubble.js failed, $!\n";
copy "$template/summary_heatmap.js",$subdirs[4] or die "Copy file=$template/summary_heatmap.js failed, $!\n";
########## Alpha_beta diversity ##########
copy "$analysisName.summary",$subdirs[1] or die "Copy file=$analysisName.summary failed, $!\n";;

for $file (
    
    #(glob "*.rabund"),
    (glob "*.groups.rarefaction"),
    (glob "*.groups.summary"),
    (glob "*.corr"),
    (glob "*.pcoa.*"),
    (glob "*.tre"),
    (glob "*.nmds.*"),
    (glob "*.groups"),
    (glob "*.shared"),
    (glob "*relabund")

    ){
    copy $file, $subdirs[1]  or die "Copy file=$file failed, $!\n";;
}

#hypothesis test
for $file (
    (glob "*.parsimony"),
    (glob "*.psummary"),
    (glob "*.weighted"),
    (glob "*.wsummary"),
    (glob "*.unweighted"),
    (glob "*.uwsummary"),
    (glob "*.amova"),
    (glob "*.homova")

    ){#print STDERR "$file\n";
    copy $file,$subdirs[2] or die "Copy file=$file failed, $!\n";
}


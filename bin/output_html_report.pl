#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Long;
use FindBin;
use lib "$FindBin::Bin";
require "util.pl";
my($analysis, $seqtype, $inputdir,$distance);

$seqtype = "paired";
$inputdir = "./results";
my $cpus = 8;
$distance = "asv";

&GetOptions("an=s" =>\$analysis,
	    "seqtype=s" =>\$seqtype,
	    "dir=s" =>\$inputdir,
	    "d=s" =>\$distance
    );

($analysis && $seqtype) or
    die "usage: $0 \n".
    "-an <one word analysis name>\n".
    "-seqtype <sequence type single|paired, default is paired>\n".
    "-d <custering distance or dissimilarity: asv|0.03>\n".
    "-dir <input directory for the metaamp packed results>\n";


my $bin = "$FindBin::RealBin";
my $template_dir = "$bin/../template";
my $type = "asv";

$type = "otu" if $distance ne "asv";

#umask 022;

open(JOB_INFO, "JobInfo.txt") or die "Can't open job information file for reading: $!";

my($email, $human_time, $machine_time, $command);

while(<JOB_INFO>){
    chomp;
    if(/Email:(\S+)/){
	$email = $1;
    }
    elsif(/Submit human time:(\S.*)/){
	$human_time = $1;
    }
    elsif(/Submit machine time:(\S.*)/){
	$machine_time = $1;
    }
    elsif(/Command:(\S.*)/){
	$command = $1;
	$command =~ s/.*?bin\///;
    }
}
close JOB_INFO;


chdir $inputdir;
vis_otu_taxon($analysis);
vis_alpha_beta_diversity($analysis);
use File::Slurper 'read_text';
my $index_str = read_text("$template_dir/index.html");
$index_str =~ s/METAAMP_ANALYSISNAME/$analysis/g;
$index_str =~ s/METAAMP_INPUT_COMMAND/$command/g;
$index_str =~ s/METAAMP_USER_EMAIL/$email/g;
$index_str =~ s/METAAMP_SUBMIT_TIME/$human_time/g;
$index_str =~ s/METAAMP_LABEL/$distance/g;
$index_str =~ s/METAAMP_TYPE/$type/g;
$index_str =~ s/!--METAAMP_MERGE_DATA_|_METAAMP_MERGE_DATA--//g if $distance ne "asv" && $seqtype ne "single";
$index_str =~ s/\/\* METAAMP_MERGE_SCRIPT|METAAMP_MERGE_SCRIPT \*\///g if $distance ne "asv" && $seqtype ne "single";
$index_str =~ s/METAAMP_MERGE_HTML--|\!--METAAMP_MERGE_HTML//g if $distance ne "asv" && $seqtype ne "single";

open(OUT, ">index.html") or die "Could not open index.html to write, $!\n";
binmode(OUT, ":utf8");
print OUT $index_str;
close(OUT);


chdir "..";
my $cmd = "zip -q results.zip -r results";
runcmds($cpus, $cmd);



sub vis_otu_taxon{
    my ($analysis) = @_;
    
    my $cmd = $type eq "otu" ? "$bin/output_tree_json.pl  -i otu_and_taxonomy/$analysis.taxonomy.summary.txt -t Taxonomy -l 3 > html/taxon.tree.json;" : "$bin/output_tree_json.pl  -i asv_and_taxonomy/$analysis.taxonomy.summary.txt -t Taxonomy -l 3 > html/taxon.tree.json;";
    runcmds($cpus, $cmd);
}


sub vis_alpha_beta_diversity{

    my($analysis) = @_;
    my $cmd = "$bin/createJson.pl -c alpha_and_beta_diversity/$analysis.groups.rarefaction -s \"\\t\" > alpha_and_beta_diversity/$analysis.groups.rarefaction.json;";
    $cmd .= "$bin/visJson.pl -i alpha_and_beta_diversity/$analysis.groups.rarefaction.json -t $template_dir/rarefaction.html -d $distance> html/$analysis.groups.rarefaction.html;";
    $cmd .= "$bin/createJson.pl -c alpha_and_beta_diversity/$analysis.subsample.groups.rarefaction -s \"\\t\" > alpha_and_beta_diversity/$analysis.subsample.groups.rarefaction.json;";
    $cmd .= "$bin/visJson.pl -i alpha_and_beta_diversity/$analysis.subsample.groups.rarefaction.json -t $template_dir/rarefaction.html -d $distance > html/$analysis.subsample.groups.rarefaction.html;";
    runcmds($cpus, $cmd);

    $cmd = "$bin/createJson.pl -c alpha_and_beta_diversity/$analysis.groups.summary -s \"\\t\" > alpha_and_beta_diversity/$analysis.groups.summary.json;";
    $cmd .= "$bin/visJson.pl -i alpha_and_beta_diversity/$analysis.groups.summary.json -t $template_dir/alpha_diversity_indexes.html -d $distance > html/$analysis.groups.summary.html;";
    $cmd .= "$bin/createJson.pl -c  alpha_and_beta_diversity/$analysis.subsample.groups.summary -s \"\\t\" > alpha_and_beta_diversity/$analysis.subsample.groups.summary.json;";
    $cmd .= "$bin/visJson.pl -i alpha_and_beta_diversity/$analysis.subsample.groups.summary.json -t $template_dir/alpha_diversity_indexes.html -d $distance > html/$analysis.subsample.groups.summary.html";
    runcmds($cpus, $cmd);
    
    #############shared file
    $cmd = "$bin/createJson_heatmap.pl -c alpha_and_beta_diversity/$analysis.shared -s \"\\t\" -ex \"0:2\"> alpha_and_beta_diversity/$analysis.shared_heatmap.json;";
    $cmd .= "$bin/visJson.pl -i alpha_and_beta_diversity/$analysis.shared_heatmap.json -t $template_dir/shared_heatmap.html -d $distance > html/$analysis.shared_heatmap.html;";
    $cmd .= "$bin/createJson_heatmap.pl -c alpha_and_beta_diversity/$analysis.subsample.shared -s \"\\t\" -ex \"0:2\"> alpha_and_beta_diversity/$analysis.subsample.shared_heatmap.json;";
    $cmd .= "$bin/visJson.pl -i alpha_and_beta_diversity/$analysis.subsample.shared_heatmap.json -t $template_dir/shared_heatmap.html -d $distance > html/$analysis.subsample.shared_heatmap.html;";
    runcmds($cpus, $cmd);

    ###########summary file
    $cmd = "$bin/summary2json_heatmap.pl -c alpha_and_beta_diversity/$analysis.summary > alpha_and_beta_diversity/$analysis.summary.json;";
    $cmd .= "$bin/visJson.pl -i alpha_and_beta_diversity/$analysis.summary.json -t $template_dir/summary_heatmap.html -d $distance > html/$analysis.summary_heatmap.html;";
    runcmds($cpus, $cmd);
    
    ##OTU relabund file
    $cmd = "$bin/createJson_heatmap.pl -c alpha_and_beta_diversity/$analysis.relabund -s \"\\t\" -ex \"0:2\"> alpha_and_beta_diversity/$analysis.relabund_heatmap.json;";
    $cmd .= "$bin/visJson.pl -i alpha_and_beta_diversity/$analysis.relabund_heatmap.json -t $template_dir/relabund_heatmap.html -d $distance > html/$analysis.relabund_heatmap.html;";
    $cmd .= "$bin/createJson_heatmap.pl -c alpha_and_beta_diversity/$analysis.subsample.relabund -s \"\\t\" -ex \"0:2\"> alpha_and_beta_diversity/$analysis.subsample.relabund_heatmap.json;";
    $cmd .= "$bin/visJson.pl -i alpha_and_beta_diversity/$analysis.subsample.relabund_heatmap.json -t $template_dir/relabund_heatmap.html -d $distance > html/$analysis.subsample.relabund_heatmap.html;";
    runcmds($cpus, $cmd);

    ##OTU relabund file bubble_plot
    $cmd = "$bin/visJson.pl -i alpha_and_beta_diversity/$analysis.relabund_heatmap.json -t $template_dir/relabund_bubble_plot.html -d $distance > html/$analysis.relabund_bubble_plot.html;";
    $cmd .= "$bin/visJson.pl -i alpha_and_beta_diversity/$analysis.subsample.relabund_heatmap.json -t $template_dir/relabund_bubble_plot.html -d $distance > html/$analysis.subsample.relabund_bubble_plot.html;";
    runcmds($cpus, $cmd);
    
    #network
    #$cmd = "$bin/filter_association.pl -a alpha_and_beta_diversity/$analysis.$distance.spearman.corr > alpha_and_beta_diversity/$analysis.$distance.spearman.filtered.corr;";
    #$cmd .= "$bin/createJson_cof.pl -c alpha_and_beta_diversity/$analysis.$distance.spearman.filtered.corr > alpha_and_beta_diversity/$analysis.$distance.spearman.filtered.json;";
    #$cmd .= "$bin/visJson.pl -i alpha_and_beta_diversity/$analysis.$distance.spearman.filtered.json -t $template_dir/spearman_network_heatmap.html -d $distance > html/$analysis.$distance.spearman_heatmap.html;";
    #runcmds($cpus, $cmd);
    
    #$cmd = "$bin/filter_association.pl -a alpha_and_beta_diversity/$analysis.$distance.pearson.corr > alpha_and_beta_diversity/$analysis.$distance.pearson.filtered.corr;";
    #$cmd .= "$bin/createJson_cof.pl -c alpha_and_beta_diversity/$analysis.$distance.pearson.filtered.corr > alpha_and_beta_diversity/$analysis.$distance.pearson.filtered.json;";
    #$cmd .= "$bin/visJson.pl -i alpha_and_beta_diversity/$analysis.$distance.pearson.filtered.json -t $template_dir/pearson_network_heatmap.html -d $distance > html/$analysis.$distance.pearson_heatmap.html;";
    #runcmds($cpus, $cmd);

    #$cmd = "$bin/filter_association.pl -a alpha_and_beta_diversity/$analysis.$distance.kendall.corr > alpha_and_beta_diversity/$analysis.$distance.kendall.filtered.corr;";
    #$cmd .= "$bin/createJson_cof.pl -c alpha_and_beta_diversity/$analysis.$distance.kendall.filtered.corr > alpha_and_beta_diversity/$analysis.$distance.kendall.filtered.json;";
    #$cmd .= "$bin/visJson.pl -i alpha_and_beta_diversity/$analysis.$distance.kendall.filtered.json -t $template_dir/kendall_network_heatmap.html -d $distance > html/$analysis.$distance.kendall_heatmap.html;";
    #runcmds($cpus, $cmd);

    #network bubble plot
    
    #$cmd = "$bin/visJson.pl -i alpha_and_beta_diversity/$analysis.$distance.spearman.filtered.json -t $template_dir/spearman_network_bubble_plot.html -d $distance > html/$analysis.$distance.spearman_bubble_plot.html;";
    #$cmd .= "$bin/visJson.pl -i alpha_and_beta_diversity/$analysis.$distance.pearson.filtered.json -t $template_dir/pearson_network_bubble_plot.html -d $distance > html/$analysis.$distance.pearson_bubble_plot.html;";
    #$cmd .= "$bin/visJson.pl -i alpha_and_beta_diversity/$analysis.$distance.kendall.filtered.json -t $template_dir/kendall_network_bubble_plot.html  -d $distance > html/$analysis.$distance.kendall_bubble_plot.html;";
    #runcmds($cpus, $cmd);
    
    #network subsample
    #$cmd = "$bin/filter_association.pl -a alpha_and_beta_diversity/$analysis.subsample.$distance.spearman.corr > alpha_and_beta_diversity/$analysis.subsample.$distance.spearman.filtered.corr;";
    #$cmd .= "$bin/createJson_cof.pl -c alpha_and_beta_diversity/$analysis.subsample.$distance.spearman.filtered.corr > alpha_and_beta_diversity/$analysis.subsample.$distance.spearman.filtered.json;";
    #$cmd .= "$bin/visJson.pl -i alpha_and_beta_diversity/$analysis.subsample.$distance.spearman.filtered.json -t $template_dir/spearman_network_heatmap.html -d $distance > html/$analysis.subsample.$distance.spearman_heatmap.html;";
    #runcmds($cpus, $cmd);

    #$cmd = "$bin/filter_association.pl -a alpha_and_beta_diversity/$analysis.subsample.$distance.pearson.corr > alpha_and_beta_diversity/$analysis.subsample.$distance.pearson.filtered.corr;";
    #$cmd .= "$bin/createJson_cof.pl -c alpha_and_beta_diversity/$analysis.subsample.$distance.pearson.filtered.corr > alpha_and_beta_diversity/$analysis.subsample.$distance.pearson.filtered.json;";
    #$cmd .= "$bin/visJson.pl -i alpha_and_beta_diversity/$analysis.subsample.$distance.pearson.filtered.json -t $template_dir/pearson_network_heatmap.html -d $distance > html/$analysis.subsample.$distance.pearson_heatmap.html;";
    #runcmds($cpus, $cmd);
    
    #$cmd = "$bin/filter_association.pl -a alpha_and_beta_diversity/$analysis.subsample.$distance.kendall.corr > alpha_and_beta_diversity/$analysis.subsample.$distance.kendall.filtered.corr;";
    #$cmd .= "$bin/createJson_cof.pl -c alpha_and_beta_diversity/$analysis.subsample.$distance.kendall.filtered.corr > alpha_and_beta_diversity/$analysis.subsample.$distance.kendall.filtered.json;";
    #$cmd .= "$bin/visJson.pl -i alpha_and_beta_diversity/$analysis.subsample.$distance.kendall.filtered.json -t $template_dir/kendall_network_heatmap.html -d $distance > html/$analysis.subsample.$distance.kendall_heatmap.html;";
    
    #runcmds($cpus, $cmd);

    #network subsample bubble plot

    
    #$cmd = "$bin/visJson.pl -i alpha_and_beta_diversity/$analysis.subsample.$distance.spearman.filtered.json -t $template_dir/spearman_network_bubble_plot.html -d $distance > html/$analysis.subsample.$distance.spearman_bubble_plot.html;";
    #$cmd .= "$bin/visJson.pl -i alpha_and_beta_diversity/$analysis.subsample.$distance.pearson.filtered.json -t $template_dir/pearson_network_bubble_plot.html -d $distance > html/$analysis.subsample.$distance.pearson_bubble_plot.html;";
    #$cmd .= "$bin/visJson.pl -i alpha_and_beta_diversity/$analysis.subsample.$distance.kendall.filtered.json -t $template_dir/kendall_network_bubble_plot.html -d $distance > html/$analysis.subsample.$distance.kendall_bubble_plot.html;";
    #runcmds($cpus, $cmd);


    #######################
    $cmd = "$bin/createJson.pl -c alpha_and_beta_diversity/$analysis.braycurtis.$distance.lt.pcoa.axes -s \"\\t\" -p axes > alpha_and_beta_diversity/$analysis.braycurtis.$distance.lt.pcoa.axes.json;";
    $cmd .= "$bin/createJson.pl -c alpha_and_beta_diversity/$analysis.braycurtis.$distance.lt.pcoa.loadings -s \"\\t\" -p loading > alpha_and_beta_diversity/$analysis.braycurtis.$distance.lt.pcoa.loadings.json;";
    $cmd .= "$bin/visPCoA.pl -a alpha_and_beta_diversity/$analysis.braycurtis.$distance.lt.pcoa.axes.json -l alpha_and_beta_diversity/$analysis.braycurtis.$distance.lt.pcoa.loadings.json -t $template_dir/pcoa.html> html/$analysis.braycurtis.$distance.lt.pcoa.html;";
    runcmds($cpus, $cmd);
    
    $cmd = "$bin/createJson.pl -c alpha_and_beta_diversity/$analysis.jclass.$distance.lt.pcoa.axes -s \"\\t\" -p axes > alpha_and_beta_diversity/$analysis.jclass.$distance.lt.pcoa.axes.json;";
    $cmd .= "$bin/createJson.pl -c alpha_and_beta_diversity/$analysis.jclass.$distance.lt.pcoa.loadings -s \"\\t\" -p loading > alpha_and_beta_diversity/$analysis.jclass.$distance.lt.pcoa.loadings.json;";
    $cmd .= "$bin/visPCoA.pl -a alpha_and_beta_diversity/$analysis.jclass.$distance.lt.pcoa.axes.json -l alpha_and_beta_diversity/$analysis.jclass.$distance.lt.pcoa.loadings.json -t $template_dir/pcoa.html> html/$analysis.jclass.$distance.lt.pcoa.html;";
    runcmds($cpus, $cmd);
    
    $cmd = "$bin/createJson.pl -c alpha_and_beta_diversity/$analysis.jest.$distance.lt.pcoa.axes -s \"\\t\" -p axes > alpha_and_beta_diversity/$analysis.jest.$distance.lt.pcoa.axes.json;";
    $cmd .= "$bin/createJson.pl -c alpha_and_beta_diversity/$analysis.jest.$distance.lt.pcoa.loadings -s \"\\t\" -p loading > alpha_and_beta_diversity/$analysis.jest.$distance.lt.pcoa.loadings.json;";
    $cmd .= "$bin/visPCoA.pl -a alpha_and_beta_diversity/$analysis.jest.$distance.lt.pcoa.axes.json -l alpha_and_beta_diversity/$analysis.jest.$distance.lt.pcoa.loadings.json -t $template_dir/pcoa.html> html/$analysis.jest.$distance.lt.pcoa.html;";
    runcmds($cpus, $cmd);
    
    $cmd = "$bin/createJson.pl -c alpha_and_beta_diversity/$analysis.thetayc.$distance.lt.pcoa.axes -s \"\\t\" -p axes > alpha_and_beta_diversity/$analysis.thetayc.$distance.lt.pcoa.axes.json;";
    $cmd .= "$bin/createJson.pl -c alpha_and_beta_diversity/$analysis.thetayc.$distance.lt.pcoa.loadings -s \"\\t\" -p loading > alpha_and_beta_diversity/$analysis.thetayc.$distance.lt.pcoa.loadings.json;";
    $cmd .= "$bin/visPCoA.pl -a alpha_and_beta_diversity/$analysis.thetayc.$distance.lt.pcoa.axes.json -l alpha_and_beta_diversity/$analysis.thetayc.$distance.lt.pcoa.loadings.json -t $template_dir/pcoa.html> html/$analysis.thetayc.$distance.lt.pcoa.html;";
    runcmds($cpus, $cmd);
    
    $cmd = "$bin/createJson.pl -c alpha_and_beta_diversity/$analysis.subsample.braycurtis.$distance.lt.pcoa.axes -s \"\\t\" -p axes > alpha_and_beta_diversity/$analysis.subsample.braycurtis.$distance.lt.pcoa.axes.json;";
    $cmd .= "$bin/createJson.pl -c alpha_and_beta_diversity/$analysis.subsample.braycurtis.$distance.lt.pcoa.loadings -s \"\\t\" -p loading > alpha_and_beta_diversity/$analysis.subsample.braycurtis.$distance.lt.pcoa.loadings.json;";
    $cmd .= "$bin/visPCoA.pl -a alpha_and_beta_diversity/$analysis.subsample.braycurtis.$distance.lt.pcoa.axes.json -l alpha_and_beta_diversity/$analysis.subsample.braycurtis.$distance.lt.pcoa.loadings.json -t $template_dir/pcoa.html> html/$analysis.subsample.braycurtis.$distance.lt.pcoa.html;";
    runcmds($cpus, $cmd);
    
    $cmd = "$bin/createJson.pl -c alpha_and_beta_diversity/$analysis.subsample.jclass.$distance.lt.pcoa.axes -s \"\\t\" -p axes > alpha_and_beta_diversity/$analysis.subsample.jclass.$distance.lt.pcoa.axes.json;";
    $cmd .= "$bin/createJson.pl -c alpha_and_beta_diversity/$analysis.subsample.jclass.$distance.lt.pcoa.loadings -s \"\\t\" -p loading > alpha_and_beta_diversity/$analysis.subsample.jclass.$distance.lt.pcoa.loadings.json;";
    $cmd .= "$bin/visPCoA.pl -a alpha_and_beta_diversity/$analysis.subsample.jclass.$distance.lt.pcoa.axes.json -l alpha_and_beta_diversity/$analysis.subsample.jclass.$distance.lt.pcoa.loadings.json -t $template_dir/pcoa.html> html/$analysis.subsample.jclass.$distance.lt.pcoa.html;";
    runcmds($cpus, $cmd);
    
    $cmd = "$bin/createJson.pl -c alpha_and_beta_diversity/$analysis.subsample.jest.$distance.lt.pcoa.axes -s \"\\t\" -p axes > alpha_and_beta_diversity/$analysis.subsample.jest.$distance.lt.pcoa.axes.json;";
    $cmd .= "$bin/createJson.pl -c alpha_and_beta_diversity/$analysis.subsample.jest.$distance.lt.pcoa.loadings -s \"\\t\" -p loading > alpha_and_beta_diversity/$analysis.subsample.jest.$distance.lt.pcoa.loadings.json;";
    $cmd .= "$bin/visPCoA.pl -a alpha_and_beta_diversity/$analysis.subsample.jest.$distance.lt.pcoa.axes.json -l alpha_and_beta_diversity/$analysis.subsample.jest.$distance.lt.pcoa.loadings.json -t $template_dir/pcoa.html> html/$analysis.subsample.jest.$distance.lt.pcoa.html;";
    runcmds($cpus, $cmd);
    
    $cmd = "$bin/createJson.pl -c alpha_and_beta_diversity/$analysis.subsample.thetayc.$distance.lt.pcoa.axes -s \"\\t\" -p axes > alpha_and_beta_diversity/$analysis.subsample.thetayc.$distance.lt.pcoa.axes.json;";
    $cmd .= "$bin/createJson.pl -c alpha_and_beta_diversity/$analysis.subsample.thetayc.$distance.lt.pcoa.loadings -s \"\\t\" -p loading > alpha_and_beta_diversity/$analysis.subsample.thetayc.$distance.lt.pcoa.loadings.json;";
    $cmd .= "$bin/visPCoA.pl -a alpha_and_beta_diversity/$analysis.subsample.thetayc.$distance.lt.pcoa.axes.json -l alpha_and_beta_diversity/$analysis.subsample.thetayc.$distance.lt.pcoa.loadings.json -t $template_dir/pcoa.html> html/$analysis.subsample.thetayc.$distance.lt.pcoa.html;";
    runcmds($cpus, $cmd);
    

    #NMDS
    $cmd = "$bin/createJson.pl -c alpha_and_beta_diversity/$analysis.thetayc.$distance.lt.nmds.axes -s \"\\t\" -p axes > alpha_and_beta_diversity/$analysis.thetayc.$distance.lt.nmds.axes.json;";
    $cmd .= "$bin/createJson.pl -c alpha_and_beta_diversity/$analysis.thetayc.$distance.lt.nmds.stress -s \"\\t\" -p stress > alpha_and_beta_diversity/$analysis.thetayc.$distance.lt.nmds.stress.json;";
    $cmd .= "$bin/visNMDs.pl -s alpha_and_beta_diversity/$analysis.thetayc.$distance.lt.nmds.stress.json -a alpha_and_beta_diversity/$analysis.thetayc.$distance.lt.nmds.axes.json -t $template_dir/nmds.html> html/$analysis.thetayc.$distance.lt.nmds.html;";
    runcmds($cpus, $cmd);

    $cmd = "$bin/createJson.pl -c alpha_and_beta_diversity/$analysis.jclass.$distance.lt.nmds.axes -s \"\\t\" -p axes > alpha_and_beta_diversity/$analysis.jclass.$distance.lt.nmds.axes.json;";
    $cmd .= "$bin/createJson.pl -c alpha_and_beta_diversity/$analysis.jclass.$distance.lt.nmds.stress -s \"\\t\" -p stress > alpha_and_beta_diversity/$analysis.jclass.$distance.lt.nmds.stress.json;";
    $cmd .= "$bin/visNMDs.pl -s alpha_and_beta_diversity/$analysis.jclass.$distance.lt.nmds.stress.json -a alpha_and_beta_diversity/$analysis.jclass.$distance.lt.nmds.axes.json -t $template_dir/nmds.html> html/$analysis.jclass.$distance.lt.nmds.html;";
    runcmds($cpus, $cmd);


    $cmd = "$bin/createJson.pl -c alpha_and_beta_diversity/$analysis.jest.$distance.lt.nmds.axes -s \"\\t\" -p axes > alpha_and_beta_diversity/$analysis.jest.$distance.lt.nmds.axes.json;";
    $cmd .= "$bin/createJson.pl -c alpha_and_beta_diversity/$analysis.jest.$distance.lt.nmds.stress -s \"\\t\" -p stress > alpha_and_beta_diversity/$analysis.jest.$distance.lt.nmds.stress.json;";
    $cmd .= "$bin/visNMDs.pl -s alpha_and_beta_diversity/$analysis.jest.$distance.lt.nmds.stress.json -a alpha_and_beta_diversity/$analysis.jest.$distance.lt.nmds.axes.json -t $template_dir/nmds.html> html/$analysis.jest.$distance.lt.nmds.html;";
    runcmds($cpus, $cmd);

    $cmd = "$bin/createJson.pl -c alpha_and_beta_diversity/$analysis.braycurtis.$distance.lt.nmds.axes -s \"\\t\" -p axes > alpha_and_beta_diversity/$analysis.braycurtis.$distance.lt.nmds.axes.json;";
    $cmd .= "$bin/createJson.pl -c alpha_and_beta_diversity/$analysis.braycurtis.$distance.lt.nmds.stress -s \"\\t\" -p stress > alpha_and_beta_diversity/$analysis.braycurtis.$distance.lt.nmds.stress.json;";
    $cmd .= "$bin/visNMDs.pl -s alpha_and_beta_diversity/$analysis.braycurtis.$distance.lt.nmds.stress.json -a alpha_and_beta_diversity/$analysis.braycurtis.$distance.lt.nmds.axes.json -t $template_dir/nmds.html> html/$analysis.braycurtis.$distance.lt.nmds.html;";
    runcmds($cpus, $cmd);
    
    
    $cmd = "$bin/createJson.pl -c alpha_and_beta_diversity/$analysis.subsample.thetayc.$distance.lt.nmds.axes -s \"\\t\" -p axes > alpha_and_beta_diversity/$analysis.subsample.thetayc.$distance.lt.nmds.axes.json;";
    $cmd .= "$bin/createJson.pl -c alpha_and_beta_diversity/$analysis.subsample.thetayc.$distance.lt.nmds.stress -s \"\\t\" -p stress > alpha_and_beta_diversity/$analysis.subsample.thetayc.$distance.lt.nmds.stress.json;";
    $cmd .= "$bin/visNMDs.pl -s alpha_and_beta_diversity/$analysis.subsample.thetayc.$distance.lt.nmds.stress.json -a alpha_and_beta_diversity/$analysis.subsample.thetayc.$distance.lt.nmds.axes.json -t $template_dir/nmds.html> html/$analysis.subsample.thetayc.$distance.lt.nmds.html;";
    runcmds($cpus, $cmd);

    $cmd = "$bin/createJson.pl -c alpha_and_beta_diversity/$analysis.subsample.jclass.$distance.lt.nmds.axes -s \"\\t\" -p axes > alpha_and_beta_diversity/$analysis.subsample.jclass.$distance.lt.nmds.axes.json;";
    $cmd .= "$bin/createJson.pl -c alpha_and_beta_diversity/$analysis.subsample.jclass.$distance.lt.nmds.stress -s \"\\t\" -p stress > alpha_and_beta_diversity/$analysis.subsample.jclass.$distance.lt.nmds.stress.json;";
    $cmd .= "$bin/visNMDs.pl -s alpha_and_beta_diversity/$analysis.subsample.jclass.$distance.lt.nmds.stress.json -a alpha_and_beta_diversity/$analysis.subsample.jclass.$distance.lt.nmds.axes.json -t $template_dir/nmds.html> html/$analysis.subsample.jclass.$distance.lt.nmds.html;";
    runcmds($cpus, $cmd);


    $cmd = "$bin/createJson.pl -c alpha_and_beta_diversity/$analysis.subsample.jest.$distance.lt.nmds.axes -s \"\\t\" -p axes > alpha_and_beta_diversity/$analysis.subsample.jest.$distance.lt.nmds.axes.json;";
    $cmd .= "$bin/createJson.pl -c alpha_and_beta_diversity/$analysis.subsample.jest.$distance.lt.nmds.stress -s \"\\t\" -p stress > alpha_and_beta_diversity/$analysis.subsample.jest.$distance.lt.nmds.stress.json;";
    $cmd .= "$bin/visNMDs.pl -s alpha_and_beta_diversity/$analysis.subsample.jest.$distance.lt.nmds.stress.json -a alpha_and_beta_diversity/$analysis.subsample.jest.$distance.lt.nmds.axes.json -t $template_dir/nmds.html> html/$analysis.subsample.jest.$distance.lt.nmds.html;";
    runcmds($cpus, $cmd);

    $cmd = "$bin/createJson.pl -c alpha_and_beta_diversity/$analysis.subsample.braycurtis.$distance.lt.nmds.axes -s \"\\t\" -p axes > alpha_and_beta_diversity/$analysis.subsample.braycurtis.$distance.lt.nmds.axes.json;";
    $cmd .= "$bin/createJson.pl -c alpha_and_beta_diversity/$analysis.subsample.braycurtis.$distance.lt.nmds.stress -s \"\\t\" -p stress > alpha_and_beta_diversity/$analysis.subsample.braycurtis.$distance.lt.nmds.stress.json;";
    $cmd .= "$bin/visNMDs.pl -s alpha_and_beta_diversity/$analysis.subsample.braycurtis.$distance.lt.nmds.stress.json -a alpha_and_beta_diversity/$analysis.subsample.braycurtis.$distance.lt.nmds.axes.json -t $template_dir/nmds.html> html/$analysis.subsample.braycurtis.$distance.lt.nmds.html;";
    runcmds($cpus, $cmd);
    
    

    $cmd = "$bin/visNewickTree.pl -i alpha_and_beta_diversity/$analysis.braycurtis.$distance.lt.tre > html/$analysis.braycurtis.$distance.lt.tre.html;";
    $cmd .= "$bin/visNewickTree.pl -i alpha_and_beta_diversity/$analysis.thetayc.$distance.lt.tre > html/$analysis.thetayc.$distance.lt.tre.html;";
    $cmd .= "$bin/visNewickTree.pl -i alpha_and_beta_diversity/$analysis.jclass.$distance.lt.tre > html/$analysis.jclass.$distance.lt.tre.html;";
    $cmd .= "$bin/visNewickTree.pl -i alpha_and_beta_diversity/$analysis.jest.$distance.lt.tre > html/$analysis.jest.$distance.lt.tre.html;";

    $cmd .= "$bin/visNewickTree.pl -i alpha_and_beta_diversity/$analysis.subsample.braycurtis.$distance.lt.tre > html/$analysis.subsample.braycurtis.$distance.lt.tre.html;";
    $cmd .= "$bin/visNewickTree.pl -i alpha_and_beta_diversity/$analysis.subsample.thetayc.$distance.lt.tre > html/$analysis.subsample.thetayc.$distance.lt.tre.html;";
    $cmd .= "$bin/visNewickTree.pl -i alpha_and_beta_diversity/$analysis.subsample.jclass.$distance.lt.tre > html/$analysis.subsample.jclass.$distance.lt.tre.html;";
    $cmd .= "$bin/visNewickTree.pl -i alpha_and_beta_diversity/$analysis.subsample.jest.$distance.lt.tre > html/$analysis.subsample.jest.$distance.lt.tre.html;";

    runcmds($cpus, $cmd);


}

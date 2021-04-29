#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Long;
use Statistics::Descriptive;
use Number::Format qw(:subs);
use FindBin;
use Benchmark;
use File::Copy qw(move);
use Data::Dumper qw(Dumper);
use File::Copy;
use Time::Piece;
use Scalar::Util qw(openhandle);
use FindBin;
use lib "$FindBin::Bin";
require "util.pl";
require "global.pl";


our ($bin, $dbdir, $usearch, $mothur, $template, $taxonomy, $fprimers, $rprimers, $proj, $cpus);
our %mapHash;
our %fname2sname;

my ($mappingFile,$analysisname,$seqtype, $oligos, $seqformat, $pdiffs, $email, $maxN, $trunclenF,$trunclenR,$maxeeF, $maxeeR, $truncqF, $truncqR, $minoverlen, $maxdiffs);

$seqtype = "paired";
$trunclenF = 250;
$trunclenR = 220;
$maxeeF = 1;
$maxeeR = 1;
$maxN = 0;
$truncqF = 2;
$truncqR = 2;
$minoverlen = 12;
$maxdiffs = 0;
$seqformat = "fastq";
$pdiffs = 0;
$email = "xdong\@ucalgary.ca";
my $quiet = 0;


&GetOptions("map=s" =>\$mappingFile,
	    "an=s" =>\$analysisname,
	    "seqtype=s" =>\$seqtype,
	    "oligos=s" =>\$oligos,
	    "pdiffs=i" =>\$pdiffs,
	    "email=s" =>\$email,
	    "maxN=i" =>\$maxN,
	    "maxeeF=s" =>\$maxeeF,
	    "maxeeR=s" =>\$maxeeR,
	    "truncqF=i" =>\$truncqF,
	    "truncqR=i" =>\$truncqR,
	    "trunclenF=i" =>\$trunclenF,
	    "trunclenR=i" =>\$trunclenR,

	    "minoverlen=i" =>\$minoverlen,
	    "maxdiffs=i" =>\$maxdiffs,
	    "q=i" =>\$quiet
    );

($mappingFile && $analysisname && $oligos) or
    die "\nusage: $0 \n".
    "-map <reads mapping file>\n".
    "-an <analysis name>\n".
    "-seqtype <paired|single, default is \"paired\">\n".
    "-oligos <oligo file>\n".
    "-pdiffs <default 0: maximum number of differences to the primer sequence>\n".
    "-email <default xdong\@ucalgary.ca>\n".
    "-maxN <number discard reads with maxn number of Ns default is: 0>\n".
    "-maxeeF <float, discard reads with expected errors greater than this value, default is: 1>\n".
    "-maxeeR <float, discard reads with expected errors greater than this value, default is: 1>\n".
    "-trunclenF <number, trim the single end sequences to the fixed length>\n".
    "-trunclenR <number, trim the single end sequences to the fixed length>\n".
    "-truncqF <default is: 2>\n".
    "-truncqR <default is: 2>\n".
    "-minoverlen <default is:12>\n".
    "-maxdiffs <max number of the mismatch in the overlap region:0>\n".
    "-q <quiet default 0>\n".
    "metaamp analysis pipeline: assume all the input reads are in fastq format\n";

my $distance = "asv";
# Record the command line
$analysisname =~ tr/[ \-]/_/;
my $params = "Parameters for running metaamp: $0 -map $mappingFile -an $analysisname -seqtype $seqtype -oligos $oligos -pdiffs $pdiffs -email $email -maxN $maxN  -maxeeF $maxeeF -maxeeR $maxeeR -trunclenF $trunclenF -trunclenR $trunclenR  -truncqF $truncqF -truncqR $truncqR -minoverlen $minoverlen -maxdiffs $maxdiffs\n";

my $t0 = Benchmark->new;

# Record process id in case we must kill this  job
msg("Command line:\n $params\n");
msg("Process ID = $$\n\n");

unless (-e "JobInfo.txt"){
save_job_info("JobInfo.txt", $analysisname, $email, $params);
}
msg("Reading mapping file\n");
#my $err = readMappingFile($mappingFile, $seqtype, \%mapHash);
my $err = readMappingFile($mappingFile, $seqtype, $seqformat, \%mapHash);
#Error message in the mapping file
if($err ne ""){
    msg("Mapping file errors: $err\n");
    exit(4);
}

my $cmd = "$bin/seqStats.pl -f $seqformat -t file -s raw_file.list.txt > metaamp.raw.stats.csv;";
$cmd .= "$bin/createJson.pl -c  metaamp.raw.stats.csv -p raw > metaamp.raw.stats.json;";
runcmds($cpus, $cmd);

msg("Making design file for the later on sample comparsion and multivariate analysis");
makeDesignFile(\%mapHash);

msg("reading primer informaiton");
$err = getPrimerInfo($oligos);
if($err ne ""){
    msg("Primer validation errors: $err\n");
    exit(4);
}

my $runids = trimPrimers_dada2(\%mapHash, $seqtype, $fprimers, $rprimers, $cpus);
#getPrimerTrimStats(\%mapHash, $seqtype, \%sampleStats);

$cmd = "$bin/seqStats.pl -f $seqformat -t file -s after_strip_primer_file.list.txt > metaamp.noprimer.stats.csv;";
$cmd .= "$bin/createJson.pl -c  metaamp.noprimer.stats.csv -p trim > metaamp.noprimer.stats.json;";
runcmds($cpus, $cmd);

my $cat_str = "cat ";

my $truncl_str = "";
if($trunclenF > 0){

    $truncl_str .= "--truncLen_F $trunclenF "
}
if($trunclenR > 0){

    $truncl_str .= "--truncLen_R $trunclenR "
}
foreach my $runid (keys %$runids){
#####TODO ASV###########################
#RScript

    $cmd = "Rscript $bin/run_dada2_paired.R -n $analysisname -p NOPRIM.$runid  -i $runid $truncl_str --maxEE_F $maxeeF --maxEE_R $maxeeR --maxN $maxN --truncqF $truncqF --truncqR $truncqR;";
    runcmds($cpus, $cmd);
    $cat_str .= "$runid\.dada2_filter_trim_file.list.txt "
}

runcmds($cpus, "$cat_str >  dada2_filter_trim_file.list.txt");

my $runids_str = join(",", keys %$runids);
$cmd = "Rscript $bin/run_dada2_merge_runs.R -r $runids_str -a $analysisname";
runcmds($cpus, $cmd);


my $asvFasta = "$analysisname.asv.fasta";
my $asvShared = "$analysisname.shared";
my @samplenames = keys %mapHash;

$cmd = "$bin/get_diversity.pl -s $asvShared -t $cpus -l $distance";
runcmds($cpus, $cmd);

$cmd = "$bin/seqStats.pl -f $seqformat -t file -s dada2_filter_trim_file.list.txt > metaamp.qc.stats.csv;";
$cmd .= "$bin/createJson.pl -c metaamp.qc.stats.csv -p qc > metaamp.qc.stats.json;";
runcmds($cpus, $cmd);


getTaxonInfo($analysisname,$cpus, $distance, $asvFasta);

$cmd ="$bin/packing.pl -an $analysisname -l $distance";
runcmds($cpus, $cmd);

$cmd = "$bin/output_html_report.pl -an $analysisname -seqtype $seqtype -dir results -d $distance;";
runcmds($cpus, $cmd);
send_job_finished_email();
delfile(glob "mothur\.*");
my $t1 = Benchmark->new;
my $td = timediff($t1, $t0);
msg("metaamp took:" . timestr($td) . " to run\n");



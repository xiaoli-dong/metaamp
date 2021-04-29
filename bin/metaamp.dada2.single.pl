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

my ($mappingFile,$analysisname,$seqtype, $oligos,$trunclen,$maxee, $seqformat, $pdiffs, $email, $maxN, $truncq);
$seqtype = "paired";
$trunclen = 150;
$maxee = 1;
$maxN = 0;
$truncq = 2;
$seqformat = "fastq";

$pdiffs = 0;
$email = "xdong\@ucalgary.ca";
my $quiet = 0;


&GetOptions("map=s" =>\$mappingFile,
	    "an=s" =>\$analysisname,
	    "seqtype=s" =>\$seqtype,
	    "seqformat=s" =>\$seqformat,
	    "maxee=s" =>\$maxee,
	    "truncq=i" =>\$truncq,
	    "trunclen=i" =>\$trunclen,
	    "maxN=i" =>\$maxN,
	    "oligos=s" =>\$oligos,
	    "pdiffs=i" =>\$pdiffs,
	    "email=s" =>\$email,
	    "q=i" =>\$quiet
    );

($mappingFile && $analysisname && $oligos) or
    die "\nusage: $0 \n".
    "-map <reads mapping file>\n".
    "-an <analysis name>\n".
    "-seqtype <paired|single, default is \"paired\">\n".
    "-seqformat <fastq|fasta, default is fastq>\n".
    "-trunclen <number, trim the single end sequences to the fixed length>\n".
    "-maxee <float, discard reads with expected errors greater than this value, default is: 1>\n".
    "-maxN <number discard reads with maxn number of Ns default is: 0>\n".
    "-truncq <number discard reads with maxn number of Ns default is: 2>\n".
    "-oligos <oligo file>\n".
    "-pdiffs <default 0: maximum number of differences to the primer sequence>\n".
    "-email <default xdong\@ucalgary.ca>\n".
    "-q <quiet default 0>\n".
    "metaamp analysis pipeline: assume all the input reads are in fastq format\n";

my $distance = "asv";
# Record the command line
$analysisname =~ tr/[ \-]/_/;
my $params = "$0 -map $mappingFile -an $analysisname -seqtype $seqtype -seqformat $seqformat -trunclen $trunclen -maxee $maxee -oligos $oligos -pdiffs $pdiffs -email $email\n";

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
foreach my $runid (keys %$runids){
#####TODO ASV###########################
#RScript
    $cmd = "Rscript $bin/run_dada2_single.R -n $analysisname -p NOPRIM.$runid -i $runid --truncLen $trunclen --maxEE $maxee  --maxN $maxN --truncq $truncq";
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

$cmd = "$bin/output_html_report.pl -an $analysisname -seqtype $seqtype -dir results -d $distance";
runcmds($cpus, $cmd);
send_job_finished_email();
delfile(glob "mothur\.*");
my $t1 = Benchmark->new;
my $td = timediff($t1, $t0);
msg("metaamp took:" . timestr($td) . " to run\n");



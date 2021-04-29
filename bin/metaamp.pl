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

our ($bin, $dbdir, $usearch, $mothur, $template, $taxonomy, $fprimers, $rprimers, $proj);
our %mapHash;
our %fname2sname;

my ($mappingFile,$analysisname,$seqtype, $scutoff, $oligos, $minoverlen, $maxdiffs,$trunclen,$maxee, $seqformat, $pdiffs, $email);
$seqtype = "paired";
$scutoff = 0.97;
$trunclen = 150;
$maxee = 1;
$seqformat = "fastq";
$minoverlen = 8;
$maxdiffs = 0;
$pdiffs = 0;
$email = "xdong\@ucalgary.ca";
my $quiet = 0;
my $cpus = 8;

&GetOptions("map=s" =>\$mappingFile,
	    "an=s" =>\$analysisname,
	    "seqtype=s" =>\$seqtype,
	    "seqformat=s" =>\$seqformat,
	    "s=f" =>\$scutoff,
	    "minoverlen=i" =>\$minoverlen,
	    "maxdiffs=i" =>\$maxdiffs,
	    "trunclen=i" =>\$trunclen,
	    "maxee=s" =>\$maxee,
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
    "-s <similarity cutoff for generating otus, default is 0.97>\n".
    "-minoverlen <minimum overlap lenght while mergeing paired end reads, default is 8)\n".
    "-maxdiffs <maximum number of mismatches allowed in the overlap regions>\n".
    "-trunclen <number, trim the single end sequences to the fixed length>\n".
    "-maxee <float, discard reads with expected errors greater than this value>\n".
    "-oligos <oligo file>\n".
    "-pdiffs <default 0: maximum number of differences to the primer sequence>\n".
    "-email <default xdong\@ucalgary.ca>\n".
    "-q <quiet default 0>\n".
    "metaamp analysis pipeline: assume all the input reads are in fastq format\n";


# Record the command line
$analysisname =~ tr/[ \-]/_/;
my $params = "$0 -map $mappingFile -an $analysisname -seqtype $seqtype -seqformat $seqformat -s $scutoff -minoverlen $minoverlen -maxdiffs $maxdiffs  -trunclen $trunclen -maxee $maxee -oligos $oligos -pdiffs $pdiffs -email $email\n";

my $t0 = Benchmark->new;

# Record process id in case we must kill this  job
msg("Command line:\n $params\n");
msg("Process ID = $$\n\n");

unless (-e "JobInfo.txt"){
save_job_info("JobInfo.txt", $analysisname, $email, $params);
}

my $distance = 1 - $scutoff;



msg("Reading mapping file\n");
my $err = readMappingFile($mappingFile, $seqtype, $seqformat, \%mapHash);
#Error message in the mapping file
if($err ne ""){
    msg("Mapping file errors: $err\n");
    exit(4);
}

my $cmd = "$bin/seqStats.pl -f $seqformat -t file -s raw_file.list.txt > metaamp.raw.stats.csv;";
runcmds($cpus, $cmd);

msg("Making design file for the later on sample comparsion and multivariate analysis");
makeDesignFile(\%mapHash);

msg("reading primer informaiton");
$err = getPrimerInfo($oligos);
if($err ne ""){
    msg("Primer validation errors: $err\n");
    exit(4);
}

###############start QC######################
#output, single-end: $samplename.NOPRIM.fastq
#paired-end: $samplename_R1.NOPRIM.fastq;";$samplename_R2.NOPRIM.fastq;";

trimPrimers_uparse(\%mapHash, $seqtype, $fprimers, $rprimers);

$cmd = "$bin/seqStats.pl -f $seqformat -t file -s raw_file.list.txt > meteaamp.raw.stats.csv;";
$cmd .= "$bin/createJson.pl -c  metaamp.raw.stats.csv -p raw > metaamp.raw.stats.json;";
$cmd .= "$bin/seqStats.pl -f $seqformat -t file -s after_strip_primer_file.list.txt > metaamp.noprimer.stats.csv;";
$cmd .= "$bin/createJson.pl -c  metaamp.noprimer.stats.csv -p trim > metaamp.noprimer.stats.json;";
runcmds($cpus, $cmd);

#merging pairs
if($seqtype eq "paired"){
    $cmd = "$usearch -fastq_mergepairs *R1.NOPRIM.fastq  -fastq_minlen 50 -fastq_maxdiffs $maxdiffs -fastq_minovlen $minoverlen -fastqout $analysisname.merged.fastq -threads $cpus -log merged.log;";
    $cmd .= "$bin/seqStats.uparse.pl -f fastq -i $analysisname.merged.fastq > metaamp.merged.stats.csv;";
    $cmd .= "$bin/createJson.pl -c  metaamp.merged.stats.csv -p merge > metaamp.merged.stats.json;";
    $cmd .= "$usearch -fastq_filter $analysisname.merged.fastq -threads $cpus -fastq_trunclen $trunclen -fastq_maxee $maxee -fastaout $analysisname.qc.fasta;";
    runcmds($cpus, $cmd);

}
#single-end, need to pool samples first
else{
    #pooling samples
    delfile("$analysisname.NOPRIM.fastq") if -e "$analysisname.NOPRIM.fastq";
    foreach my $samplename (keys %mapHash){
	$cmd = "cat $samplename.NOPRIM.fastq >> $analysisname.NOPRIM.fastq";
	runcmds($cpus, $cmd);
    }

    $cmd = "$usearch -fastq_filter $analysisname.NOPRIM.fastq -threads $cpus -fastq_trunclen $trunclen -fastq_maxee $maxee -fastaout $analysisname.qc.fasta;";
    runcmds($cpus, $cmd);
}

$cmd = "$bin/seqStats.uparse.pl -f fasta -i $analysisname.qc.fasta > metaamp.qc.stats.csv;";
$cmd .= "$bin/createJson.pl -c  metaamp.qc.stats.csv -p qc > metaamp.qc.stats.json;";
runcmds($cpus, $cmd);

###############END QC######################

makeCluster("$analysisname.qc.fasta", $cpus, $analysisname, $scutoff);
$cmd = "$bin/get_diversity.pl -s $analysisname.shared -t $cpus -l $distance";
runcmds($cpus, $cmd);
getTaxonInfo($analysisname,$cpus, $distance,  "$analysisname.otu.fasta");

$cmd ="$bin/packing.pl -an $analysisname -l $distance;";
runcmds($cpus, $cmd);

$cmd = "$bin/output_html_report.pl -an $analysisname -seqtype $seqtype -dir results -d $distance;";
runcmds($cpus, $cmd);
send_job_finished_email();
delfile(glob "mothur\.*");
my $t1 = Benchmark->new;
my $td = timediff($t1, $t0);
msg("metaamp took:" . timestr($td) . " to run\n");


sub makeCluster{

    my($fasta, $cpus, $analysisname, $scutoff) = @_;
# dereplicate seqs
#my $cmd = "$usearch -fastx_uniques $fasta -threads 8 -fastaout $analysisname.derep.fasta -sizeout;";
    my $cmd = "$usearch -fastx_uniques $fasta -threads $cpus -fastaout $analysisname.derep.fasta -sizeout;";
    runcmds($cpus, $cmd);

# cluster OTUs and discarding singletons
    $cmd = "$usearch -cluster_otus $analysisname.derep.fasta  -threads $cpus -otus $analysisname.rep.fasta -minsize 2 -relabel OTU_;";
    runcmds($cpus, $cmd);

#$cmd = "$usearch -usearch_global $fasta -threads 8 -db $analysisname.otus.fasta -strand plus -id 0.97 -uc otu.map.uc -otutabout otu_table.txt -biomout otu_table.json";
    $cmd = "$usearch -usearch_global $fasta -threads $cpus -db $analysisname.rep.fasta -strand plus -id $scutoff -uc $analysisname.otu.map.uc -otutabout $analysisname.otu_table.txt -biomout $analysisname.otu_table.json";
    runcmds($cpus, $cmd);

    #make shared file
    open TABLE, "$analysisname.otu_table.txt" or die "Could not open $analysisname.otu_table.txt to read, $!\n";
        my %table = ();
    my $total_otu_count = 0;
    my @samples = ();
    my @otuids = ();
    my %otu2count = ();
    while(<TABLE>){
	chomp;
	next if !/\S/;

	if(/^#/){
	    @samples = split(/\t/, $_);
	    shift @samples;
	    next;
	}
	my @l = split(/\t/, $_);
	my $otuid = shift @l;
	$otuid =~ s/OTU_//;
	for (my $i=0; $i < @samples; $i++) {
	    $table{$samples[$i]}->{$otuid} = $l[$i];
	    $otu2count{$otuid} += $l[$i];
	}
	push(@otuids, $otuid);
	$total_otu_count++;
    }
    close(TABLE);
    my @sorted_otuids = sort {$a <=> $b} @otuids;
    open SHARED, ">$analysisname.shared" or die "Could not open $analysisname.shared to write, $!\n";
    my @header = ("label", "Group", "numOtus");
    push(@header, map{"OTU_" . $_} @sorted_otuids);
    print SHARED join("\t", @header), "\n";
    foreach my $sample (@samples){
	my @row = (1-$scutoff, $sample, $total_otu_count);
	foreach my $otuid (@sorted_otuids){
	    push(@row, $table{$sample}->{$otuid});
	}
	print SHARED join("\t", @row), "\n";
    }
    close(SHARED);

    #make OTU fasta file
    open REP, "$analysisname.rep.fasta" or die "Could not open $analysisname.rep.fasta to read, $!\n";
    open OTUS, ">$analysisname.otu.fasta" or die "Could not open $analysisname.otu.fasta to write, $!\n";
    while(<REP>){

	if(/^>(\S+)/){
	    my $id = $1;
	    $id =~ s/OTU_//;
	    if(exists $otu2count{$id}){
		print OTUS ">OTU\_$id\|", $otu2count{$id}, "\n";
	    }
	    else{
		msg("OTU\_$id has not count");
	    }

	}
	else{
	    print OTUS $_;
	}
    }
    close(REP);
    close(OTUS);
}



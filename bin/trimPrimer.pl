#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Long;
use File::Copy qw(move);
use FindBin;
use FindBin;
use lib "$FindBin::Bin";
require "util.pl";

my ($seqtype, $read1, $read2, $samplename,$fprimer_strs, $rprimer_strs, $cpus);
$seqtype = "paired";
$cpus = 1;

&GetOptions(
    "seqtype=s" =>\$seqtype,
    "r1=s"  =>\$read1,
    "r2=s"  =>\$read2,
    "s=s" =>\$samplename,
    "fp=s" =>\$fprimer_strs,
    "rp=s" =>\$rprimer_strs,
    "j=s" =>\$cpus
    );

($seqtype && $samplename && $read1) or 
    die "\nusage: $0 \n".
    "-seqtype <single|paired, default: paired>\n".
    "-s <sample name>\n".
    "-r1 <read1 file name>\n".
    "-r2 <read2 file name>\n".
    "-fp <default \"\": forward primers in the format of: primer1:primer2...>\n".
    "-rp <default \"\": reverse primers in the format of: primer1:primer2...>\n".
    "-j <default 1: number of cpu cores to use>\n".
    "trimming off amplicon data primers\n";

#perl ../../bin/trimPrimer.pl -seqtype paired -s test -r1 rep1_C1_S244_L001_R1_001.fastq -r2 rep1_C1_S244_L001_R2_001.fastq -fp CCTACGGGAGGCAGCAG -rp GACTACHVGGGTATCTAATCC
my $bin = "$FindBin::RealBin";

$fprimer_strs ||= "";
$rprimer_strs ||= "";
$read2 ||= "";


if($seqtype eq "paired"){
    if($read2 eq "" || ! -e $read2){
	print STDERR "read2 file $read2 not exists";
	exit;
    }
}

my @f_primers = ();
my @r_primers = ();

push(@f_primers, split(":", $fprimer_strs)) if $fprimer_strs ne "";
push(@r_primers, split(":", $rprimer_strs)) if $rprimer_strs ne "";


my $cmd = "cutadapt -e 0 --quiet -j $cpus ";

if($seqtype eq "paired"){
    
    my $r1_linked_adapt = "";
    my $r2_linked_adapt = "";
   
    
    
    if(@f_primers > 0 && @r_primers > 0){
	foreach my $f (@f_primers){
	    my $f_rc = reverse_complement_IUPAC($f);
	    foreach my $r (@r_primers){
		my $r_rc = reverse_complement_IUPAC($r);
		$r1_linked_adapt .= " -a \^$f...$r_rc";
		$r2_linked_adapt .= " -A \^$r...$f_rc";
	    }
	}
	
	$cmd .= $r1_linked_adapt;
	$cmd .= $r2_linked_adapt;
	$cmd .= " --discard-untrimmed";
    }
    elsif(@f_primers > 0){
	foreach my $f (@f_primers){
	    $cmd .= " -g \^$f";
	}
	$cmd .= " --discard-untrimmed";
    }
    elsif(@r_primers > 0){
	foreach my $r (@r_primers){
	    $cmd .= " -G  \^$r";
	}
	$cmd .= " --discard-untrimmed";
    }
    $cmd .= " -o $samplename.R1.trim.fastq -p $samplename.R2.trim.fastq $read1 $read2";
    
}
elsif($seqtype eq "single"){
    
    my $r1_linked_adapt = "";
    if(@f_primers > 0 && @r_primers > 0){
	foreach my $f (@f_primers){
	    my $f_rc = reverse_complement_IUPAC($f);
	    foreach my $r (@r_primers){
		my $r_rc = reverse_complement_IUPAC($r);
		$r1_linked_adapt .= " -a \^$f...$r_rc";
	    }
	}
	$cmd .= $r1_linked_adapt;
	$cmd .= " --discard-untrimmed";
    }
    elsif(@f_primers > 0){
	foreach my $f (@f_primers){
	    $cmd .= " -g \^$f";
	}
	$cmd .= " --discard-untrimmed";
    }
    
    $cmd .= " -o $samplename.trim.fastq $read1";
}
runcmds(1, $cmd);




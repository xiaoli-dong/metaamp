#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Long;
use File::Copy qw(move);
use FindBin;
my ($fastq, $samplename, $oligos, $pdiffs, $cpus);
$pdiffs = 0;
$cpus = 1;
&GetOptions(
    "fastq=s" =>\$fastq,
    "s=s" =>\$samplename,
    "oligos=s" =>\$oligos,
    "pdiffs=i" =>\$pdiffs,
    "t=i" =>\$cpus
    );

($fastq && $samplename) or 
    die "\nusage: $0 \n".
    "-fastq <fastq file>\n".
    "-oligos <oligo file>\n".
    "-s <sample name>\n".
    "-pdiffs <default 0: maximum number of differences to the primer sequence>\n".
    "-t <number of threads>\n".
    "strip off no-biological sequences (primer) from reads and produce fastq file\n";

my $bin = "$FindBin::RealBin";
my $mothur = "$bin/../programs/mothur/mothur";
my ($prefix) = $fastq =~ /(.*)\.[^.]+$/;

#my $cmd = "$^X $bin/fastq2fasta $fastq $prefix 0 $samplename";
#print STDERR "$cmd\n";
#system($cmd) >> 8 and  die "Could not execute cmd=$cmd, $!\n";
#mothur has bug for trim.seqs when using multiple processors
#$cmd = "$mothur \"#trim.seqs(fasta=$prefix.fasta, qfile=$prefix.qual, oligos=oligos.txt, pdiffs=$pdiffs, processors=1);make.fastq();\"";
my $cmd = "$mothur \"#fastq.info(fastq=$fastq,fasta=T, qfile=T);trim.seqs(fasta=current, qfile=current, oligos=oligos.txt, pdiffs=$pdiffs, processors=$cpus);make.fastq();\"";
print STDERR "$cmd\n";
system($cmd) >> 8 and  die "Could not execute cmd=$cmd, $!\n";
#move "$prefix.trim.fastq", "$samplename.strip.fastq";






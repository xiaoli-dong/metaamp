#!/usr/bin/env perl
use strict;	
use warnings;
use List::Util qw(sum);
use Statistics::Descriptive;
use Getopt::Long;

my ($input, $format);
$format = "fasta";

&GetOptions(
    "f=s" =>\$format,
    "i=s" => \$input
    );

($input) or
    die "Usage: $0  OPTIONS\n".
    " -f <input file format, fasta|fastq, default: fasta, it can be compressed gz file>\n".
    " -i <list of sequence file name or patterns seperated by:[pattern1:pattern2..]\n".
    "     or a file with each sequence file name or pattern in a separate line> \n";



my %statInfo = ();


my $f2lens = [];
if($format eq "fasta"){
    $f2lens = readFasta("$input");
}
elsif($format eq "fastq"){
    $f2lens = readFastq("$input");
}

foreach my $sampleid (keys %$f2lens){
    $Statistics::Descriptive::Tolerance = 1e-24;
my $stat = Statistics::Descriptive::Full->new();
    my @lens = @{$f2lens->{$sampleid}};
    print STDERR "$sampleid=", scalar(@lens), "\n";
    if(@lens){
	$stat->add_data(@lens);
	$statInfo{$sampleid}->{min} = $stat->min();
	$statInfo{$sampleid}->{max} = $stat->max();
	$statInfo{$sampleid}->{mean} = $stat->mean();
	$statInfo{$sampleid}->{median} = $stat->median();
	$statInfo{$sampleid}->{mode} = $stat->mode();
	$statInfo{$sampleid}->{count} = $stat->count();
	$statInfo{$sampleid}->{total_bp} = $stat->sum();
	$statInfo{$sampleid}->{std} = $stat->standard_deviation();
	$statInfo{$sampleid}->{n50} = get_N50(\@lens);
    }
}


my @header = ("count", "total_bp", "min", "max", "mean", "median", "mode", "n50", "std");
print "#filename,", join(",", @header), "\n";
foreach my $f (sort keys %statInfo){
    my @values = ();
    push(@values, $f);
    foreach my $col (@header){
	if($col =~ /mean|median|std/){
	    push(@values, sprintf("%0.2f", $statInfo{$f}->{$col}));
	}
	else{
	    push(@values, $statInfo{$f}->{$col});
	}
    }
    print join(",", @values),"\n";
}

sub readFasta{
    my($fastaFile) = @_;
    my %lens = ();
    if ($fastaFile =~ /.gz$/) {
	open(FASTA, "gunzip -c $fastaFile |") or die "Could not open $fastaFile file to read, $!\n";
	$/ = "\n>";
	while (<FASTA>){
	    chomp;
	    if(! />?(\S+?);barcodelabel=(\S+?)\n(.+)/s){
		die "Could not read fasta record #$.: $_\n";
	    }
	    my $id = $1;
	    my $sampleid = $2;
	    my $sequence = $3;
	    $sequence =~ s/\s+//g;
	    my $len = length ($sequence);
	    push(@{$lens{$sampleid}}, $len);
	}
	$/ = "\n";
	close(FASTA);
    }
    else{
	open(FASTA, $fastaFile) or die "Could not open $fastaFile file to read, $!\n";
	$/ = "\n>";
	while (<FASTA>){
	    chomp;
	    if(! />?(\S+?);barcodelabel=(\S+?)\n(.+)/s){
		die "Could not read fasta record #$.: $_\n";
	    }
	    my $sid = $1;
	    my $sampleid = $2;
	    my $sequence = $3;
	    $sequence =~ s/\s+//g;
	    my $len = length ($sequence);
	    push(@{$lens{$sampleid}}, $len);
	}
	$/ = "\n";
	close(FASTA);
    }
    return \%lens;
    
}
sub readFastq{
    my($fastqFile) = @_;
    my %lens = ();
    
    
    if ($fastqFile =~ /.gz$/) {
	
	open(FASTQ, "gunzip -c $fastqFile |") or die "Could not open $fastqFile file to read, $!\n";
	
	while (<FASTQ>){
	    
	    chomp;
	    my($sampleid) = $_ =~ /^\S+?;barcodelabel=(\S+)/;
	    my $sequence = <FASTQ>;
	    <FASTQ>;
	    <FASTQ>;
	    $sequence =~ s/\s+//g;
	    my $len = length ($sequence);
	    push(@{$lens{$sampleid}}, $len);
	}
	close(FASTQ);
	
    }
    else {
	open(FASTQ, $fastqFile) or die "Could not open $fastqFile file to read, $!\n";
	
	while (<FASTQ>){
	    chomp;
	    my($sampleid) = $_ =~ /^\S+?;barcodelabel=(\S+)/;
	    my $sequence = <FASTQ>;
	    <FASTQ>;
	    <FASTQ>;
	    $sequence =~ s/\s+//g;
	    my $len = length ($sequence);
	    push(@{$lens{$sampleid}}, $len);
	}
	close(FASTQ);
    }

    
    return \%lens;
}
sub get_N50{
    my ($arrRef) = @_;
    my @sort = sort {$b <=> $a} @$arrRef;
    my $totalLength = sum(@sort);
    my $n50 = 0;
    my $n50_value = 0;
    foreach my $val(@sort){
	$n50+=$val;
	if($n50 >= $totalLength/2){
	    #print "N50 length is $n50 and N50 value is: $val\n";
	    $n50_value = $val;
	    last;
	}
    }
    return $n50_value;
}


#!/usr/bin/env perl
use strict;	
use warnings;
use List::Util qw(sum);
use Statistics::Descriptive;
use Getopt::Long;

my ($input, $format, $type, $path);
$format = "fasta";
$type = "name";
$path = ".";
&GetOptions(
    "f=s" =>\$format,
    "t=s" =>\$type,
    "s=s" => \$input,
    "p=s" =>\$path
    );

($input) or
    die "Usage: $0  OPTIONS\n".
    " -f <input file format, fasta|fastq, default: fasta, it can be compressed gz file>\n".
    " -t <input type:name|file, default: name>\n" .
    " -p <path of the input file directory, default is \".\">\n" .
    " -s <list of sequence file name or patterns seperated by:[pattern1:pattern2..]\n".
    "     or a file with each sequence file name or pattern in a separate line> \n";


my %files = ();

if($type eq "name"){
    %files = map { $_ => 1 } split(/:/, $input);
    
}
elsif($type eq "file"){
    
    open(INPUT, $input) or die "Could not open $input for read, $!\n";
    while(<INPUT>){
	chomp;
	next if /^#/;
	next if !/\S/;
	$files{$_} = 1;
    }
}

if(! -d $path){
    print STDERR "Input file directory does not exist\n";
    exit;
    
}

my %statInfo = ();
foreach my $file (sort keys %files){
    
    $Statistics::Descriptive::Tolerance = 1e-24;
    my $stat = Statistics::Descriptive::Full->new();
    my $flens = [];
    if($format eq "fasta"){
	$flens = readFasta("$path/$file");
    }
    elsif($format eq "fastq"){
	$flens = readFastq("$path/$file");
    }
    my $c = scalar @$flens;
    
    if(@$flens){
	$stat->add_data(@$flens);
	$statInfo{$file}->{min} = $stat->min();
	$statInfo{$file}->{max} = $stat->max();
	$statInfo{$file}->{mean} = $stat->mean();
	$statInfo{$file}->{median} = $stat->median();
	$statInfo{$file}->{mode} = $stat->mode();
	$statInfo{$file}->{count} = $stat->count();
	$statInfo{$file}->{total_bp} = $stat->sum();
	$statInfo{$file}->{std} = $stat->standard_deviation();
	$statInfo{$file}->{n50} = get_N50($flens);
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
    my @lens = ();
    if ($fastaFile =~ /.gz$/) {
	open(FASTA, "gunzip -c $fastaFile |") or die "Could not open $fastaFile file to read, $!\n";
	$/ = "\n>";
	while (<FASTA>){
	    chomp;
	    if(! />?(\S+).*?\n(.+)/s){
		die "Could not read fasta record #$.: $_\n";
	    }
	    my $id = $1;
	    my $sequence = $2;
	    $sequence =~ s/\s+//g;
	    my $len = length ($sequence);
	    push(@lens, $len);
	}
	$/ = "\n";
	close(FASTA);
    }
    else{
	open(FASTA, $fastaFile) or die "Could not open $fastaFile file to read, $!\n";
	$/ = "\n>";
	while (<FASTA>){
	    chomp;
	    if(! />?(\S+).*?\n(.+)/s){
		die "Could not read fasta record #$.: $_\n";
	    }
	    my $id = $1;
	    my $sequence = $2;
	    $sequence =~ s/\s+//g;
	    my $len = length ($sequence);
	    push(@lens, $len);
	}
	$/ = "\n";
	close(FASTA);
    }
    return \@lens;
    
}
sub readFastq{
    my($fastqFile) = @_;
    my @lens = ();
    
    
    if ($fastqFile =~ /.gz$/) {
	
	open(FASTQ, "gunzip -c $fastqFile |") or die "Could not open $fastqFile file to read, $!\n";
	
	while (<FASTQ>){
	    chomp;
	    my $sequence = <FASTQ>;
	    <FASTQ>;
	    <FASTQ>;
	    $sequence =~ s/\s+//g;
	    my $len = length ($sequence);
	    push(@lens, $len);
	}
	close(FASTQ);
	
    }
    else {
	open(FASTQ, $fastqFile) or die "Could not open $fastqFile file to read, $!\n";
	
	while (<FASTQ>){
	    chomp;
	    my $sequence = <FASTQ>;
	    <FASTQ>;
	    <FASTQ>;
	    $sequence =~ s/\s+//g;
	    my $len = length ($sequence);
	    push(@lens, $len);
	}
	close(FASTQ);
    }

    
    return \@lens;
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


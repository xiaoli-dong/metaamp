#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Long;

my($fastaFile, $qualFile, $format);

$format = "Phred+33";
$qualFile = "";
&GetOptions("fasta=s" =>\$fastaFile,
	    "qual=s" =>\$qualFile,
	    "format=s" =>\$format
    );

($fastaFile) or
    die "\nusage: $0 \n".
    "-fasta <fasta format input file>\n".
    "-qual <default: \"\", fasta format quality file name>\n".
    "-format <default: Phred+33, indicate whether your sequences are Phred+33 or Phred+64>\n".
    "this program convert the input fasta file to phred+33 or phred+64 format fastq file\n";

if(length($qualFile)){

    fasta2fastq($fastaFile, $qualFile, $format);
}
else{

    fasta2fastq_noqual($fastaFile, $format);
}


sub fasta2fastq_noqual{
    my($fastaFile, $format) = @_;
    
    $/ = "\n>";
    open(FASTA, $fastaFile) or die "Could not open $fastaFile to read, $!\n";
    
    while(<FASTA>){
	chomp;
	if(!/>?(\S.*?)\n(.+)/s){
	    die "Could not read FastA record #$.: $_\n";
	}
	my $head = $1;
	my $seq = $2;
	$seq =~ s/[ \r\n\t]//g;
	my $len = length($seq);
	print "\@$head\n$seq\n+\n";
	
	for (my $i=0; $i < $len ; $i++) {
	    print "I";
	}
	print  "\n";
	
	}
    
    $/ = "\n";
    close(FASTA);
}
sub fasta2fastq{

    my($fastaFile, $qualFile, $format) = @_;
    
    open F, "<$fastaFile" or die "Could not open $fastaFile to read, $!\n";
    open Q, "<$qualFile" or die "Could not open $qualFile to read, $!\n";

    my ($fh, $fs, $qh, $qs);
    my @quals;  my $qv;

    my $qline = "";
    $/ = "\n>";

    while (<F>) {
	chomp;
	if(($fh,$fs) =  /^>?(\S+).*?\n(.*)/s){
	    #print "head=$fh, seq=$fs\n";
	    $fs =~ s/\s+//g;
	}
	$qline = <Q>;
	chomp($qline);
	if(($qh,$qs) =  $qline =~ /^>?(\S+).*?\n(.*)/s){
	    #print "qhead=$qh, seq=$qs\n";
	    $qs =~ s/\s+$//g;
	}
	
	if (($fh eq $qh) && ($fh =~ /^(\S+)/)) {
	    
	    print "\@$fh\n$fs\n";
	    print "+$qh\n";
	    
	    @quals = split(/\s+/,$qs);
	    for (my $i=0; $i<=$#quals; $i++) {
		if($format eq "Phred+33"){
		    $qv = $quals[$i] + 33;
		}
		elsif($format eq "Phred+64"){
		    $qv = $quals[$i] + 64;
		}
		else{
		    print STDERR "unknow fastq format $format\n";
		    exit;
		}
		print chr($qv);
	    }
	    print "\n";
	}
	else { print "sequence id is not matching the id in the quality file\n" }
    }
    $/ = "\n";

    close F; close Q;
}

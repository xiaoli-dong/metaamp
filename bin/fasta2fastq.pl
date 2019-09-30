#!/usr/bin/perl

# AUTHOR: Joseph Fass
# LAST REVISED: August 2009
# 
# The Bioinformatics Core at UC Davis Genome Center
# http://bioinformatics.ucdavis.edu
# Copyright (c) 2008 The Regents of University of California, Davis Campus.
# All rights reserved.

use strict;

my $usage = "\nusage: $0 <reads (fasta format)> <qualities (fasta format)>\n\n".
            "-- All sequences and qualities must be on a single line each. --.\n".
            "-- All quality values MUST be followed by a space character. --\n\n".
            "Merges fasta and qual files into a single Sanger fastq file (STDOUT).";

my $seqfile = shift or die $usage;
my $qualfile = shift or die $usage;

open F, "<$seqfile" or die $usage;
open Q, "<$qualfile" or die $usage;

my ($fh, $fs, $qh, $qs);
my @quals;  my $qv;

while (<F>) {
    ($fh) = $_ =~ /^(>\S+)/;
    $fs = <F>;
    ($qh) = <Q> =~  /^(>\S+)/;
    $qs = <Q>;
    if (($fh eq $qh) && ($fh =~ /^(\S+)/)) {
	my $hs = $1;
	my $hq = $1;
	$hs =~ s/>/\@/;
	print "$hs\n$fs";
	$hq =~ s/>/+/;
	print "$hq\n";
	@quals = split(/\s+/,$qs);
	for (my $i=0; $i<=$#quals; $i++) {
	    $qv = $quals[$i] + 33;
	    print chr($qv);
	}
	print "\n";
    }
    else { print "PROBLEM\n" }
}

close F; close Q;

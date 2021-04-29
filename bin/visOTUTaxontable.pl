#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;

my ($input1, $input2);
&GetOptions("i1=s" =>\$input1,
            "i2=s" =>\$input2
    );
($input1 && $input2) ||
    die "\nusage: $0\n".
    "-i1 <taxonomic annotation of each otu: *.taxonomy>\n".
    "-i2 <otu table file: *.OTU-table.txt>\n".
    "generate otu table summary file with taxonomic inforation\n";

#OTU_1|7275 Archaea(100);Euryarchaeota(100);Methanomicrobia(100);Methanomicrobiales(100);Methanoregulaceae(100);Methanoregula
open TAXON, "<$input1" or die "Could not open $input1 to read, $!\n";
my %otus = ();

while (<TAXON>) {
    chomp;
    my ($otuid, $count, $taxon) = $_ =~ /^OTU_?(\d+?)\|(\d+?)\s+?(\w.*)$/;
    $otus{$otuid}->{count} = $count;
    $otus{$otuid}->{taxon} = $taxon;
}
close(TAXON);

#OTUID 27 28 29
#OTU_1 2834 2553 1888

open TABLE, "<$input2" or die "Could not open $input2 to read, $!\n";

my %table = ();
my %sample_totals = ();
my @samples = ();
while (<TABLE>) {
    next if /^label/;
    chomp;
    my @line = split(/\s+/, $_);
    my $distance = shift @line;
    my $samplename = shift @line;
    my $numOtus = shift @line;
    my $id = 1;
    for(@line){
	my $num = shift @line;
	$table{$id}->{$samplename} = $num;
	$sample_totals{$samplename} +=  $num;
	$id++;
    }
    push(@samples, $samplename);
}

close(TABLE);


print_web_start();
print "<div>\n";
print "<h4>OTU abundant table (count): header is sortable</h4>\n";
print_countTable(\@samples, \%table, \%otus);
print "</div>\n";
print "<div style=\"height:30px;\"></div>\n";
print "<div>\n";
print "<hr>\n";
print "<h4>OTU abundant table (percentage): header is sortable</h4>\n";
print_percentTable(\@samples, \%sample_totals, \%table, \%otus);
print "</div>\n";
print_web_end();
sub print_web_end{
print  <<END;

 <div id="sep"></div>
 <div id="foot" align="center">Problems? Questions? Suggestions? Please contact <a href="mailto:xdong\@ucalgary.ca">Xiaoli Dong</a> or <a HREF="mailto:mstrous\@ucalgary.ca">Marc Strous</a><br><a href="http://www.ucalgary.ca/ebg/">Energy Bioengineering Group</a> in <a href="http://geoscience.ucalgary.ca/">Department of Geoscience</a> at <a href="http://www.ucalgary.ca">Univeristy of Calgary</a></div>
    
    </div>

    </body>
    </html>
END
}
sub print_web_start{
	print  <<HEAD;	
	<html lang=”en”>
   <head>
      <meta charset="utf-8">
      <title>MetaAmp Results Page</title>
      <link rel="stylesheet" href="../css/metaamp.css">
      <link rel="stylesheet" href="../css/jquery.dataTables.min.css">
      <link rel="stylesheet" href="../css/jquery-ui.css">
      <script src="../js/jquery-2.1.3.min.js"></script>
      <script src="../js/jquery-ui-1.11.4/jquery-ui.min.js"></script>
      <script src="../js/jquery.dataTables.min.js"></script>

      <script>
         \$(document).ready(function() {
         
        \$('#count').dataTable(
         { 
         bJQueryUI:false,
	"scrollX": true,
         sPaginationType: "full_numbers"
         
         });
	\$('#percent').dataTable(
         {
         bJQueryUI:false,
"scrollX": true,
         sPaginationType: "full_numbers"

         });

         
         });
      </script>
   </head>
<body bgcolor="#ffffff">
<div id="outform">
	<a href="index.html" id="logo">MetaAmp Logo</a>

    <h1>MetaAmp Version 2.0 </h1>
    <div id="sep"></div>

HEAD

}

sub print_countTable{

my ($samples, $table, $otus) = @_;
#print STDERR "this is @{$samples}\n";
print "<table id=\"count\" class=\"cell-border\">\n";
print  "<thead>\n";
#print "<caption>OTU-table.taxonomy (count)in sortable table</caption>\n";
print	"<tr>\n";
print "<th>OTUID</th>\n";
print "<th>Total</th>\n"; 
foreach (@$samples){
    print "<th>$_</th>\n";	
}
print "<th>Taxonomy</th>\n";
print "</tr>\n";
print "	</thead>\n";


print "	<tbody>\n";
print "<tr>\n";
foreach my $otuid (sort {$a <=> $b} keys %$table){
    print "<td>OTU_$otuid</td>\n" ;

#total of each otu
    print "<td>";
    print $otus->{$otuid}->{count} if exists  $otus->{$otuid}->{count}; 
    print "</td>\n";

    foreach  my $sample (@$samples){
        print "<td>";
	print $table->{$otuid}->{$sample} if exists $table->{$otuid}->{$sample};
	print "</td>\n";
    }
    print "<td>";
    print $otus->{$otuid}->{taxon} if exists $otus->{$otuid}->{taxon};
    print "</td>\n";
   print "</tr>\n";
}

print "</tbody>\n</table>\n";

}

sub print_percentTable{

my ($samples, $sample_totals, $table, $otus) = @_;

print "<table id=\"percent\" class=\"cell-border\">\n";
print  "<thead>\n";
#print "<caption>OTU-table.taxonomy (Percent)in sortable table</caption>\n";
print   "<tr>\n";
print "<th>OTUID</th>\n";
foreach (@$samples){
    print "<th>$_(\%)</th>\n";
}
print "<th>Taxonomy</th>\n";
print "</tr>\n";
print " </thead>\n";


print " <tbody>\n";
print "<tr>\n";
foreach my $otuid (sort {$a <=> $b} keys %$table){
    print "<td>OTU_$otuid</td>\n" ;

#total of each otu

    foreach  my $sample (@$samples){
	print "<td>", sprintf("%.4f", 100 * $table->{$otuid}->{$sample}/$sample_totals->{$sample}), "</td>\n";
    }
    print "<td>";
    print $otus->{$otuid}->{taxon} if exists $otus->{$otuid}->{taxon};
    print "</td>\n";
  print "</tr>\n";
}

print "</tbody>\n</table>\n";

}


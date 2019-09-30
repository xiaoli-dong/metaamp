#!/usr/local/bin/env perl

BEGIN {
    # MAGPIE is required for taxonomic rank determination
  $ENV{MAGPIEHOME} = "/export/home/xdong/deve/metagenome/magpie";
  push @INC, "$ENV{MAGPIEHOME}/lib";
}

use Taxonomy;

$#ARGV == -1
    and die "Usage: $0  <silva.nr_v123.fasta.tax>\n";

open INPUT, "$ARGV[0]" or die "Could not open $ARGV[0] to read, $!\n";

my %name2rank = ();
my $tax = new Taxonomy();
while(<INPUT>){
     if(/^#/){
	 print $_;
	 next;
     }
     s/\s+\[\S+\];/;/g;
     my @line = split("\t", $_);
     my $taxon = shift @line;
     my @names = split(/;/, $taxon);
     my $lineage = "";
     foreach my $name (@names){
	 my $rank = "";
	 if(not exists $name2rank{$name}){
	     my $taxida = $tax->name2taxid($name);
	     $rank = $tax->rank($taxida);
	     if(length($rank) == 0){
		 $rank = "NA";
	     }
	     else{
		 $name2rank{$name} = $rank;
	     }
	     $lineage .= "$rank\::$name;";
	 }
	 else{
	     $rank = $name2rank{$name};
	     $lineage .= "$rank\::$name;";
	 }
     }
     #my $c = () = $lineage =~ /phylum/g; 
     #if($c > 1){
     if($lineage =~ /phylum::Spirochaetae;phylum::Spirochaetes;/){
	 $lineage =~ s/phylum::Spirochaetae;phylum::Spirochaetes;/phylum::Spirochaetes;class::Spirochaetia;/g;
     }
     if($lineage =~ /phylum::Gemmatimonadetes;phylum::Gemmatimonadetes;/){
	 print STDERR $lineage, "\n";
	 $lineage =~ s/phylum::Gemmatimonadetes;phylum::Gemmatimonadetes;/phylum::Gemmatimonadetes;class::Gemmatimonadetes;/;
     }
     if($lineage =~ /phylum::Actinobacteria;phylum::Actinobacteria;/){
	 $lineage =~ s/phylum::Actinobacteria;phylum::Actinobacteria;/phylum::Actinobacteria;class::Actinobacteria;/;
	 
     }
     if($lineage =~ /phylum::Thermomicrobia;/){
	 #print STDERR $lineage, "\n";
	 $lineage =~ s/phylum::Thermomicrobia;/class\::Thermomicrobia;/g;
     }
     
     
     if($lineage =~ /phylum::Aquificae;phylum::Aquificae;/){
	 $lineage =~ s/phylum::Aquificae;phylum::Aquificae;/phylum::Aquificae;class::Aquificae;/;
     }
     if($lineage =~ /phylum::Chrysiogenetes;phylum::Chrysiogenetes;/){
	 $lineage =~ s/phylum::Chrysiogenetes;phylum::Chrysiogenetes;/phylum::Chrysiogenetes;class::Chrysiogenetes;/;
     }
     if($lineage =~ /phylum::Deferribacteres;phylum::Deferribacteres;/){
	 $lineage =~ s/phylum::Deferribacteres;phylum::Deferribacteres;/phylum::Deferribacteres;class::Deferribacteres;/;
     }
     if($lineage =~ /phylum::Elusimicrobia;phylum::Elusimicrobia;/){
	 $lineage =~ s/phylum::Elusimicrobia;phylum::Elusimicrobia;/phylum::Elusimicrobia;class::Elusimicrobia;/;
     }
     if($lineage =~ /phylum::Thermodesulfobacteria;phylum::Thermodesulfobacteria;/){
	 $lineage =~ s/phylum::Thermodesulfobacteria;phylum::Thermodesulfobacteria;/phylum::Thermodesulfobacteria;class::Thermodesulfobacteria;/;
     }
     if($lineage =~ /phylum::Thermotogae;phylum::Thermotogae;/){
	 
	 $lineage =~ s/phylum::Thermotogae;phylum::Thermotogae;/phylum::Thermotogae;class::Thermotogae;/;
     }
     
     #}
     print "$lineage\t", join("\t", @line);
}

close(INPUT);







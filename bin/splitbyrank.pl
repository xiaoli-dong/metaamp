#!/usr/local/bin/env perl


$#ARGV == -1
    and die "Usage: $0  <*.rank.txt><superkingdom|phylum|class|order|family|genus|species> <abund cutoff>\n";

open INPUT, "$ARGV[0]" or die "Could not open $ARGV[0] to read, $!\n";
my $rank = $ARGV[1];
my $cutoff = $ARGV[2];
my %total = ();
my @header = ();
my $header_line = "";
while(<INPUT>){
    my $domain = "NA";
    chomp;
    if(/^#/){
	print "$_\n";
	$header_line = $_;
	next;
    }
    if(scalar @header == 0){
	@header =  split(/\t/,$header_line);
	shift @header;
    }
    my @line = split(/\t/, $_);
      
    if(/superkingdom\:\:(\w+?);/){
	$domain = $1;
    
    }
   
    my $taxon = shift @line;

    if($taxon =~ /$rank\::(\S.*?);/){
	my $i = 0;
	foreach my $sp (@header){
	    if($rank ne $domain){
		#$total{"$domain\_$1"}->{$sp} += $line[$i++];
		$total{"$1"}->{$sp} += $line[$i++];
	    }
	    else{
		$total{$1}->{$sp} += $line[$i++];
	    }
	}
    }
}

close(INPUT);
#print "#taxon\t", join("\t", @header), "\n";
foreach my $term (keys %total){
    my $output = $term;
    my $pass = 0;
    foreach my $sp (@header){
	$output .= "\t". sprintf("%.2f", $total{$term}->{$sp});
	if($total{$term}->{$sp} >= $cutoff){
	    $pass = 1;
	}
    }
    if($pass == 1){
	print "$output\n";
    }
}






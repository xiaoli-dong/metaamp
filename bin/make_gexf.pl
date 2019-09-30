#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Long;

my ($assoFile, $sharedFile, $coef, $sig, $abound);
$coef = 0;
$sig = 0.05;
$abound = 0;

&GetOptions("a=s" =>\$assoFile,
	    "s=s" =>\$sharedFile,
	    "c=f" =>\$coef,
	    "p=f" =>\$sig,
	    "oa=f" =>\$abound
    );
($assoFile && $sharedFile) or
    die "\nusage: $0 \n".
    "-a <otu.association output file>\n".
    "-s <shared file producted from mothur>\n".
    "-c <correlation coefficient cutoff, the value should be between 0 and 1>\n".
    "-p <p value cutoff>\n".
    "-oa <otu size cutoff>\n";
    
    
open(SHARE, $sharedFile) or die "Could not open $sharedFile to read, $!\n";
my %otus = ();
my $total_reads = 0;
my %all_samples = ();
my %otu2samples = ();
my $total_test = 0;
while(<SHARE>){
    chomp;
    if(/label/){
	my @header = split(/\s+/, $_);
	shift @header;
	shift @header;
	shift @header;
	
	foreach my $a (0..$#header) {
	    foreach my $b ($a+1..$#header) {
		$total_test++;
	    }
	}
	next;
    }
    
    my @line = split(/\s+/, $_);
    my $dist = shift @line;
    my $sample = shift @line;
    shift @line;
    my $count = 1;
    foreach my $otu (@line){
	$otus{"Otu$count"} += $otu;
	$otu2samples{"Otu$count"}->{$sample} = $otu;
	if(not exists $all_samples{$sample}){
	    $all_samples{$sample} = 1;
	}
	$count++;
	$total_reads += $otu;
    }
}
close(SHARE);

my $assoFile_with_q = calc_qvalue($assoFile, $total_test, $sig);

open(ASSO, $assoFile_with_q) or die "Could not open $assoFile_with_q to read, $!\n";
my %edges = ();
my %filteredList = ();
while(<ASSO>){
    next if /Significance/;
    chomp;
    my @line = split(/\t/, $_);
    my $otuA = shift @line;
    my $otuB = shift @line;
    my $ab_coef = shift @line;
    my $ab_sig = shift @line;
    my $rank = shift @line;
    my $ab_qvalue = shift @line;
    
    if(exists $otus{$otuA} && exists $otus{$otuB}){
	#my $a_size = sprintf("%.3f", 100 * $otus{$otuA}/$total_reads);
	#my $b_size = sprintf("%.3f", 100 * $otus{$otuB}/$total_reads);
	#if(abs($ab_coef) >= $coef && $ab_qvalue <= $sig && $a_size >= $abound && $b_size >= $abound ){ 
	if(abs($ab_coef) >= $coef && $ab_qvalue <= $sig){ 
	    $edges{$otuA}->{$otuB}->{weight} = $ab_coef;
	    $edges{$otuA}->{$otuB}->{sig} = $ab_sig;
	    $edges{$otuA}->{$otuB}->{qvalue} = $ab_qvalue;
	    if(not exists $filteredList{$otuA}){
		$filteredList{$otuA} = 1;
	    }
	    if(not exists $filteredList{$otuB}){
		$filteredList{$otuB} = 1; 
	    }
	}
    }
}
close(ASSO);
    


print "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
print "<gexf xmlns=\"http://www.gexf.net/1.2draft\" xmlns:viz=\"http://www.gexf.net/1.2draft/viz\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:schemaLocation=\"http://www.gexf.net/1.2draft http://www.gexf.net/1.2draft/gexf.xsd\" version=\"1.2\">\n";
#print "<gexf xmlns=\"http://www.gexf.net/1.2draft\" version=\"1.2\">\n";

print "<graph mode=\"static\" defaultedgetype=\"undirected\">\n";
print "<attributes class=\"edge\">\n";
print "<attribute id=\"0\" title=\"q-value\" type=\"float\"/>\n";
print "<\/attributes>\n";

#nodes section

my %filtered_otus = ();
print "<nodes>\n";
my $nodes = "{\n\"nodes\":[\n";

foreach my $otu (sort{$filteredList{$b} <=> $filteredList{$a}} keys %filteredList){
    
    #my $size = sprintf("%.3f", 100 * $otus{$otu}/$total_reads);
    my $size = $otus{$otu};
    if(exists $filteredList{$otu}){
	print "\t<node id=\"$otu\" label=\"$otu\">\n";
	print "\t\t<viz:size value=\"$size\"/>\n";
	print "<viz:color r=\"255\" g=\"0\" b=\"255\"></viz:color>\n";
	print "\t</node>\n";
	$nodes .= "{\"id\": \"$otu\", \"group\": 1},\n";
    }
}
$nodes =~ s/,$//;
$nodes .= "],";

my $links = "\"links\": [\n";

foreach my $sample (keys %all_samples){
    
    #print "\t<node id=\"$sample\" label=\"$sample\" />\n";
}

print "</nodes>\n";

#edge section
print "<edges>\n";
my $edge_count = 0;
foreach my $otuA (sort keys %edges){
    
    foreach my $otuB (sort keys %{$edges{$otuA}}){
	
	my $weight = $edges{$otuA}->{$otuB}->{weight};
	my $size = sprintf("%.2f", abs($edges{$otuA}->{$otuB}->{weight}));
	my $sig = $edges{$otuA}->{$otuB}->{sig};
	my $qvalue = $edges{$otuA}->{$otuB}->{qvalue};
	print "\t<edge id=\"$edge_count\" source=\"$otuA\" target=\"$otuB\" weight=\"$size\">\n";
	if($weight > 0){
	    print "\t<viz:color r=\"0\" g=\"255\" b=\"0\"></viz:color>\n"
	}
	else{
	    print "\t<viz:color r=\"0\" g=\"255\" b=\"255\"></viz:color>\n"
	}
	print "\t\t<attvalues><attvalue for=\"0\" value=\"$qvalue\"/></attvalues>\n";
	print "\t</edge>\n";
	
	$links .= "{\"source\": \"$otuA\", \"target\": \"$otuB\", \"value\": $size},\n";
	
	$edge_count++;
	
    }
}
$links =~ s/,$//;
$links .= "]\n}\n";

foreach my $otu (sort keys %otu2samples){
    
    foreach my $sample (sort keys %{$otu2samples{$otu}}){
	
	#my $size = abs($edges{$otuA}->{$otuB}->{weight});
	#my $sig = $edges{$otuA}->{$otuB}->{sig};
	
	#print "\t<edge id=\"$edge_count\" source=\"$sample\" target=\"$otu\" />\n";
	
	$edge_count++;
	
    }
}




print "</edges>\n";

#end of the graph
print "</graph>\n";
print "</gexf>\n";


open JSON, ">test.json";
print JSON $nodes, $links;
close(JSON);

#P-values were corrected for multiple testing using the 
#False Discovery Rate approach (Benjamini and Hochberg, 1995)
#same p-value will give same rank and very importantly-increment 
#the rank by one for the following p-values
sub calc_qvalue{
    
    my ($assoFile, $total_test, $sig) = @_;
    my ($prefix) = $assoFile =~ /(^\S+?)\.otu.corr/;
    open(CORRECT, ">$prefix.qvalue.otu.corr");
    open(ASSO, $assoFile) or die "Could not open $assoFile to read, $!\n";
    my %matrix = ();
    my $header_line = "";
    while(<ASSO>){
	chomp;
	if(/Significance/){
	    $header_line = $_;
	    next;
	}
	my @line = split(/\t/, $_);
	$matrix{$_} = $line[3];
	
    }
    close(ASSO);
    my $prev_sig = -1;
    my $rank = 0;
    print CORRECT "$header_line\tRank\tq-value\n";
    foreach my $key (sort {$matrix{$a} <=> $matrix{$b}} keys %matrix){
	#print STDERR $key, "\n";
	my $value = $matrix{$key};
	if($value > $prev_sig){
	    $prev_sig = $value;
	    $rank++;
	}
	my $qvalue = sprintf("%.6f", $rank/$total_test);
	#if($qvalue > $sig){
	 #   last;
	#}
	print CORRECT "$key\t$rank\t$qvalue\n";
    }
    close(CORRECT);
    return "$prefix.qvalue.otu.corr";
}

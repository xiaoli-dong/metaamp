#!/usr/bin/env perl
use Getopt::Long;

my ($input, $depth);
$depth = 4;

&GetOptions("i=s" =>\$input,
	    "d=i" =>\$depth
            
    );

($input) || die "usage: $0 OPTIONS where options are: -i input file -d <tree depth, default is 4>\n";


my $leading = 1;
open(INPUT, "<$input") or die $!;

my %kid_parent = ();
my %taxon_terms = ();
my %level_1 = ();
my $root_count = 0;

#Taxonomy Total sample1 sample2
my $head = <INPUT>;
chomp $head;
my @samples = ();
if($head =~ /^#Taxonomy/){
    @samples = split(/\s+/, $head);
    
}
shift @samples;
shift @samples;

my %taxons = ();
while (<INPUT>){
    chmop;
    
    my @line = split(/\t/, $_);
    
    my $str = shift @line;
    my $count = shift @line;
        
    $taxons{$str}->{count} = @line;
    
    my @terms = split(/;/, $str);
    my $len = @terms;
    
    $root_count += $count;
    
    if(@terms && not exists $level_1{$terms[0]}){
	
	$level_1{$terms[0]} = 1;
	
	#print STDERR "$terms[0]\n";
    }
    #print STDERR "root_count = $root_count\n";
    
    for(my $i = 0; $i < $len - 1; $i++){
	
	$kid_parent {trim($terms[$i+1])} = trim($terms[$i]);
	
	for (my $index = 0; $index <@samples; $index++){
	    #print STDERR $samples[$index], "\t", $line[$index], "\t", $terms[$i], "\n";
	    $taxon_terms{trim($terms[$i])}->{$samples[$index]} += $line[$index];
	}
    }
    for (my $index = 0; $index <@samples; $index++){
	$taxon_terms{trim($terms[$len -1])}->{$samples[$index]} += $line[$index];
    }
}
close(INPUT);

foreach my $kid (keys %kid_parent){

    #print "kid=$kid, parent=", $kid_parent{$kid}, "\n";
}

foreach my $term (sort {$taxon_terms{$b} <=> $taxon_terms{$a}} keys %taxon_terms){
    #print "$term\t", $taxon_terms{$term}, "\n";
}

my %parent_kid = ();


foreach my $kid (keys %kid_parent){
    
    my $parent = $kid_parent{$kid};
    
    if(exists $parent_kid{$parent} && $parent_kid{$parent} !~ /$kid/){
	
	$parent_kid{$parent} .= "\t$kid";
    }
    else{
	$parent_kid{$parent} = "$kid";
    }
    #print STDERR "$parent\t$kid\n";
}


foreach (keys %parent_kid){

    #print "parent=$_\tkid=$parent_kid{$_}\t$taxon_terms{$_}\n";
}


print "{\"text\":\".\", \"children\": [\n";

#my $root_count = keys %level_1;
my $count = 0;
foreach (keys %level_1){
    $count++;
    #print STDERR "before getKids\n";
    getKids(1, $_, \%parent_kid, \%taxon_terms);
    if($count < $root_count){
	print "},\n";
    }
    else{
	print "}\n";
    }
}
print "]}\n";
sub getKids{

    my ($leading, $current, $parent_kid, $taxon_terms) = @_;
    
    #foreach my $current (keys %$level_1){
	#print STDERR "xiaoli $current\n";
    my $local_leading = $leading;
    
    print "{\nterm: \'$current\',\n";
    #print "desc: ", $taxon_terms->{$current}, ",\n";
    print "desc: \'desc\',\n";
    #my $root_count = $taxon_terms->{"1_root"};
    #my $current_count = $taxon_terms->{$current};
    #my $percent = sprintf("%.2f", $current_count*100/$root_count);
    #print STDERR "current=$current_count, root=$root_count, percent=$percent\n";
    #print "percent: $percent,\n";
    foreach (@samples){
	#print STDERR "current=$current\n";
	print "$_: ", $taxon_terms{trim($current)}->{$_}, ",\n"; 
    }
    
    my @kids = split (/\t/, $parent_kid->{$current});
    
    if (!@kids){
	print "iconCls:\'term\',\n";
	print "leaf:true\n";
	
	return;
    }
    if($leading < $depth){
	print "iconCls:\'term-folder\',\nexpanded: true,\n";
    }
    else{
	print "iconCls:\'term-folder\',\nexpanded: false,\n";
	
    }
    print "children:[\n";
    my $num = 0;
    foreach my $kid (@kids){
	$num++;
	getKids($leading+1, $kid, $parent_kid, $taxon_terms);
	if($num < @kids){
	    print "},";
	}
	else{
	    print "}";
	}
    }
    
    
    print "]\n";
    #}
}

#trim leading and trailing space
sub trim($)
{
    my $string = shift;
    $string =~ s/^\s+//;
    $string =~ s/\s+$//;
    $string =~ s/\'//g;
    $string =~ s/uncult ured/uncultured/g;
    $string =~ s/uncul tured/uncultured/g;
    $string =~ s/uncultu red/uncultured/g;
    $string =~ s/unculture d/uncultured/g;
    $string =~ s/;uncultured\s?;/;/g;
    $string =~ s/uncultured bacterium//g;
    
    
    return $string;
}

#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Long;
use Statistics::Descriptive;
use Number::Format qw(:subs);
use FindBin;
use Benchmark;
use File::Copy qw(move);
use Data::Dumper qw(Dumper);
use File::Copy;
use Time::Piece;
use Scalar::Util qw(openhandle);
my ($map,$geneType, $analysisname,$seqtype, $scutoff, $oligos, $minoverlen, $maxdiffs,$trunclen,$maxee, $seqformat, $pdiffs, $email);
$geneType = "rRNA";
$seqtype = "single";
$scutoff = 0.97;
$trunclen = 150;
$maxee = 1;
$seqformat = "fastq";
$minoverlen = 8;
$maxdiffs = 0;
$pdiffs = 0;
$email = "xdong\@ucalgary.ca";
my $quiet = 0;
my $cpus = 8;

&GetOptions("map=s" =>\$map,
	    "g=s" =>\$geneType,
	    "an=s" =>\$analysisname,
	    "seqtype=s" =>\$seqtype,
	    "seqformat=s" =>\$seqformat,
	    "s=f" =>\$scutoff,
	    "minoverlen=i" =>\$minoverlen,
	    "maxdiffs=i" =>\$maxdiffs,
	    "trunclen=i" =>\$trunclen,
	    "maxee=s" =>\$maxee,
	    "oligos=s" =>\$oligos,
	    "pdiffs=i" =>\$pdiffs,
	    "email=s" =>\$email,
	    "q=i" =>\$quiet
    );

($map && $analysisname && $oligos) or
    die "\nusage: $0 \n".
    "-map <reads mapping file>\n".
    "-g <maker gene type: rRNA|non-rRNA, default is rRNA\n".
    "-an <analysis name>\n".
    "-seqtype <paired|single>\n".
    "-seqformat <fastq|fasta>\n".
    "-s <similarity cutoff for generating otus, default is 0.97>\n".
    "-minoverlen <minimum overlap lenght while mergeing paired end reads, default is 8)\n".
    "-maxdiffs <maximum number of mismatches allowed in the overlap regions>\n".
    "-trunclen <number, trim the single end sequences to the fixed length>\n".
    "-maxee <float, discard reads with expected errors greater than this value>\n".
    "-oligos <oligo file>\n".
    "-pdiffs <default 0: maximum number of differences to the primer sequence>\n".
    "-email <default xdong\@ucalgary.ca>\n".
    "-q <quiet default 0>\n".
    "metaamp analysis pipeline: assume all the input reads are in fastq format\n";


# Record the command line
$analysisname =~ tr/[ \-]/_/;
my $params = "Parameters for running metaamp: $0 -map $map -g $geneType -an $analysisname -seqtype $seqtype -seqformat $seqformat -s $scutoff -minoverlen $minoverlen -maxdiffs $maxdiffs  -trunclen $trunclen -maxee $maxee -oligos $oligos -pdiffs $pdiffs -email $email\n";

my $logfile = "MetaAmp.$analysisname.log.txt";
open LOG, '>', $logfile or err("Can't open logfile");

my $qc_status_file = "MetaAmp.qcstatus.txt";
open STATUS, '>', $qc_status_file or err("Can't open $qc_status_file");

# Record process id in case we must kill this  job
msg("Command line:\n $params\n");
msg("Process ID = $$\n\n");


save_job_info("JobInfo.txt", $analysisname, $email);


my $bin = "$FindBin::RealBin";
my $projDir = $bin;
$projDir =~ s/\/bin//;
my $dbdir = "$projDir/database";
my $usearch = "$projDir/programs/usearch64";
my $mothur = "$projDir/programs/mothur/mothur";

my $t0 = Benchmark->new;
my $template = "$dbdir/silva.nr_v132.fasta";
my $taxonomy = "$dbdir/silva.nr_v132.tax";
my $dissimilarity = 1 - $scutoff;

#---------------- verify the mapping file input format---------------
my ($err, $mapHash) = verifyMappingFile($map, $seqtype);

#Error message in the mapping file
if($err ne ""){
    msg("Mapping file errors: $err\n");
    exit(4);
}

my $sampleStats;
($mapHash,$sampleStats) = qc($mapHash, $seqtype, $seqformat, $analysisname, $oligos);

if(keys %$mapHash == 0){
    msg("No sample passed QC stage in $analysisname job, exit\n");
    status("No sample passed QC stage in $analysisname job, exit\n");
    exit(3);
}

my @samplenames = keys %$mapHash;

my $cmd = "";
#capture stderr but discard stdout
$cmd = "$bin/make_clusters.pl -fasta $analysisname.qc.fasta -an $analysisname -s $scutoff -t $cpus";
runcmds($cpus, $cmd);

$cmd = "$bin/getMothurList.pl -i otu.map.uc -l $analysisname.list -g $analysisname.groups -ds $dissimilarity -fasta $analysisname.otus.fasta -rep $analysisname.OTUs";
runcmds($cpus, $cmd);

#using RDP classifier to classifier otus.
$cmd = "$mothur \"\#classify.seqs(fasta=$analysisname.OTUs.fasta,template=$template,taxonomy=$taxonomy,method=wang,ksize=8, iters=100, cutoff=80, processors=$cpus)\"  ";
runcmds($cpus, $cmd);

######################generate alpha and beta diversity################
$cmd = "$bin/get_diversity.pl -l $analysisname.list -g $analysisname.groups -t $cpus -p $mothur";
runcmds($cpus, $cmd);

if($geneType eq "rRNA"){
    get_taxon_info();
}



$cmd ="$bin/packing.pl -an $analysisname";
runcmds($cpus, $cmd);

build_result_web(\@samplenames,$sampleStats, $seqtype);
delfile(glob "mothur\.*");
my $t1 = Benchmark->new;
my $td = timediff($t1, $t0);
msg("metaamp took:" . timestr($td) . " to run\n");


sub verifyMappingFile{
    my ($map, $seqtype) = @_;
    #for the 454 reads using mothur make.fastq to covert them into fastq file
    open(MAPPING, "$map") or die "Could not open $map to read, $!\n";
    my %maps = ();

#Samplename             Treatment       Strand          Read1             Read2
    my %names = ();
    my %read1s = ();
    my %read2s = ();
    my $err = "";

    while(<MAPPING>){
	chomp;
	#skip the header of the mapping file
	next if /^#/;
	next if !/\S/;
	my @line = split(/\s+/, $_);

	if($seqtype eq "single" && @line < 4){
	    $err = "The mapping file for single-end data must have 4 colums";
	    last;
	}
	elsif($seqtype eq "paired" && @line < 5){
	    $err = "The mapping file for paired-end data must have 5 colums";
	    last;
	}

	if(exists $names{$line[0]}){
	    $err = "Samplename in mapping file must be unique: $line[0] has duplicate entries";
	    last;
	}

	if($line[2] !~ /[+-]/){
	    $err = "Sequnce strand can only be positive strand (+) or negative strand (-)";
	    last;
	}

	if(exists $read1s{$line[3]}){
	    $err = "Read entries in mapping file must be unique: $line[3] has duplicate entries in the mapping file";
	    last;
	}

	if($seqtype eq "paired" && exists $read2s{$line[4]}){
	    $err = "Read entries in mapping file must be unique: $line[4] has duplicate entries in the mapping file";
	    last;
	}

	if($err eq ""){

	    my $name = shift @line;
	    $name =~ tr/[ \-]/_/;
	    push(@{$maps{$name}}, @line);
	}
    }

    return ($err, \%maps);
}

sub qc{

    my ($mapHash, $seqtype, $seqformat, $analysisname, $oligos) = @_;
    my %sampleStats = ();

    #produce design file,rename input files, check whether the file existing ------------
    $mapHash = prepare_verify_inputs($mapHash, $seqtype, $seqformat);
    getRawStats($mapHash, $seqtype, \%sampleStats);

    if($seqtype eq "paired"){
	merging($mapHash, $seqtype, \%sampleStats);
    }
    if(keys %$mapHash > 0){
	trim_techSeq_addBarcodelabel($mapHash, $seqtype, $oligos, \%sampleStats);
    }

    if(keys %$mapHash > 0){
	qcFilter($mapHash, $seqtype ,$sampleStats);
    }

    if(keys %$mapHash > 0){
	pool_samples($seqtype, $analysisname, $mapHash, \%sampleStats);
    }
    return ($mapHash,\%sampleStats);
}

#----- check input file exists, s/ -/_/g in mapping, rename input file ----
sub prepare_verify_inputs{

    my ($mapHash, $seqtype, $seqformat) = @_;

    #generate design file for the sample comparsion and multivariate analysis
    open(DESIGN, ">design.txt") or die "Could not open design.txt to write, $!\n";
    my %new_mapHash = ();

    foreach my $samplename (keys %$mapHash){
	my $fastq = "";
	my $flip = "false";

	my @line = @{$mapHash->{$samplename}};
	#$samplename =~ tr/[ \-]/_/;
	my $treatment = $line[0];
	$treatment =~ tr/[ \-]/_/;
	push(@{$new_mapHash{$samplename}}, $treatment);

	my $strand = $line[1];
	push(@{$new_mapHash{$samplename}}, $strand);
	my $read1 = $line[2];
	my $read2 = "";
	die "Read1 $read1 does not exist, $!\n" if ! -e $read1;
	my $read1_old = $read1;
	$read1 =~ tr/[ \-]/_/;
	move $read1_old, $read1;
	push(@{$new_mapHash{$samplename}}, $read1);
	if($seqtype eq "paired"){
	    $read2 = $line[3];
	    die "Read2 $read2 does not exist, $!\n" if ! -e $read2;
	    my $read2_old = $read2;
	    $read2 =~ tr/[ \-]/_/;
	    move $read2_old, $read2;
	    push(@{$new_mapHash{$samplename}}, $read2);
	}

        #convert fasta to fastq
	if($seqformat eq "fasta"){
	    $read1 = fasta2fastq($read1);
	    $new_mapHash{$samplename}[2] = $read1;
	    if($seqtype eq "paired"){
		$read2 = fasta2fastq($read2);
		$new_mapHash{$samplename}[3] = $read2;
	    }
	}
	print DESIGN "$samplename\t$treatment\n";
    }
    close(DESIGN);

    return (\%new_mapHash);
}


sub getRawStats{

    my ($mapHash, $seqtype, $sampleStats) = @_;

    foreach my $samplename (keys %$mapHash){
	##############stat for raw input reads###################
	if($seqtype eq "paired"){
	    my $read1 = $mapHash->{$samplename}[2];
	    my $read2 = $mapHash->{$samplename}[3];
	    $sampleStats->{$samplename}->{raw}->{read1} = fastqStats($read1);
	    $sampleStats->{$samplename}->{raw}->{read2} = fastqStats($read2);
	}
	else{
	    my $read1 = $mapHash->{$samplename}[2];
	    #do nothing for unpaired sample
	    $sampleStats->{$samplename}->{raw}->{read1} = fastqStats($read1);
	}
    }

}


##################### do merging for paired end ################
sub merging{

    my ($mapHash, $seqtype, $sampleStats) = @_;

    foreach my $samplename (keys %$mapHash){
	my $read1 = $mapHash->{$samplename}[2];
	my $read2 = $mapHash->{$samplename}[3];
	my $cmd = "$usearch -fastq_mergepairs $read1 -reverse $read2  -fastq_minlen 50 -fastq_maxdiffs $maxdiffs -fastq_minovlen $minoverlen -fastqout $samplename.merged.fastq -threads $cpus -log $samplename.merged.log;";
	msg("Running:", $cmd);
	system($cmd)==0 or err("Could not run command:", $cmd);
    }

    my @excludes = ();

    foreach my $samplename (keys %$mapHash){
	my $fastq = "$samplename.merged.fastq";
	$sampleStats->{$samplename}->{merged} = fastqStats("$samplename.merged.fastq");
	if(-z $fastq || ! -e $fastq){
	    $sampleStats->{$samplename}->{merged} = dummyFastqStats();
	    push(@excludes, $samplename);
	}
    }

    foreach (@excludes){
	msg("******Sample $_: no reads passing QC_merging, exclude it from further analysis\n");
	status("******Sample $_: no reads can be merged together (overlap length shorter than $minoverlen cutoff or the number of the mismatches in the overlap regions are more than $maxdiffs cutoff, exclude it from further analysis");
	delete($mapHash->{$_});
    }


}

sub trim_techSeq_addBarcodelabel{

    my ($mapHash, $seqtype, $oligos, $sampleStats) = @_;

    my @excludes = ();
    foreach my $samplename (keys %$mapHash){
	my $fastq = "";
	if($seqtype eq "paired"){
	    $fastq = "$samplename.merged.fastq";
	}
	else{
	    $fastq = $mapHash->{$samplename}[2];
	}
	if(-s $fastq){
	    my ($prefix) = $fastq =~ /(.*)\.[^.]+$/;

	    if( -s $oligos){
		my $cmd = "$bin/fastq_trimPrimer.pl -fastq $fastq -oligos $oligos -s $samplename -pdiffs $pdiffs -t $cpus;";
		$cmd .= "$bin/addBarcodelabel.pl  $prefix.trim.fastq $samplename > $samplename.strip.fastq;";
		msg("Running:", $cmd);
		system($cmd)==0 or err("Could not run command:", $cmd);
	    }
	    else{
		my $cmd = "$bin/addBarcodelabel.pl $fastq $samplename > $samplename.strip.fastq";
		msg("Running:", $cmd);
		system($cmd)==0 or err("Could not run command:", $cmd);
	    }
	}
    }


    @excludes = ();
    foreach my $samplename (keys %$mapHash){
	my $fastq = "$samplename.strip.fastq";
	$sampleStats->{$samplename}->{nonbarcodeprimer} = fastqStats("$samplename.strip.fastq");
	if(-s $fastq){
	    my $cmd = "$usearch -fastq_stats $fastq -threads $cpus -log $samplename.strip.stats;";
	    msg("Running:", $cmd);
	    system($cmd)==0 or err("Could not run command:", $cmd);
	}
	else{
	    $sampleStats->{$samplename}->{nonbarcodeprimer} = dummyFastqStats();
	    push(@excludes, $samplename);
	}
    }

    foreach (@excludes){
	msg("******Sample $_: no reads passing QC primer matching stage, exclude it from further analysis\n");
	status("******Sample $_: no reads passed primer matching and trimming stage, exclude it from further analysis");
	delete($mapHash->{$_});
    }

}


sub qcFilter{
    my ($mapHash, $seqtype, $sampleStats) = @_;

    my @excludes = ();

    foreach my $samplename (keys %$mapHash){
	my $fastq = "$samplename.strip.fastq";
	my $cmd = "$usearch -fastq_filter $fastq -threads $cpus -fastq_trunclen $trunclen -fastq_maxee $maxee -fastaout $samplename.filtered.fasta;";
	msg("Running:", $cmd);
	system($cmd)==0 or err("Could not run command:", $cmd);
    }

    foreach my $samplename (keys %$mapHash){
	my $strand = $mapHash->{$samplename}[1];
	my $flip = "false";
	if($strand eq "-"){
	    $flip = "T";
	}
	my $prefix = "$samplename.filtered";
	if(-s "$prefix.fasta"){
	    ###################do not do flip till finish qc#################
	    #$cmd = "$mothur \"#fastq.info(fastq=$prefix.fastq);trim.seqs(fasta=$prefix.fasta, qfile=$prefix.qual, flip=$flip)\" ";
	    my $cmd = "$mothur \"#trim.seqs(fasta=$prefix.fasta, flip=$flip,processors=$cpus)\" ";
	    msg("Running:", $cmd);
	    system($cmd)==0 or err("Could not run command:", $cmd);
	}
	else{
	    push(@excludes, $samplename);
	    next;
	}
    }

    foreach (@excludes){
	msg("******Sample $_: no reads passing QCFilter, exclude it from further analysis\n");
	status("******Sample $_: no reads passed quality filter stage (It can be the Maximum number of expected errors is greater than $maxee cutoff or trimmed amplicon length is shorter than $trunclen cutoff) exclude it from further analysis");
	delete $mapHash->{$_};
    }
    #return $mapHash;


    @excludes = ();
    foreach my $samplename (keys %$mapHash){
	my $fasta = "$samplename.filtered.trim.fasta";
	$sampleStats->{$samplename}->{qc} = fastaStats($fasta);
	if(-z $fasta){
	    $sampleStats->{$samplename}->{qc} = dummyFastqStats($fasta);
	    push(@excludes, $samplename);
	}
    }
    foreach (@excludes){
	msg("******Sample $_: no reads passing qcFilter , exclude it from further analysis\n");
	status("******Sample $_: no reads passed quality filter stage (It can be the Maximum number of expected errors is greater than $maxee cutoff or the primer trimmed amplicon length did not reach the minimum lenght cutoff $trunclen) exclude it from further analysis");
	delete($mapHash->{$_});
    }

}

sub pool_samples{
    my ($seqType, $analysisname, $mapHash, $sampleStats) = @_;
    print STDERR "pool samples\n";
# Erase the contents if already exists
    if(-e "$analysisname.qc.fasta"){
        unlink("$analysisname.qc.fasta") or die "Could not delete the file: $analysisname.qc.fasta, $!\n";
    }
    if(-e "$analysisname.qc.qual"){
        unlink("$analysisname.qc.qual") or die "Could not delete the file: $analysisname.qc.qual, $!\n";
    }
    if(-e "$analysisname.qc.fastq"){
        unlink("$analysisname.qc.fastq") or die "Could not delete the file: $analysisname.qc.fastq, $!\n";
    }

    foreach my $samplename (keys %$mapHash){
	if($seqtype eq "paired"){
	    $sampleStats->{$samplename}->{merged} = fastqStats("$samplename.merged.fastq");
	}
	$sampleStats->{$samplename}->{nonbarcodeprimer} = fastqStats("$samplename.strip.fastq");
	$sampleStats->{$samplename}->{qc} = fastaStats("$samplename.filtered.trim.fasta");
	print STDERR "pool samples outside $samplename.filtered_reads.trim.fasta\n";

	if(-s "$samplename.filtered.trim.fasta"){
	    print STDERR "pool samples inside\n";
            my $cmd =  "cat $samplename.filtered.trim.fasta >> $analysisname.qc.fasta;";
	    msg("Running:", $cmd);
	    system($cmd)==0 or err("Could not run command:", $cmd);

        }

    }
    $sampleStats->{"POOL"}->{qc} = fastaStats("$analysisname.qc.fasta");
    return $sampleStats;
}

sub fasta2fastq {

    my($fasta) = @_;
    my ($fasta_prefix) = $fasta =~ /(.*)\.[^.]+$/;
    my $fastq = "$fasta_prefix.fastq";

    if( -e "$fasta_prefix\.qual"){
	my $fasta_qual = "$fasta_prefix.qual";
	my $cmd = "$mothur \"#make.fastq(fasta=$fasta, qfile=$fasta_qual)\"";
	runcmds($cpus, $cmd);
    }
    #no quality value, make psedo quality value 40
    else{

	$/ = "\n>";
	open(FASTA, $fasta) or die "Could not open $fasta to read, $!\n";
	open(FASTQ, ">$fastq") or die "Could not open $fastq file to write, $!\n";

	while(<FASTA>){
	    chomp;
	    if(!/>?(\S.*?)\n(.+)/s){
		die "Could not read FastA record #$.: $_\n";
	    }
	    my $head = $1;
	    my $seq = $2;
	    $seq =~ s/[ \r\n\t]//g;

	    my $len = length($seq);
	    print FASTQ "\@$head\n$seq\n+\n";

	    for (my $i=0; $i < $len ; $i++) {
		print FASTQ "I";
	    }
	    print FASTQ "\n";

	}

	$/ = "\n";
	close(FASTQ);
	close(FASTA);
    }
    return $fastq;
}

sub dummyFastqStats{
    my %stat_info = ();
    $stat_info{min} = 0;
    $stat_info{max} = 0;
    $stat_info{mean} = 0;
    $stat_info{median} = 0;
    $stat_info{mode} = 0;
    $stat_info{count} = 0;
    $stat_info{sum} = 0;
    $stat_info{std} = 0;

    return \%stat_info;

}
sub fastaStats{

    my ($fasta) = @_;
    my %stat_info = ();

    open(FASTA, "$fasta") or die "Could not open $fasta file to read, $!\n";
    #open(HIST, ">$fasta.hist.txt") or die "Could not open $fasta.hist.txt file to write, $!\n";

    my @lens = ();

    $/ = "\n>";
    while(<FASTA>){
	chomp;
	if(my ($seq_name,$seq) =  /^>?(\S+.*?)\n(.*)/s){
	    $seq =~ tr/ \r\n\t//d;
	    push(@lens, length($seq));
	}
    }

    $/ = "\n";

    close(FASTA);

    $Statistics::Descriptive::Tolerance = 1e-10;
    my $stat = Statistics::Descriptive::Full->new();
    #$stat->add_data(1,2,3,4);
    $stat->add_data(@lens);


    $stat_info{min} = $stat->min();
    $stat_info{max} = $stat->max();
    $stat_info{mean} = $stat->mean();
    $stat_info{median} = $stat->median();
    $stat_info{mode} = $stat->mode();
    $stat_info{count} = $stat->count();
    $stat_info{sum} = $stat->sum();
    $stat_info{std} = $stat->standard_deviation();



#get histogram of the sequence length distribution
    my @bins = ();

    for(my $i = $stat_info{min}; $i <= $stat_info{max}; $i++){
	push(@bins, $i);
    }

    #$stat->frequency_distribution_ref(\@bins);
    #my $f = $stat->frequency_distribution_ref(\@bins);

    #print HIST "Len,Count\n";
    #for (sort {$a <=> $b} keys %$f) {
    #if($f->{$_} > 0){
#	    print HIST "$_,$f->{$_}\n";
#	}
 #   }
  #  close(HIST);
    return \%stat_info;
}


sub fastqStats{

    my ($fastq) = @_;
    my %stat_info = ();

    open(FASTQ, "$fastq") or die "Could not open $fastq file to read, $!\n";
    #open(HIST, ">$fastq.hist.txt") or die "Could not open $fastq.hist.txt file to write, $!\n";

    my @lens = ();
    while (<FASTQ>){
	my $line1 = $_;
	my $line2 = <FASTQ>;
	my $line3 = <FASTQ>;
	my $line4 = <FASTQ>;

	if(defined $line2){
	    chomp $line2;
	    $line2 =~ s/\s+//g;
	    push(@lens, length($line2));
	}
    }
    close(FASTQ);


    $Statistics::Descriptive::Tolerance = 1e-10;
    my $stat = Statistics::Descriptive::Full->new();
    #$stat->add_data(1,2,3,4);
    $stat->add_data(@lens);


    $stat_info{min} = $stat->min();
    $stat_info{max} = $stat->max();
    $stat_info{mean} = $stat->mean();
    $stat_info{median} = $stat->median();
    $stat_info{mode} = $stat->mode();
    $stat_info{count} = $stat->count();
    $stat_info{sum} = $stat->sum();
    $stat_info{std} = $stat->standard_deviation();



#get histogram of the sequence length distribution
    my @bins = ();
    if(exists $stat_info{min} && exists $stat_info{max}){
	for(my $i = $stat_info{min}; $i <= $stat_info{max}; $i++){
	    push(@bins, $i);
	}
    }
    #$stat->frequency_distribution_ref(\@bins);
    #my $f = $stat->frequency_distribution_ref(\@bins);

    #print HIST "Len,Count\n";
    #for (sort {$a <=> $b} keys %$f) {
	#if($f->{$_} > 0){
	    #print HIST "key = $_, count = $f->{$_}\n";
	    #print HIST "$_,$f->{$_}\n";
	#}
    #}
    #close(HIST);
    return \%stat_info;
}

sub get_taxon_info{

    open(WANG, "$analysisname.OTUs.nr_v132.wang.taxonomy") or die "Could not open $analysisname.OTUs.nr_v132.wang.taxonomy to read, $!\n";
    open(MWANG, ">$analysisname.OTUs.nr_v132.mwang.taxonomy") or die "Could not open $analysisname.OTUs.nr_v132.mwang.taxonomy to write, $!\n";
    open(MWANGNR, ">$analysisname.OTUs.nr_v132.mwang.norank.taxonomy") or die "Could not open $analysisname.OTUs.nr_v132.mwang.norank.taxonomy to write, $!\n";
    while(<WANG>){
	my ($id, $tax) = $_ =~ /^(\S+?)\s+(\S.*)$/;
	my $new_tax = "";
	foreach my $tmp (split(/;/, $tax)){
	    next if $tmp =~ /unclassified/i;
	    $new_tax .= "$tmp;";
	}
	print MWANG "$id\t$new_tax\n";
	$new_tax =~ s/\;?\w+\__/;/g;
	$new_tax =~ s/^;//;
	print MWANGNR "$id\t$new_tax\n";

	#s/unclassified;//g;
	#print MWANG $_;
	#s/\t\S+?\__/\t/g;
	#s/\;?\w+\__/;/g;
	#print MWANGNR $_;
    }
    close(WANG);
    close(MWANG);
    close(MWANGNR);


######################get analysis otu table and taxonomy summary table################
    $cmd = "$bin/getOTUTaxonTable.pl -i1 $analysisname.OTUs.nr_v132.mwang.norank.taxonomy  -i2 $analysisname.shared > $analysisname.OTU-table.taxonomy";
    runcmds($cpus, $cmd);

######################get analysis otu table and taxonomy summary table################
    $cmd = "$bin/getTaxonSummary.pl -i1 $analysisname.OTUs.nr_v132.mwang.norank.taxonomy  -i2 $analysisname.shared > $analysisname.taxonomy.norank.summary";
    runcmds($cpus, $cmd);


######################get analysis otu table and taxonomy summary table without rank ################
    $cmd = "$bin/getTaxonSummary.pl -i1 $analysisname.OTUs.nr_v132.mwang.taxonomy  -i2 $analysisname.shared > $analysisname.taxonomy.summary";
    runcmds($cpus, $cmd);
}



sub build_result_web{
    my ($samplenames, $stats, $seqtype) = @_;

    open(JOB_INFO, "JobInfo.txt")
	or die "Can't open job information file for reading: $!";

    chomp(my $analysis = <JOB_INFO>);
    chomp(my $email = <JOB_INFO>);
    chomp(my $human_time = <JOB_INFO>);
    chomp(my $machine_time = <JOB_INFO>);
    close JOB_INFO;
    umask 022;
    my $results_dir = "results";
    my $results_link = "$results_dir";
    chdir $results_dir
	or die "Can't change to $results_dir directory: $!";
    open my $RESULTS, ">index.html"
	or die "Can't open $results_dir/index.html to write: $!";
    print_web_head($RESULTS);
    print_web_common($RESULTS);

    print_web_download($RESULTS, $analysis, $email, $human_time, $results_link);
    print_web_QC($RESULTS, $stats, $seqtype, $analysis);

    print_web_otus($RESULTS, $analysis);
    print_web_alpha_diversity($RESULTS, $analysis, $samplenames);
    #print_web_beta_diversity($RESULTS, $analysis);
    print_hypothesis_testing($RESULTS, $analysis, $samplenames);
    print_web_common_end($RESULTS, $analysis);
    chdir "..";
    my $cmd = "zip -q $results_link.zip -r $results_link";

    runcmds($cpus, $cmd);
    #chdir $results_dir;
    send_job_finished_email($analysis, $email, $human_time, $machine_time, $results_dir);

}
########## Subroutines ##########

# Sends an email to the user notifying the job is finished
sub send_job_finished_email {
  	my ($analysis, $email, $human_time, $machine_time, $results_dir) = @_;
	open(MAIL, "|/usr/sbin/sendmail -t");
  print MAIL "To: $email\n";
  #print MAIL "Bcc: xdong\@ucalgary.ca\n";
  print MAIL "From: metaAmp\@ebg.ucalgary.ca\n";
  print MAIL "Subject: metaAmp job done: $analysis\n\n";
  print MAIL "Dear metaAmp user,\n\nYour metaAmp job with the analysis name: $analysis has been finished. Please visit the link below to view or download your analysis results:\n\nhttp://ebg.ucalgary.ca/metaamp/tmp/$machine_time-$analysis/$results_dir/index.html\n\nThank you for using metaAmp.\n";
  close(MAIL);
}


sub print_web_head {
    my ($RESULTS) = @_;
    print $RESULTS <<RESULTS_HEAD;
<html lang=”en”>
   <head>
      <meta charset="utf-8">
      <title>MetaAmp Results Page</title>
      <link rel="stylesheet" href="http://ebg.ucalgary.ca/metaamp/css/metaamp.css">
      <script src="http://code.jquery.com/jquery-3.1.1.min.js"></script>
      <script src="http://code.jquery.com/ui/1.12.1/jquery-ui.min.js"></script>
      <link rel="stylesheet" href="http://code.jquery.com/ui/1.12.1/themes/base/jquery-ui.css">
      <script src="http://cdn.datatables.net/1.10.13/js/jquery.dataTables.min.js"></script>
      <link rel="stylesheet" href="http://cdn.datatables.net/1.10.13/css/jquery.dataTables.min.css">
      <script>
      function help(divid, btnid, closeid){
	  // Get the modal
	      var modal = document.getElementById(divid);
	  modal.style.display = "block";

	  // Get the button that opens the modal
	      var btn = document.getElementById(btnid);

	  // Get the <span> element that closes the modal
	      var span = document.getElementById(closeid);

	  // When the user clicks the button, open the modal
	      btn.onclick = function() {
		  modal.style.display = "block";
	  }

	  // When the user clicks on <span> (x), close the modal
	      span.onclick = function() {
		  modal.style.display = "none";
	  }

	  // When the user clicks anywhere outside of the modal, close it
	      window.onclick = function(event) {
		  if (event.target == modal) {
			  modal.style.display = "none";
		  }
	  }
  }
      </script>

      <script>
         \$(document).ready(function() {
         \$('article.tabs').tabs();

         \$('#qcstats').dataTable(
         {
         bJQueryUI:false,
         sPaginationType: "full_numbers"

         });

	 \$('#trimstats').dataTable(
         {

         bJQueryUI:false,
         sPaginationType: "full_numbers"

         });
	 \$('#mergestats').dataTable(
         {

         bJQueryUI:false,
         sPaginationType: "full_numbers"

         });

         \$('#rawstats').dataTable(
         {


         bJQueryUI:false,
         sPaginationType: "full_numbers"

         });
         });
      </script>
   </head>

RESULTS_HEAD
}

sub print_web_common{

    my($RESULTS) = @_;
    print $RESULTS <<RESULTS_COMMON;
    <body>
	<div id="outform">
	<a href="/metaamp/index.html" id="logo">MetaAmp Logo</a>

	<h1>MetaAmp Version 2.0 Results</h1>
	<div id="sep"></div>
	<div style="height: 20px; "></div>
	<div>
	<p>Click the Download tab to download the results to browse it offline. We only keeps the results on the server for a week</p>
	<p>Or click the other tabs to browse the result online</b></p>
      </div>
      <section>
      <article id="tabs-nohdr" class="tabs">
      <ul>
      <li><a href="#tabs-nohdr-1">Download</a></li>
      <li><a href="#tabs-nohdr-2">QC Summary</a></li>
      <li><a href="#tabs-nohdr-3">OTU & Taxonomy</a></li>
      <li><a href="#tabs-nohdr-5">&alpha;,&beta;-Diversity</a></li>
      <li><a href="#tabs-nohdr-6">Hypothesis Testing</a></li>

      </ul>
RESULTS_COMMON
}

sub print_web_download{
    my ($RESULTS, $analysis, $email, $human_time, $results_link) = @_;
    print $RESULTS <<RESULTS_DOWNLOAD;
    <div id="tabs-nohdr-1">
	<p><b>MetaAmp analysis results:</b></p>
	<p> Analysis name: $analysis</p>
	<p>$params</p>
	<p>Submitted by $email at $human_time</p>
	<p class="thick"><a href="../$results_link.zip">Download a packaged file containing all results in this page</a></p>
	</div>

RESULTS_DOWNLOAD
}

sub print_web_QC{
    my($RESULTS, $stats, $seqtype, $analysisname) = @_;

    open(QCSUMMARY, ">QC/$analysisname.qc.summary.txt") or die "Could not open QC/$analysisname.qc.summary.txt to write, $!\n";
    if($seqtype eq "paired"){
	print QCSUMMARY "#samplename\tR1_count\tR1_min\tR1_max\tR1_mean\tR1_median\tR1_std\tR1_sum\tR2_count\tR2_min\tR2_max\tR2_mean\tR2_median\tR2_std\tR2_sum\tmerged_count\tmerged_percent\%\tmerged_min\tmerged_max\tmerged_mean\tmerged_median\tmerged_std\tmerged_sum\tnonbarcodeprimer_count\tnonbarcodeprimer_percent\%\tnonbarcodeprimer_min\tnonbarcodeprimer_max\tnonbarcodeprimer_mean\tnonbarcodeprimer_median\tnonbarcodeprimer_std\tnonbarcodeprimer_sum\tQC_count\tQC_percent\%\tQC_min\tQC_max\tQC_mean\tQC_median\tQC_std\tQC_sum\n";}
    else{

	print QCSUMMARY "#samplename\tR1_count\tR1_min\tR1_max\tR1_mean\tR1_median\tR1_std\tR1_sum\tnonbarcodeprimer_count\tnonbarcodeprimer_percent\%\tnonbarcodeprimer_min\tnonbarcodeprimer_max\tnonbarcodeprimer_mean\tnonbarcodeprimer_median\tnonbarcodeprimer_std\tnonbarcodeprimer_sum\tQC_count\tQC_percent\%\tQC_min\tQC_max\tQC_mean\tQC_median\tQC_std\tQC_sum\n";}

    my $totalRawCount_read1 = 0;
    foreach my $key (%$stats){
	if(exists $stats->{$key} && $key ne "POOL"){
	    my $rawStat = $stats->{$key}->{raw};
	    $totalRawCount_read1 += $sampleStats->{$key}->{raw}->{read1}->{count};

	    if($seqtype eq "paired"){
		print QCSUMMARY "$key\t", $sampleStats->{$key}->{raw}->{read1}->{count}, "\t";
		print QCSUMMARY $sampleStats->{$key}->{raw}->{read1}->{min}, "\t";
		print QCSUMMARY $sampleStats->{$key}->{raw}->{read1}->{max}, "\t";
		print QCSUMMARY $sampleStats->{$key}->{raw}->{read1}->{mean}, "\t";
		print QCSUMMARY $sampleStats->{$key}->{raw}->{read1}->{median}, "\t";
		print QCSUMMARY $sampleStats->{$key}->{raw}->{read1}->{std}, "\t";
		print QCSUMMARY $sampleStats->{$key}->{raw}->{read1}->{sum}, "\t";

		print QCSUMMARY $sampleStats->{$key}->{raw}->{read2}->{count}, "\t";
		print QCSUMMARY $sampleStats->{$key}->{raw}->{read2}->{min}, "\t";
		print QCSUMMARY $sampleStats->{$key}->{raw}->{read2}->{max}, "\t";
		print QCSUMMARY $sampleStats->{$key}->{raw}->{read2}->{mean}, "\t";
		print QCSUMMARY $sampleStats->{$key}->{raw}->{read2}->{median}, "\t";
		print QCSUMMARY $sampleStats->{$key}->{raw}->{read2}->{std}, "\t";
		print QCSUMMARY $sampleStats->{$key}->{raw}->{read2}->{sum}, "\t";

		print QCSUMMARY $sampleStats->{$key}->{merged}->{count}, "\t";
		my $percent = sprintf("%.2f", 100*$sampleStats->{$key}->{merged}->{count}/$sampleStats->{$key}->{raw}->{read1}->{count});
		print QCSUMMARY "$percent\t";
		print QCSUMMARY $sampleStats->{$key}->{merged}->{min}, "\t";
		print QCSUMMARY $sampleStats->{$key}->{merged}->{max}, "\t";
		print QCSUMMARY $sampleStats->{$key}->{merged}->{mean}, "\t";
		print QCSUMMARY $sampleStats->{$key}->{merged}->{median}, "\t";
		print QCSUMMARY $sampleStats->{$key}->{merged}->{std}, "\t";
		print QCSUMMARY $sampleStats->{$key}->{merged}->{sum}, "\t";

		$percent = sprintf("%.2f", 100*$sampleStats->{$key}->{nonbarcodeprimer}->{count}/$sampleStats->{$key}->{raw}->{read1}->{count});

		print QCSUMMARY $sampleStats->{$key}->{nonbarcodeprimer}->{count}, "\t";
		print QCSUMMARY "$percent\t";
		print QCSUMMARY $sampleStats->{$key}->{nonbarcodeprimer}->{min}, "\t";
		print QCSUMMARY $sampleStats->{$key}->{nonbarcodeprimer}->{max}, "\t";
		print QCSUMMARY $sampleStats->{$key}->{nonbarcodeprimer}->{mean}, "\t";
		print QCSUMMARY $sampleStats->{$key}->{nonbarcodeprimer}->{median}, "\t";
		print QCSUMMARY $sampleStats->{$key}->{nonbarcodeprimer}->{std}, "\t";
		print QCSUMMARY $sampleStats->{$key}->{nonbarcodeprimer}->{sum}, "\t";

		print QCSUMMARY $sampleStats->{$key}->{qc}->{count}, "\t";
		$percent = sprintf("%.2f", 100*$sampleStats->{$key}->{qc}->{count}/$sampleStats->{$key}->{raw}->{read1}->{count});
		print QCSUMMARY "$percent\t";
		print QCSUMMARY $sampleStats->{$key}->{qc}->{min}, "\t";
		print QCSUMMARY $sampleStats->{$key}->{qc}->{max}, "\t";
		print QCSUMMARY $sampleStats->{$key}->{qc}->{mean}, "\t";
		print QCSUMMARY $sampleStats->{$key}->{qc}->{median}, "\t";
		print QCSUMMARY $sampleStats->{$key}->{qc}->{std}, "\t";
		print QCSUMMARY $sampleStats->{$key}->{qc}->{sum}, "\t";
		print QCSUMMARY "\n";
	    }

	    elsif($seqtype eq "single"){

		print QCSUMMARY "$key\t", $sampleStats->{$key}->{raw}->{read1}->{count}, "\t";
		print QCSUMMARY $sampleStats->{$key}->{raw}->{read1}->{min}, "\t";
		print QCSUMMARY $sampleStats->{$key}->{raw}->{read1}->{max}, "\t";
		print QCSUMMARY $sampleStats->{$key}->{raw}->{read1}->{mean}, "\t";
		print QCSUMMARY $sampleStats->{$key}->{raw}->{read1}->{median}, "\t";
		print QCSUMMARY $sampleStats->{$key}->{raw}->{read1}->{std}, "\t";
		print QCSUMMARY $sampleStats->{$key}->{raw}->{read1}->{sum}, "\t";



		my $percent = sprintf("%.2f", 100*$sampleStats->{$key}->{nonbarcodeprimer}->{count}/$sampleStats->{$key}->{raw}->{read1}->{count});

		print QCSUMMARY $sampleStats->{$key}->{nonbarcodeprimer}->{count}, "\t";
		print QCSUMMARY "$percent\t";
		print QCSUMMARY $sampleStats->{$key}->{nonbarcodeprimer}->{min}, "\t";
		print QCSUMMARY $sampleStats->{$key}->{nonbarcodeprimer}->{max}, "\t";
		print QCSUMMARY $sampleStats->{$key}->{nonbarcodeprimer}->{mean}, "\t";
		print QCSUMMARY $sampleStats->{$key}->{nonbarcodeprimer}->{median}, "\t";
		print QCSUMMARY $sampleStats->{$key}->{nonbarcodeprimer}->{std}, "\t";
		print QCSUMMARY $sampleStats->{$key}->{nonbarcodeprimer}->{sum}, "\t";

		print QCSUMMARY $sampleStats->{$key}->{qc}->{count}, "\t";
		$percent = sprintf("%.2f", 100*$sampleStats->{$key}->{qc}->{count}/$sampleStats->{$key}->{raw}->{read1}->{count});
		print QCSUMMARY "$percent\t";
		print QCSUMMARY $sampleStats->{$key}->{qc}->{min}, "\t";
		print QCSUMMARY $sampleStats->{$key}->{qc}->{max}, "\t";
		print QCSUMMARY $sampleStats->{$key}->{qc}->{mean}, "\t";
		print QCSUMMARY $sampleStats->{$key}->{qc}->{median}, "\t";
		print QCSUMMARY $sampleStats->{$key}->{qc}->{std}, "\t";
		print QCSUMMARY $sampleStats->{$key}->{qc}->{sum}, "\t";
		print QCSUMMARY "\n";
	    }
	}
    }
close(QCSUMMARY);

    print $RESULTS <<RESULTS_QC_QC;
    <div id="tabs-nohdr-2">
	<p class="thick"><a href="QC/$analysisname.qc.summary.txt">Download reads quality control summary in text format</a></p>
	<h4>QC Reads Statistical Summary</h4>
	<hr/>
	<div style="height: 20px; "></div>
	<table id="qcstats" class="cell-border" width="85%" cellspacing="0">
	<thead>
                     <tr>
                        <th>Name</th>
                         <th>Count</th>
			 <th>Percent\%</th>
                        <th>Min</th>
                        <th>Max</th>
                        <th>Mean</th>
                        <th>Median</th>
			<th>STD</th>
                        <th>Total bases</th>
			</tr>
                  </thead>
                  <tbody>
RESULTS_QC_QC

my $poolStat = $stats->{"POOL"}->{qc};
    my $percentPool = sprintf("%.2f", 100*$poolStat->{count}/$totalRawCount_read1);
    print $RESULTS "<tr class=\"tr_result\">\n";
    print $RESULTS "<td> Pooled Sample</td>\n";
    print $RESULTS "<td>", format_number($poolStat->{count}),"<\/td>\n";
    print $RESULTS "<td>$percentPool<\/td>\n";
    print $RESULTS "<td>", format_number($poolStat->{min}),"<\/td>\n";
    print $RESULTS "<td>", format_number($poolStat->{max}),"<\/td>\n";
    print $RESULTS "<td>", format_number($poolStat->{mean}),"<\/td>\n";
    print $RESULTS "<td>", format_number($poolStat->{median}),"<\/td>\n";
    print $RESULTS "<td>", format_number($poolStat->{std}),"<\/td>\n";
    print $RESULTS "<td>", format_number($poolStat->{sum}),"<\/td>\n";
    print $RESULTS "<\/tr>\n";
foreach my $name (%$stats){
    if(exists $stats->{$name} && $name ne "POOL"){
	my $qcStat = $stats->{$name}->{qc};
	my $read1Stat = $sampleStats->{$name}->{raw}->{read1};
	my $percent = sprintf("%.2f", 100*$qcStat->{count}/$read1Stat->{count});
	print $RESULTS "<tr class=\"tr_result\">\n";
	print $RESULTS "<td>$name</td>\n";
	print $RESULTS "<td>", format_number($qcStat->{count}),"<\/td>\n";
	print $RESULTS "<td>$percent<\/td>\n";
	print $RESULTS "<td>", format_number($qcStat->{min}),"<\/td>\n";
	print $RESULTS "<td>", format_number($qcStat->{max}),"<\/td>\n";
	print $RESULTS "<td>", format_number($qcStat->{mean}),"<\/td>\n";
	print $RESULTS "<td>", format_number($qcStat->{median}),"<\/td>\n";
	print $RESULTS "<td>", format_number($qcStat->{std}),"<\/td>\n";
	print $RESULTS "<td>", format_number($qcStat->{sum}),"<\/td>\n";
	print $RESULTS "<\/tr>\n";

    }
}

    print $RESULTS "<\/tbody>\n";
    print $RESULTS "<\/table>\n";

    #if($seqtype eq "single"){

	print $RESULTS <<RESULTS_QC_NONTECHSEQ;
	<h4>Reads Statistical Summary After Stripping Off The Technical Sequences</h4>
	    <hr/>
	    <div style="height: 20px; "></div>
	    <table id="trimstats" class="cell-border" width="85%" cellspacing="0">
	    <thead>

	    <tr>
	    <th>Name</th>
	    <th>Count</th>
	    <th>Percent\%</th>
	    <th>Min</th>
	    <th>Max</th>
	    <th>Mean</th>
	    <th>Median</th>
	    <th>STD</th>
	    <th>Total bases</th>
	    </tr>

	    </thead>
	    <tbody>
RESULTS_QC_NONTECHSEQ
foreach my $name (%$stats){
    if(exists $stats->{$name} && $name ne "POOL"){
	#my $cmd = "$bin/visStats.pl -i stats/non_barcode_primer_reads/$name.strip.stats > stats/non_barcode_primer_reads/$name.strip.stats.html";
	#runcmds($cpus, $cmd);
	my $nonbpStat = $stats->{$name}->{nonbarcodeprimer};
	my $read1Stat = $sampleStats->{$name}->{raw}->{read1};
	my $percent_passprimer = sprintf("%.2f",100*$nonbpStat->{count}/$read1Stat->{count});
	print $RESULTS "<tr class=\"tr_result\">\n";
	print $RESULTS "<td>$name</td>\n";
	print $RESULTS "<td>", format_number($nonbpStat->{count}),"<\/td>\n";
	print $RESULTS "<td>$percent_passprimer<\/td>\n";
	print $RESULTS "<td>", format_number($nonbpStat->{min}),"<\/td>\n";
	print $RESULTS "<td>", format_number($nonbpStat->{max}),"<\/td>\n";
	print $RESULTS "<td>", format_number($nonbpStat->{mean}),"<\/td>\n";
	print $RESULTS "<td>", format_number($nonbpStat->{median}),"<\/td>\n";
	print $RESULTS "<td>", format_number($nonbpStat->{std}),"<\/td>\n";
	print $RESULTS "<td>", format_number($nonbpStat->{sum}),"<\/td>\n";
	print $RESULTS "<\/tr>\n";

    }
}
	print $RESULTS "<\/tbody>\n";
	print $RESULTS "<\/table>\n";

    #}

    if($seqtype eq "paired"){

	print $RESULTS <<RESULTS_QC_NONTECHSEQ;
	<h4>Reads Statistical Summary After Merging</h4>
	    <hr/>
	    <div style="height: 20px; "></div>
        <table id="mergestats" class="cell-border" width="85%" cellspacing="0">
        <thead>
	<tr >
	<th>Name</th>
	<th>Count</th>
	<th>Percent\%</th>
	<th>Min</th>
	<th>Max</th>
	<th>Mean</th>
	<th>Median</th>
	<th>STD</th>
	<th>Total bases</th>


	</thead>
	<tbody>
RESULTS_QC_NONTECHSEQ
foreach my $name (%$stats){
    if(exists $stats->{$name} && $name ne "POOL"){
	#my $cmd = "$bin/visStats.pl -i stats/non_barcode_primer_reads/$name.strip.stats > stats/non_barcode_primer_reads/$name.strip.stats.html";
	#runcmds($cpus, $cmd);
	#$cmd = "$bin/visStats.pl -i stats/merged/$name.merged.stats > stats/merged/$name.merged.stats.html";
	#runcmds($cpus, $cmd);

	my $nonbpStat = $stats->{$name}->{nonbarcodeprimer};
	my $mergedStat =  $stats->{$name}->{merged};
	my $read1Stat = $sampleStats->{$name}->{raw}->{read1};
	my $percent_merged = sprintf("%.2f",100*$mergedStat->{count}/$read1Stat->{count});
	my $percent_passprimer = sprintf("%.2f",100*$nonbpStat->{count}/$read1Stat->{count});
	print $RESULTS "<tr class=\"tr_result\">\n";
	print $RESULTS "<td>$name</td>\n";

	print $RESULTS "<td>", format_number($mergedStat->{count}),"<\/td>\n";
	print $RESULTS "<td>$percent_merged<\/td>\n";
	print $RESULTS "<td>", format_number($mergedStat->{min}),"<\/td>\n";
	 print $RESULTS "<td>", format_number($mergedStat->{max}),"<\/td>\n";
	print $RESULTS "<td>", format_number($mergedStat->{mean}),"<\/td>\n";
	print $RESULTS "<td>", format_number($mergedStat->{median}),"<\/td>\n";
	print $RESULTS "<td>", format_number($mergedStat->{std}),"<\/td>\n";
	print $RESULTS "<td>", format_number($mergedStat->{sum}),"<\/td>\n";

	print $RESULTS "<\/tr>\n";

    }
}
	print $RESULTS "<\/tbody>\n";
	print $RESULTS "<\/table>\n";

    }



if($seqtype eq "single"){
    print $RESULTS <<RESULTS_QC_RAW;

    <h4>Raw Reads Statistical Summary</h4>
	<hr/>
	<div style="height: 20px; "></div>
	<table id="rawstats" class="cell-border" width="85%" cellspacing="0">
	<thead>
	<tr>
	<th>Name</th>
	<th>Count</th>
	<th>Min</th>
	<th>Max</th>
	<th>Mean</th>
	<th>Median</th>
	<th>STD</th>
	<th>Total bases</th>
	</tr>
	</thead>
	<tbody>
RESULTS_QC_RAW

foreach my $name (%$stats){
    if(exists $stats->{$name} && $name ne "POOL"){
	#my $cmd = "$bin/visStats.pl -i stats/rawreads/$name.raw.stats > stats/rawreads/$name.raw.stats.html";
	#runcmds($cpus, $cmd);
	my $rawStat = $stats->{$name}->{raw}->{read1};

	print $RESULTS "<tr class=\"tr_result\">\n";
	print $RESULTS "<td >$name</td>\n";
	print $RESULTS "<td>", format_number($rawStat->{count}),"<\/td>\n";
	print $RESULTS "<td>", format_number($rawStat->{min}),"<\/td>\n";
	print $RESULTS "<td>", format_number($rawStat->{max}),"<\/td>\n";
	print $RESULTS "<td>", format_number($rawStat->{mean}),"<\/td>\n";
	print $RESULTS "<td>", format_number($rawStat->{median}),"<\/td>\n";
	print $RESULTS "<td>", format_number($rawStat->{std}),"<\/td>\n";
	print $RESULTS "<td>", format_number($rawStat->{sum}),"<\/td>\n";
	print $RESULTS "<\/tr>\n";

    }
}

    print $RESULTS "<\/tbody>\n";
    print $RESULTS "<\/table>\n";
    print $RESULTS "<\/div>\n";

}

elsif($seqtype eq "paired"){
    print $RESULTS <<RESULTS_QC_RAW;

    <h4>Raw Reads Statistical Summary</h4>
        <hr/>
        <div style="height: 20px; "></div>
        <table id="rawstats" class="cell-border" width="85%" cellspacing="0">
        <thead>
	<tr>
	<th rowspan=2>Name</th>
	<th colspan=2>R1</th>
	<th colspan=2>R2</th>
	</tr>
        <tr>
	<th>Count</th>
        <th>Total bases</th>
	<th>Count</th>
        <th>Total bases</th>
        </tr>
        </thead>
        <tbody>
RESULTS_QC_RAW

foreach my $name (%$stats){
    if(exists $stats->{$name} && $name ne "POOL"){
        #my $cmd = "$bin/visStats.pl -i stats/rawreads/$name.R1.raw.stats > stats/rawreads/$name.R1.raw.stats.html;";
	#$cmd .= "$bin/visStats.pl -i stats/rawreads/$name.R2.raw.stats > stats/rawreads/$name.R2.raw.stats.html;";
	#runcmds($cpus, $cmd);
        my $rawStat = $stats->{$name}->{raw};

        print $RESULTS "<tr class=\"tr_result\">\n";
	print $RESULTS "<td >$name</td>\n";
	print $RESULTS "<td>", format_number($rawStat->{read1}->{count}),"<\/td>\n";
        print $RESULTS "<td>", format_number($rawStat->{read1}->{sum}),"<\/td>\n";
	print $RESULTS "<td>", format_number($rawStat->{read2}->{count}),"<\/td>\n";
        print $RESULTS "<td>", format_number($rawStat->{read2}->{sum}),"<\/td>\n";
	print $RESULTS "<\/tr>\n";

    }
}

    print $RESULTS "<\/tbody>\n";
    print $RESULTS "<\/table>\n";
    print $RESULTS "<\/div>\n";



}
}

sub print_web_otus{
    my ($RESULTS, $analysis) = @_;

    my $cmd = "$bin/buildTree.pl -idir OTU_and_taxonomy -i $analysis.taxonomy.norank.summary -odir html";
    runcmds($cpus, $cmd);


    print $RESULTS <<RESULTS_OTUS;
    <div id="tabs-nohdr-3">
	<p class="thick">» OTUs</p>
    <ul id="lfile">
    <li><a href="OTU_and_taxonomy/$analysis.OTUs.fasta" title="Representative sequence for each OTU:otuid|sequence count" target="vis">OTU representative sequences</a></li>

    </ul>
RESULTS_OTUS

    if($geneType eq "rRNA"){
	print $RESULTS <<RESULTS_TAXON;
	<p class="thick">» Taxonomy annotation</p>
    <ul id="lfile">

    <li>OTU taxonomic assignment <a href="OTU_and_taxonomy/$analysis.OTUs.nr_v132.mwang.taxonomy" title="The consensus taxonomy of each OTU" target="vis">with</a> and <a href="OTU_and_taxonomy/$analysis.OTUs.nr_v132.mwang.norank.taxonomy" title="The consensus taxonomy of each OTU without rank information" target="vis">without</a> rank</li>
    <li><a href="OTU_and_taxonomy/$analysis.OTU-table.taxonomy" title="OTU table with taxonomic classifiion" target="vis">OTU table & taxonomic assignment</a></li>

    <li>Taxonomic distribution table <a href="OTU_and_taxonomy/$analysis.taxonomy.summary" title="Taxonomic summary table of the analysis" target="vis">with</a> and <a href="OTU_and_taxonomy/$analysis.taxonomy.norank.summary" title="Taxonomic summary table of the analysis" target="vis">without</a> rank
	<a href="html/$analysis.taxonomy.norank.summary.html" target="vis"><img src="images/chart.gif" height="20"></a></li>
    </ul>

RESULTS_TAXON
}
print $RESULTS "</div>\n";

}

sub print_web_alpha_diversity{

    my($RESULTS, $analysis, $samplenames) = @_;

    my $cmd = "$bin/visTSV.pl -i alpha_and_beta_diversity/$analysis.groups.rarefaction -t $projDir/template/rarefaction.html> html/$analysis.groups.rarefaction.html;";
    $cmd .= "$bin/visTSV.pl -i alpha_and_beta_diversity/$analysis.$dissimilarity.subsample.groups.rarefaction -t $projDir/template/rarefaction.html> html/$analysis.$dissimilarity.subsample.groups.rarefaction.html;";
    runcmds($cpus, $cmd);
    $cmd = "$bin/visTSV.pl -i alpha_and_beta_diversity/$analysis.groups.summary -t $projDir/template/alpha_diversity_indexes.html> html/$analysis.groups.summary.html;";
    $cmd .= "$bin/visTSV.pl -i alpha_and_beta_diversity/$analysis.$dissimilarity.subsample.groups.summary -t $projDir/template/alpha_diversity_indexes.html> html/$analysis.$dissimilarity.subsample.groups.summary.html";
    runcmds($cpus, $cmd);


    $cmd = "$bin/visPCoA.pl -a alpha_and_beta_diversity/$analysis.braycurtis.$dissimilarity.lt.pcoa.axes -l alpha_and_beta_diversity/$analysis.braycurtis.$dissimilarity.lt.pcoa.loadings -t $projDir/template/pcoa.html> html/$analysis.braycurtis.$dissimilarity.lt.pcoa.html;";
    $cmd .= "$bin/visPCoA.pl -a alpha_and_beta_diversity/$analysis.jclass.$dissimilarity.lt.pcoa.axes -l alpha_and_beta_diversity/$analysis.jclass.$dissimilarity.lt.pcoa.loadings -t $projDir/template/pcoa.html> html/$analysis.jclass.$dissimilarity.lt.pcoa.html;";
    $cmd .= "$bin/visPCoA.pl -a alpha_and_beta_diversity/$analysis.jest.$dissimilarity.lt.pcoa.axes -l alpha_and_beta_diversity/$analysis.jest.$dissimilarity.lt.pcoa.loadings -t $projDir/template/pcoa.html> html/$analysis.jest.$dissimilarity.lt.pcoa.html;";
    $cmd .= "$bin/visPCoA.pl -a alpha_and_beta_diversity/$analysis.thetayc.$dissimilarity.lt.pcoa.axes -l alpha_and_beta_diversity/$analysis.thetayc.$dissimilarity.lt.pcoa.loadings -t $projDir/template/pcoa.html> html/$analysis.thetayc.$dissimilarity.lt.pcoa.html;";

    $cmd .= "$bin/visPCoA.pl -a alpha_and_beta_diversity/$analysis.$dissimilarity.subsample.braycurtis.$dissimilarity.lt.pcoa.axes -l alpha_and_beta_diversity/$analysis.$dissimilarity.subsample.braycurtis.$dissimilarity.lt.pcoa.loadings -t $projDir/template/pcoa.html> html/$analysis.$dissimilarity.subsample.braycurtis.$dissimilarity.lt.pcoa.html;";
    $cmd .= "$bin/visPCoA.pl -a alpha_and_beta_diversity/$analysis.$dissimilarity.subsample.jclass.$dissimilarity.lt.pcoa.axes -l alpha_and_beta_diversity/$analysis.$dissimilarity.subsample.jclass.$dissimilarity.lt.pcoa.loadings -t $projDir/template/pcoa.html> html/$analysis.$dissimilarity.subsample.jclass.$dissimilarity.lt.pcoa.html;";
    $cmd .= "$bin/visPCoA.pl -a alpha_and_beta_diversity/$analysis.$dissimilarity.subsample.jest.$dissimilarity.lt.pcoa.axes -l alpha_and_beta_diversity/$analysis.$dissimilarity.subsample.jest.$dissimilarity.lt.pcoa.loadings -t $projDir/template/pcoa.html> html/$analysis.$dissimilarity.subsample.jest.$dissimilarity.lt.pcoa.html;";
    $cmd .= "$bin/visPCoA.pl -a alpha_and_beta_diversity/$analysis.$dissimilarity.subsample.thetayc.$dissimilarity.lt.pcoa.axes -l alpha_and_beta_diversity/$analysis.$dissimilarity.subsample.thetayc.$dissimilarity.lt.pcoa.loadings -t $projDir/template/pcoa.html> html/$analysis.$dissimilarity.subsample.thetayc.$dissimilarity.lt.pcoa.html;";
    runcmds($cpus, $cmd);


      $cmd = "$bin/visNMDs.pl -s alpha_and_beta_diversity/$analysis.braycurtis.$dissimilarity.lt.nmds.stress -a alpha_and_beta_diversity/$analysis.braycurtis.$dissimilarity.lt.nmds.axes -t $projDir/template/nmds.html> html/$analysis.braycurtis.$dissimilarity.lt.nmds.html;";
    $cmd .= "$bin/visNMDs.pl -s alpha_and_beta_diversity/$analysis.thetayc.$dissimilarity.lt.nmds.stress -a alpha_and_beta_diversity/$analysis.thetayc.$dissimilarity.lt.nmds.axes -t $projDir/template/nmds.html> html/$analysis.thetayc.$dissimilarity.lt.nmds.html;";
    $cmd .= "$bin/visNMDs.pl -s alpha_and_beta_diversity/$analysis.jclass.$dissimilarity.lt.nmds.stress -a alpha_and_beta_diversity/$analysis.jclass.$dissimilarity.lt.nmds.axes -t $projDir/template/nmds.html> html/$analysis.jclass.$dissimilarity.lt.nmds.html;";
    $cmd .= "$bin/visNMDs.pl -s alpha_and_beta_diversity/$analysis.jest.$dissimilarity.lt.nmds.stress -a alpha_and_beta_diversity/$analysis.jest.$dissimilarity.lt.nmds.axes -t $projDir/template/nmds.html> html/$analysis.jest.$dissimilarity.lt.nmds.html;";

    $cmd .=  "$bin/visNMDs.pl -s alpha_and_beta_diversity/$analysis.$dissimilarity.subsample.braycurtis.$dissimilarity.lt.nmds.stress -a alpha_and_beta_diversity/$analysis.$dissimilarity.subsample.braycurtis.$dissimilarity.lt.nmds.axes -t $projDir/template/nmds.html> html/$analysis.$dissimilarity.subsample.braycurtis.$dissimilarity.lt.nmds.html;";
    $cmd .= "$bin/visNMDs.pl -s alpha_and_beta_diversity/$analysis.$dissimilarity.subsample.thetayc.$dissimilarity.lt.nmds.stress -a alpha_and_beta_diversity/$analysis.$dissimilarity.subsample.thetayc.$dissimilarity.lt.nmds.axes -t $projDir/template/nmds.html> html/$analysis.$dissimilarity.subsample.thetayc.$dissimilarity.lt.nmds.html;";
    $cmd .= "$bin/visNMDs.pl -s alpha_and_beta_diversity/$analysis.$dissimilarity.subsample.jclass.$dissimilarity.lt.nmds.stress -a alpha_and_beta_diversity/$analysis.$dissimilarity.subsample.jclass.$dissimilarity.lt.nmds.axes -t $projDir/template/nmds.html> html/$analysis.$dissimilarity.subsample.jclass.$dissimilarity.lt.nmds.html;";
    $cmd .= "$bin/visNMDs.pl -s alpha_and_beta_diversity/$analysis.$dissimilarity.subsample.jest.$dissimilarity.lt.nmds.stress -a alpha_and_beta_diversity/$analysis.$dissimilarity.subsample.jest.$dissimilarity.lt.nmds.axes -t $projDir/template/nmds.html> html/$analysis.$dissimilarity.subsample.jest.$dissimilarity.lt.nmds.html;";

    runcmds($cpus, $cmd);


    $cmd = "$bin/visNewickTree.pl -i alpha_and_beta_diversity/$analysis.braycurtis.$dissimilarity.lt.tre > html/$analysis.braycurtis.$dissimilarity.lt.tre.html;";
    $cmd .= "$bin/visNewickTree.pl -i alpha_and_beta_diversity/$analysis.thetayc.$dissimilarity.lt.tre > html/$analysis.thetayc.$dissimilarity.lt.tre.html;";
    $cmd .= "$bin/visNewickTree.pl -i alpha_and_beta_diversity/$analysis.jclass.$dissimilarity.lt.tre > html/$analysis.jclass.$dissimilarity.lt.tre.html;";
    $cmd .= "$bin/visNewickTree.pl -i alpha_and_beta_diversity/$analysis.jest.$dissimilarity.lt.tre > html/$analysis.jest.$dissimilarity.lt.tre.html;";

    $cmd .= "$bin/visNewickTree.pl -i alpha_and_beta_diversity/$analysis.$dissimilarity.subsample.braycurtis.$dissimilarity.lt.tre > html/$analysis.$dissimilarity.subsample.braycurtis.$dissimilarity.lt.tre.html;";
    $cmd .= "$bin/visNewickTree.pl -i alpha_and_beta_diversity/$analysis.$dissimilarity.subsample.thetayc.$dissimilarity.lt.tre > html/$analysis.$dissimilarity.subsample.thetayc.$dissimilarity.lt.tre.html;";
    $cmd .= "$bin/visNewickTree.pl -i alpha_and_beta_diversity/$analysis.$dissimilarity.subsample.jclass.$dissimilarity.lt.tre > html/$analysis.$dissimilarity.subsample.jclass.$dissimilarity.lt.tre.html;";
    $cmd .= "$bin/visNewickTree.pl -i alpha_and_beta_diversity/$analysis.$dissimilarity.subsample.jest.$dissimilarity.lt.tre > html/$analysis.$dissimilarity.subsample.jest.$dissimilarity.lt.tre.html;";

    runcmds($cpus, $cmd);

    $cmd = "$bin/make_gexf.pl -a alpha_and_beta_diversity/network.$dissimilarity.spearman.otu.corr -s alpha_and_beta_diversity/network.shared -c 0.6 -p 0.05 >  alpha_and_beta_diversity/network.$dissimilarity.spearman.otu.corr.p05.gexf;";
    $cmd .= "$bin/make_gexf.pl -a alpha_and_beta_diversity/network.$dissimilarity.spearman.otu.corr -s alpha_and_beta_diversity/network.shared -c 0.6 -p 0.01 >  alpha_and_beta_diversity/network.$dissimilarity.spearman.otu.corr.p01.gexf;";
    $cmd .= "$bin/make_gexf.pl -a alpha_and_beta_diversity/network.$dissimilarity.kendall.otu.corr -s alpha_and_beta_diversity/network.shared -c 0.6 -p 0.05 >  alpha_and_beta_diversity/network.$dissimilarity.kendall.otu.corr.p05.gexf;";
    $cmd .= "$bin/make_gexf.pl -a alpha_and_beta_diversity/network.$dissimilarity.kendall.otu.corr -s alpha_and_beta_diversity/network.shared -c 0.6 -p 0.01 >  alpha_and_beta_diversity/network.$dissimilarity.kendall.otu.corr.p01.gexf;";
    $cmd .= "$bin/make_gexf.pl -a alpha_and_beta_diversity/network.$dissimilarity.pearson.otu.corr -s alpha_and_beta_diversity/network.shared -c 0.6 -p 0.05 >  alpha_and_beta_diversity/network.$dissimilarity.pearson.otu.corr.p05.gexf;";
    $cmd .= "$bin/make_gexf.pl -a alpha_and_beta_diversity/network.$dissimilarity.pearson.otu.corr -s alpha_and_beta_diversity/network.shared -c 0.6 -p 0.01 >  alpha_and_beta_diversity/network.$dissimilarity.pearson.otu.corr.p01.gexf;";


    $cmd .= "$bin/make_gexf.pl -a alpha_and_beta_diversity/network.$dissimilarity.subsample.$dissimilarity.spearman.otu.corr -s alpha_and_beta_diversity/network.$dissimilarity.subsample.shared -c 0.6 -p 0.05 >  alpha_and_beta_diversity/network.$dissimilarity.subsample.$dissimilarity.spearman.otu.corr.p05.gexf;";
    $cmd .= "$bin/make_gexf.pl -a alpha_and_beta_diversity/network.$dissimilarity.subsample.$dissimilarity.spearman.otu.corr -s alpha_and_beta_diversity/network.$dissimilarity.subsample.shared -c 0.6 -p 0.01 >  alpha_and_beta_diversity/network.$dissimilarity.subsample.$dissimilarity.spearman.otu.corr.p01.gexf;";
    $cmd .= "$bin/make_gexf.pl -a alpha_and_beta_diversity/network.$dissimilarity.subsample.$dissimilarity.kendall.otu.corr -s alpha_and_beta_diversity/network.$dissimilarity.subsample.shared -c 0.6 -p 0.05 >  alpha_and_beta_diversity/network.$dissimilarity.subsample.$dissimilarity.kendall.otu.corr.p05.gexf;";
    $cmd .= "$bin/make_gexf.pl -a alpha_and_beta_diversity/network.$dissimilarity.subsample.$dissimilarity.kendall.otu.corr -s alpha_and_beta_diversity/network.$dissimilarity.subsample.shared -c 0.6 -p 0.01 >  alpha_and_beta_diversity/network.$dissimilarity.subsample.$dissimilarity.kendall.otu.corr.p01.gexf;";
    $cmd .= "$bin/make_gexf.pl -a alpha_and_beta_diversity/network.$dissimilarity.subsample.$dissimilarity.pearson.otu.corr -s alpha_and_beta_diversity/network.$dissimilarity.subsample.shared -c 0.6 -p 0.05 >  alpha_and_beta_diversity/network.$dissimilarity.subsample.$dissimilarity.pearson.otu.corr.p05.gexf;";
    $cmd .= "$bin/make_gexf.pl -a alpha_and_beta_diversity/network.$dissimilarity.subsample.$dissimilarity.pearson.otu.corr -s alpha_and_beta_diversity/network.$dissimilarity.subsample.shared -c 0.6 -p 0.01 >  alpha_and_beta_diversity/network.$dissimilarity.subsample.$dissimilarity.pearson.otu.corr.p01.gexf;";

    runcmds($cpus, $cmd);

    print $RESULTS <<RESULTS_ALPHA;



    <div id="tabs-nohdr-5">
	<div id="wrap">
	<div id="left_col">
	<h3> Without subsampling or rarefying</h3>
	</div>
	<div id="right_col">
	<h3> With subsampling or rarefying &nbsp;&nbsp; <input TYPE="button" id="rarefy" class="btn" name="HELP" value="HELP" onclick="help('rarefyModal', 'rarefy', 'closeRarefy');"></h3>
	</div>

	<div id="rarefyModal" class="modal">
	<!-- Modal content -->
	<div class="modal-content">
	<span class="close" id="closeRarefy">&times;</span>
	<p>There\'s a subtle, but very important, difference between subsampling and rarefaction. Subsampling will draw a desired number of sequences from each sample - once. This is what Mothur do for things like metastats or machine learning models. For alpha and beta diversity, Mothur rarefy the data. This is subsampling 1000 times and averaging the 1000 subsamples. With the large number of randomizations, you average out the effects of getting lucky/unlucky.</p>
	</div>
	</div>

	</div>
	<p class="thick">»  OTU and &alpha;-diversity measurements</p>

    <div id="wrap">
    <div id="left_col">
    <ul id="lfile">
    <li><a href="alpha_and_beta_diversity/$analysis.list" title="Contains data indicating the sequences that cluster together within an OTU" target="vis">Mothur format List file</a><a href="https:\//www.mothur.org/wiki/List_file" target="vis"><img src="images/help.jpg" height="20" class="help"></li>
    <li><a href="alpha_and_beta_diversity/$analysis.groups" title="Sequence-to-sample mapping (useful for mothur's beta diversity analysis)" target="vis">Mothur format group file</a><a href="http://www.mothur.org/wiki/Group_file" target="vis"><img src="images/help.jpg" height="20" class="help"></li>
    <li><a href="alpha_and_beta_diversity/$analysis.shared" title="The data in a shared file represents the number of times that an OTU is observed in multiple samples." target="vis" >Mothur format shared file</a><a href="http://www.mothur.org/wiki/Shared_file" target="vis"><img src="images/help.jpg" height="20" class="help"></li>
    <li><a href="alpha_and_beta_diversity/$analysis.relabund" title="relabund file represent relative abundance of an OTU over the samples" target="vis" >OTU relabund file</a><a href="http://www.mothur.org/wiki/Shared_file" target="vis"><img src="images/help.jpg" height="20" class="help"></a></li>
    <li><a href="alpha_and_beta_diversity/$analysis.groups.rarefaction" title="Rarefaction curve describing the number of OTUs observed as a function of sampling effort" target="vis">Rarefaction curve</a><a href="http://www.mothur.org/wiki/Rarefaction" target="vis"><img src="images/help.jpg" height="20" class="help"></a><a href="html/$analysis.groups.rarefaction.html" target="vis"><img src="images/chart.gif" height="20" class="chart"></a></li>
    <li class=""><a href="alpha_and_beta_diversity/$analysis.groups.summary" title="Alpha diversity indices using different calculator" target="vis">Alpha diversity index</a><a href="http://www.mothur.org/wiki/Summary.single" target="vis"><img src="images/help.jpg" height="20" class="help"></a><a href="html/$analysis.groups.summary.html" target="vis"><img src="images/chart.gif" height="20" class="chart"></a></li>
    </ul>
    </div>
    <div id="right_col">
    <ul id="lfile">
    <li><a href="alpha_and_beta_diversity/$analysis.$dissimilarity.subsample.list" title="Contains data indicating the sequences that cluster together within an OTU" target="vis">Mothur format List file</a></li>
    <li><a href="alpha_and_beta_diversity/$analysis.$dissimilarity.subsample.groups" title="Sequence-to-sample mapping (useful for mothur\'s beta diversity analysis)" target="vis">Mothur format group file</a></li>
    <li><a href="alpha_and_beta_diversity/$analysis.$dissimilarity.subsample.shared" title="The data in a shared file represents the number of times that an OTU is observed in multiple samples." target="vis" >Mothur format shared file</a></li>
    <li><a href="alpha_and_beta_diversity/$analysis.$dissimilarity.subsample.relabund" title="relabund file represent relative abundance of an OTU over the samples" target="vis" >OTU relabund file</a></li>
    <li><a href="alpha_and_beta_diversity/$analysis.$dissimilarity.subsample.groups.rarefaction" title="Rarefaction curve describing the number of OTUs observed as a function of sampling effort" target="vis">Rarefaction curve</a><a href="html/$analysis.$dissimilarity.subsample.groups.rarefaction.html" target="vis"><img src="images/chart.gif" height="20" class="chart"></a></li>
    <li class=""><a href="alpha_and_beta_diversity/$analysis.$dissimilarity.subsample.groups.summary" title="Alpha diversity indices using different calculator" target="vis">Alpha diversity index</a><a href="html/$analysis.$dissimilarity.subsample.groups.summary.html" target="vis"><img src="images/chart.gif" height="20" class="chart"></a></li>
    </ul>
    </div>
    </div>

 <p class="thick">» &beta;-diversity measurements: OTU association &nbsp;&nbsp; <input TYPE="button" id="otuasso" class="btn" value="HELP" onClick="help('otuassoModal', 'otuasso', 'closeOTUasso');"></p>
     <div id="otuassoModal" class="modal">
        <!-- Modal content -->
        <div class="modal-content">
        <span class="close" id="closeOTUasso">&times;</span>
        <p>1) generate matrix: use mothur generated shared file and filtered out otus whose abundance in the sample were less than 0.1% or the otus showing up only in less than two samples. 2)use filered shared file to create a matrix of correlation values indicating the association between different OTUs. Just remember that association is not the same as causation :). When we generate the gexf file,the correlation cofficient must be greater than 0.6 or less than -0.6. FDR q-values based on p-values by multiple testing correction (Benjamini and Hochberg, 1995).</p>
        </div>
        </div>

    <div id="wrap">
    <div id="left_col">
    <ul id="lfile">
    <li>Spearman <a href="alpha_and_beta_diversity/network.$dissimilarity.spearman.qvalue.otu.corr" target="vis">OTU correlation coefficient</a></li>
    <li>Spearman gexf files at <a href="alpha_and_beta_diversity/network.$dissimilarity.spearman.otu.corr.p05.gexf" target="vis">0.05</a> and <a href="alpha_and_beta_diversity/network.$dissimilarity.spearman.otu.corr.p01.gexf" target="vis">0.01</a> p-value cutoff</a></li>
    <li>Pearson <a href="alpha_and_beta_diversity/network.$dissimilarity.pearson.qvalue.otu.corr" target="vis">OTU correlation coefficient</a></li>
    <li>Pearson gexf files at <a href="alpha_and_beta_diversity/network.$dissimilarity.pearson.otu.corr.p05.gexf" target="vis">0.05</a> and <a href="alpha_and_beta_diversity/network.$dissimilarity.pearson.otu.corr.p01.gexf" target="vis">0.01</a> p-value cutoff</a></li>

    <li>Kendall <a href="alpha_and_beta_diversity/network.$dissimilarity.kendall.qvalue.otu.corr" target="vis">OTU correlation coefficient</a></li>
    <li>Kendall gexf files at <a href="alpha_and_beta_diversity/network.$dissimilarity.kendall.otu.corr.p05.gexf" target="vis">0.05</a> and <a href="alpha_and_beta_diversity/network.$dissimilarity.kendall.otu.corr.p01.gexf" target="vis">0.01</a> p-value cutoff</a></li>
    </ul>
    </div>
    <div id="right_col">
     <ul id="lfile">
    <li>Spearman <a href="alpha_and_beta_diversity/network.$dissimilarity.subsample.$dissimilarity.spearman.qvalue.otu.corr" target="vis">OTU correlation coefficient</a></li>
    <li>Spearman gexf files at <a href="alpha_and_beta_diversity/network.$dissimilarity.subsample.$dissimilarity.spearman.otu.corr.p05.gexf" target="vis">0.05</a> and <a href="alpha_and_beta_diversity/network.$dissimilarity.subsample.$dissimilarity.spearman.otu.corr.p01.gexf" target="vis">0.01</a> p-value cutoff</a></li>
    <li>Pearson <a href="alpha_and_beta_diversity/network.$dissimilarity.subsample.$dissimilarity.pearson.qvalue.otu.corr" target="vis">OTU correlation coefficient</a></li>
    <li>Pearson gexf files at <a href="alpha_and_beta_diversity/network.$dissimilarity.subsample.$dissimilarity.pearson.otu.corr.p05.gexf" target="vis">0.05</a> and <a href="alpha_and_beta_diversity/network.$dissimilarity.subsample.$dissimilarity.pearson.otu.corr.p01.gexf" target="vis">0.01</a> p-value cutoff</a></li>
    <li>Kendall <a href="alpha_and_beta_diversity/network.$dissimilarity.subsample.$dissimilarity.kendall.qvalue.otu.corr" target="vis">OTU correlation coefficient</a></li>
    <li>Kendall gexf files at <a href="alpha_and_beta_diversity/network.$dissimilarity.subsample.$dissimilarity.kendall.otu.corr.p05.gexf" target="vis">0.05</a> and <a href="alpha_and_beta_diversity/network.$dissimilarity.subsample.$dissimilarity.kendall.otu.corr.p01.gexf" target="vis">0.01</a> p-value cutoff</a></li>
    </ul>

    </div>
    </div>

    <p class="thick">» &beta;-diversity measurements: PCoA analysis &nbsp;&nbsp; <input TYPE="button" class="btn" id="pcoa" value="HELP" onClick="help('pcoaModal', 'pcoa', 'closePcoa');"></p>

    <div id="pcoaModal" class="modal">
    <!-- Modal content -->
    <div class="modal-content">
    <span class="close" id="closePcoa">&times;</span>
    <p>Principal Coordinates (PCoA) uses an eigenvector-based approach to represent multidimensional data in as few dimesnsions as possible.loadings file will tell you what fraction of the total variance in the data are represented by each of the axes.The axes file contains the plotting coordinates for graphing the results</p>
    </div>
    </div>

    <div id="wrap">
    <div id="left_col">
    <ul id="lfile">
    <li>Jclass: <a href="alpha_and_beta_diversity/$analysis.jclass.$dissimilarity.lt.pcoa.loadings" target="vis" title="Principal Coordinates (PCoA) uses an eigenvector-based approach to represent multidimensional data in as few dimesnsions as possible. Loadings file will tell you what fraction of the total variance in the data are represented by each of the axes.">pcoa.loadings</a> & <a href="alpha_and_beta_diversity/$analysis.jclass.$dissimilarity.lt.pcoa.axes"  target="vis">pcoa.axes </a> &nbsp;&nbsp;<a href="html/$analysis.jclass.$dissimilarity.lt.pcoa.html" target="vis"><img src="images/chart.gif" height="20" class="chart"></a>
    <li>Jest: <a href="alpha_and_beta_diversity/$analysis.jest.$dissimilarity.lt.pcoa.loadings" target="vis" title="Principal Coordinates (PCoA) uses an eigenvector-based approach to represent multidimensional data in as few dimesnsions as possible. Loadings file will tell you what fraction of the total variance in the data are represented by each of the axes.">pcoa.loadings</a> & <a href="alpha_and_beta_diversity/$analysis.jest.$dissimilarity.lt.pcoa.axes" target="vis">pcoa.axes</a> &nbsp;&nbsp;<a href="html/$analysis.jest.$dissimilarity.lt.pcoa.html" target="vis"><img src="images/chart.gif" height="20" class="chart"></a>
    <li>Thetayc: <a href="alpha_and_beta_diversity/$analysis.thetayc.$dissimilarity.lt.pcoa.loadings" target="vis" title="Principal Coordinates (PCoA) uses an eigenvector-based approach to represent multidimensional data in as few dimesnsions as possible. Loadings file will tell you what fraction of the total variance in the data are represented by each of the axes.">pcoa.loadings</a> & <a href="alpha_and_beta_diversity/$analysis.thetayc.$dissimilarity.lt.pcoa.axes" target="vis">pcoa.axes</a> &nbsp;&nbsp;<a href="html/$analysis.thetayc.$dissimilarity.lt.pcoa.html" target="vis"><img src="images/chart.gif" height="20" class="chart"></a>
    <li>Braycurtis: <a href="alpha_and_beta_diversity/$analysis.braycurtis.$dissimilarity.lt.pcoa.loadings" target="vis" title="Principal Coordinates (PCoA) uses an eigenvector-based approach to represent multidimensional data in as few dimesnsions as possible. Loadings file will tell you what fraction of the total variance in the data are represented by each of the axes.">pcoa.loadings</a> & <a href="alpha_and_beta_diversity/$analysis.braycurtis.$dissimilarity.lt.pcoa.axes" target="vis">pcoa.axes</a> &nbsp;&nbsp;<a href="html/$analysis.braycurtis.$dissimilarity.lt.pcoa.html" target="vis"><img src="images/chart.gif" height="20" class="chart"></a>
    </ul>
    </div>
    <div id="right_col">
    <ul id="lfile">
    <li>Jclass: <a href="alpha_and_beta_diversity/$analysis.$dissimilarity.subsample.jclass.$dissimilarity.lt.pcoa.loadings" target="vis" title="Principal Coordinates (PCoA) uses an eigenvector-based approach to represent multidimensional data in as few dimesnsions as possible. Loadings file will tell you what fraction of the total variance in the data are represented by each of the axes.">pcoa.loadings</a> & <a href="alpha_and_beta_diversity/$analysis.$dissimilarity.subsample.jclass.$dissimilarity.lt.pcoa.axes">pcoa.axes </a> &nbsp;&nbsp;<a href="html/$analysis.$dissimilarity.subsample.jclass.$dissimilarity.lt.pcoa.html" target="vis"><img src="images/chart.gif" height="20" class="chart"></a>
    <li>Jest: <a href="alpha_and_beta_diversity/$analysis.$dissimilarity.subsample.jest.$dissimilarity.lt.pcoa.loadings" target="vis" title="Principal Coordinates (PCoA) uses an eigenvector-based approach to represent multidimensional data in as few dimesnsions as possible. Loadings file will tell you what fraction of the total variance in the data are represented by each of the axes.">pcoa.loadings</a> & <a href="alpha_and_beta_diversity/$analysis.$dissimilarity.subsample.jest.$dissimilarity.lt.pcoa.axes" target="vis">pcoa.axes</a> &nbsp;&nbsp;<a href="html/$analysis.$dissimilarity.subsample.jest.$dissimilarity.lt.pcoa.html" target="vis"><img src="images/chart.gif" height="20" class="chart"></a>
    <li>Thetayc: <a href="alpha_and_beta_diversity/$analysis.$dissimilarity.subsample.thetayc.$dissimilarity.lt.pcoa.loadings" target="vis" title="Principal Coordinates (PCoA) uses an eigenvector-based approach to represent multidimensional data in as few dimesnsions as possible. Loadings file will tell you what fraction of the total variance in the data are represented by each of the axes.">pcoa.loadings</a> & <a href="alpha_and_beta_diversity/$analysis.$dissimilarity.susample.thetayc.$dissimilarity.lt.pcoa.axes" target="vis">pcoa.axes</a> &nbsp;&nbsp;<a href="html/$analysis.$dissimilarity.subsample.thetayc.$dissimilarity.lt.pcoa.html" target="vis"><img src="images/chart.gif" height="20" class="chart"></a>
    <li>Braycurtis: <a href="alpha_and_beta_diversity/$analysis.$dissimilarity.subsample.braycurtis.$dissimilarity.lt.pcoa.loadings" target="vis" title="Principal Coordinates (PCoA) uses an eigenvector-based approach to represent multidimensional data in as few dimesnsions as possible. Loadings file will tell you what fraction of the total variance in the data are represented by each of the axes.">pcoa.loadings</a> & <a href="alpha_and_beta_diversity/$analysis.$dissimilarity.subsample.braycurtis.$dissimilarity.lt.pcoa.axes" target="vis">pcoa.axes</a> &nbsp;&nbsp;<a href="html/$analysis.$dissimilarity.subsample.braycurtis.$dissimilarity.lt.pcoa.html" target="vis"><img src="images/chart.gif" height="20" class="chart"></a>
    </ul>
    </div>
    </div>

    <p class="thick">» &beta;-diversity measurements: Tree analysis &nbsp;&nbsp; <input TYPE="button" class="btn" id="tree" value="HELP" onClick="help('treeModal', 'tree', 'closeTree');"></p>

    <div id="treeModal" class="modal">
    <!-- Modal content -->
    <div class="modal-content">
    <span class="close" id="closeTree">&times;</span>
    <p>generate a newick-formatted tree file that describes the dissimilarity (1-similarity) among multiple samples. Samples are clustered using the UPGMA algorithm using the distance between communities as calculated using any of the calculators describing the similarity in community membership or structure. Dissimilarity is calculated as one minus the similarity.The newick format tree file can be visualized in softwre like Dendroscope.</p>
    </div>
    </div>



    <div id="wrap">
    <div id="left_col">
    <ul id="lfile">
    <li>Jclass <a href="alpha_and_beta_diversity/$analysis.jclass.$dissimilarity.lt.tre" target="vis" >newick format sample relation tree</a> &nbsp;&nbsp;<a href="html/$analysis.jclass.$dissimilarity.lt.tre.html" target="vis"><img src="images/chart.gif" height="20" class="chart"></a>
    <li>Jest <a href="alpha_and_beta_diversity/$analysis.jest.$dissimilarity.lt.tre" target="vis" >newick format sample relation tree</a>  &nbsp;&nbsp;<a href="html/$analysis.jest.$dissimilarity.lt.tre.html" target="vis"><img src="images/chart.gif" height="20" class="chart"></a>
    <li>Thetayc <a href="alpha_and_beta_diversity/$analysis.thetayc.$dissimilarity.lt.tre" target="vis">newick format sample relation tree</a>&nbsp;&nbsp;<a href="html/$analysis.thetayc.$dissimilarity.lt.tre.html" target="vis"><img src="images/chart.gif" height="20" class="chart"></a>
    <li>Braycurtis <a href="alpha_and_beta_diversity/$analysis.braycurtis.$dissimilarity.lt.tre" target="vis">newick format sample relation tree</a>&nbsp;&nbsp;<a href="html/$analysis.braycurtis.$dissimilarity.lt.tre.html" target="vis"><img src="images/chart.gif" height="20" class="chart"></a>
    </ul>
    </div>
    <div id="right_col">
    <ul id="lfile">
    <li>Jclass <a href="alpha_and_beta_diversity/$analysis.$dissimilarity.subsample.jclass.$dissimilarity.lt.tre" target="vis">newick format sample relation tree</a> &nbsp;&nbsp;<a href="html/$analysis.$dissimilarity.subsample.jclass.$dissimilarity.lt.tre.html" target="vis"><img src="images/chart.gif" height="20" class="chart"></a>
    <li>Jest <a href="alpha_and_beta_diversity/$analysis.$dissimilarity.subsample.jest.$dissimilarity.lt.tre" target="vis">newick format sample relation tree</a>  &nbsp;&nbsp;<a href="html/$analysis.$dissimilarity.subsample.jest.$dissimilarity.lt.tre.html" target="vis"><img src="images/chart.gif" height="20" class="chart"></a>
    <li>Thetayc <a href="alpha_and_beta_diversity/$analysis.$dissimilarity.subsample.thetayc.$dissimilarity.lt.tre" target="vis">newick format sample relation tree</a>&nbsp;&nbsp;<a href="html/$analysis.$dissimilarity.subsample.thetayc.$dissimilarity.lt.tre.html" target="vis"><img src="images/chart.gif" height="20" class="chart"></a>
    <li>Braycurtis <a href="alpha_and_beta_diversity/$analysis.$dissimilarity.subsample.braycurtis.$dissimilarity.lt.tre" target="vis">newick format sample relation tree</a>&nbsp;&nbsp;<a href="html/$analysis.$dissimilarity.subsample.braycurtis.$dissimilarity.lt.tre.html" target="vis"><img src="images/chart.gif" height="20" class="chart"></a>
    </ul>
    </div>
    </div>


    <p class="thick">» &beta;-diversity measurements: NMDs analysis &nbsp;&nbsp; <input TYPE="button" class="btn" id="nmd" value="HELP" onClick="help('nmdModal', 'nmd', 'closeNmd');"></p>
    <div id="nmdModal" class="modal">
    <!-- Modal content -->
    <div class="modal-content">
    <span class="close" id="closeNmd">&times;</span>
    <p>non-metric multidimensional scaling (NMDS) tries to preserve the distance between samples using a user defined number of dimensions.The stress file contains the stress and R^2 values, which describe the quality of the ordination. Each line in this file represents a different iteration. The configuration obtained in the iteration with the lowest stress is reported in the axes file. Usually, they suggest that  a stress value below 0.20 is good (less than0.10 is very good). NMDS plots attempt to show ordinal distances between samples as accurately as possible in two dimensions. It is important to report the stress of these plots, because a high stress value means that the algorithm had a hard time representing the distances between samples in 2 dimensions.</p>
    </div>
    </div>



    <div id="wrap">
    <div id="left_col">
    <ul id="lfile">
    <li> Jclass <a href="alpha_and_beta_diversity/$analysis.jclass.$dissimilarity.lt.nmds.stress" target="vis">stress</a>, <a href="alpha_and_beta_diversity/$analysis.jclass.$dissimilarity.lt.nmds.iters" target="vis">iters</a>, and <a href="alpha_and_beta_diversity/$analysis.jclass.$dissimilarity.lt.nmds.axes" target="vis">axes</a> files &nbsp;&nbsp;<a href="html/$analysis.jclass.$dissimilarity.lt.nmds.html" target="vis"><img src="images/chart.gif" height="20" class="chart"></a>
    <li> Jest <a href="alpha_and_beta_diversity/$analysis.jest.$dissimilarity.lt.nmds.stress" target="vis">stress</a>, <a href="alpha_and_beta_diversity/$analysis.jest.$dissimilarity.lt.nmds.iters" target="vis">iters</a>, and <a href="alpha_and_beta_diversity/$analysis.jest.$dissimilarity.lt.nmds.axes" target="vis">axes</a> files &nbsp;&nbsp;<a href="html/$analysis.jest.$dissimilarity.lt.nmds.html" target="vis"><img src="images/chart.gif" height="20" class="chart"></a>
    <li>Thetayc <a href="alpha_and_beta_diversity/$analysis.thetayc.$dissimilarity.lt.nmds.stress" target="vis">stress</a>, <a href="alpha_and_beta_diversity/$analysis.thetayc.$dissimilarity.lt.nmds.iters" target="vis">iters</a>, and <a href="alpha_and_beta_diversity/$analysis.thetayc.$dissimilarity.lt.nmds.axes" target="vis">axes</a> files &nbsp;&nbsp;<a href="html/$analysis.thetayc.$dissimilarity.lt.nmds.html" target="vis"><img src="images/chart.gif" height="20" class="chart"></a>
    <li>Braycurtis <a href="alpha_and_beta_diversity/$analysis.braycurtis.$dissimilarity.lt.nmds.stress" target="vis">stress</a>, <a href="alpha_and_beta_diversity/$analysis.braycurtis.$dissimilarity.lt.nmds.iters" target="vis">iters</a>, and <a href="alpha_and_beta_diversity/$analysis.braycurtis.$dissimilarity.lt.nmds.axes" target="vis">axes</a> files &nbsp;&nbsp;<a href="html/$analysis.jclass.$dissimilarity.lt.nmds.html" target="vis"><img src="images/chart.gif" height="20" class="chart"></a>
    </ul>
    </div>
    <div id="right_col">
     <ul id="lfile">
    <li> Jclass <a href="alpha_and_beta_diversity/$analysis.$dissimilarity.subsample.jclass.$dissimilarity.lt.nmds.stress" target="vis">stress</a>, <a href="alpha_and_beta_diversity/$analysis.$dissimilarity.subsample.jclass.$dissimilarity.lt.nmds.iters" target="vis">iters</a>, and <a href="alpha_and_beta_diversity/$analysis.$dissimilarity.subsample.jclass.$dissimilarity.lt.nmds.axes" target="vis">axes</a> files &nbsp;&nbsp;<a href="html/$analysis.$dissimilarity.subsample.jclass.$dissimilarity.lt.nmds.html" target="vis"><img src="images/chart.gif" height="20" class="chart"></a>
    <li> Jest <a href="alpha_and_beta_diversity/$analysis.$dissimilarity.subsample.jest.$dissimilarity.lt.nmds.stress" target="vis">stress</a>, <a href="alpha_and_beta_diversity/$analysis.$dissimilarity.subsample.jest.$dissimilarity.lt.nmds.iters" target="vis">iters</a>, and <a href="alpha_and_beta_diversity/$analysis.$dissimilarity.subsample.jest.$dissimilarity.lt.nmds.axes" target="vis">axes</a> files &nbsp;&nbsp;<a href="html/$analysis.$dissimilarity.subsample.jest.$dissimilarity.lt.nmds.html" target="vis"><img src="images/chart.gif" height="20" class="chart"></a>
    <li>Thetayc <a href="alpha_and_beta_diversity/$analysis.$dissimilarity.subsample.thetayc.$dissimilarity.lt.nmds.stress" target="vis">stress</a>, <a href="alpha_and_beta_diversity/$analysis.$dissimilarity.subsample.thetayc.$dissimilarity.lt.nmds.iters" target="vis">iters</a>, and <a href="alpha_and_beta_diversity/$analysis.$dissimilarity.subsample.thetayc.$dissimilarity.lt.nmds.axes" target="vis">axes</a> files &nbsp;&nbsp;<a href="html/$analysis.$dissimilarity.subsample.thetayc.$dissimilarity.lt.nmds.html" target="vis"><img src="images/chart.gif" height="20" class="chart"></a>
    <li>Braycurtis <a href="alpha_and_beta_diversity/$analysis.$dissimilarity.subsample.braycurtis.$dissimilarity.lt.nmds.stress" target="vis">stress</a>, <a href="alpha_and_beta_diversity/$analysis.$dissimilarity.subsample.braycurtis.$dissimilarity.lt.nmds.iters" target="vis">iters</a>, and <a href="alpha_and_beta_diversity/$analysis.$dissimilarity.subsample.braycurtis.$dissimilarity.lt.nmds.axes" target="vis">axes</a> files &nbsp;&nbsp;<a href="html/$analysis.$dissimilarity.subsample.braycurtis.$dissimilarity.lt.nmds.html" target="vis"><img src="images/chart.gif" height="20" class="chart"></a>
    </ul>


	</div>
      </div>

    </div>


RESULTS_ALPHA


	}

sub print_hypothesis_testing{
    my($RESULTS, $analysis, $samplenames) = @_;
    print $RESULTS <<RESULTS_END;
    <div id="tabs-nohdr-6">
	<div id="wrap">
        <div id="left_col">
        <h3> Without subsampling or rarefying</h3>
        </div>
        <div id="right_col">
        <h3> With subsampling or rarefying</h3>
        </div>
        </div>

	<p  class="thick">» Parsimony testing &nbsp;&nbsp; <input TYPE="button" class="btn" id="parsimony" value="HELP" onClick="help('parsimonyModal', 'parsimony', 'closeParsimony');"> </p>
    <div id="parsimonyModal" class="modal">
    <!-- Modal content -->
    <div class="modal-content">
    <span class="close" id="closeParsimony">&times;</span>
    <p>In mothur, the parsimony method (aka p-test) is a generic test that describes whether two or more communities have the same structure. The significance of the test statistic can only indicate the probability that the communities have the same structure by chance. The value does not indicate a level of similarity.</p>
    </div>
    </div>



    <div id="wrap">
    <div id="left_col">

    <ul id="lfile">
    <li>Jclass <a href="hypothesis_testing/$analysis.jclass.$dissimilarity.lt.tre.parsimony" target="vis">tre.parsimony</a> & <a href="hypothesis_testing/$analysis.jclass.$dissimilarity.lt.tre.psummary" target="vis">tre.psummary</a> <a href="https://www.mothur.org/wiki/Parsimony" target="vis"><img src="images/help.jpg" height="20" class="help">
    <li>Jest <a href="hypothesis_testing/$analysis.jest.$dissimilarity.lt.tre.parsimony" target="vis">tre.parsimony</a> & <a href="hypothesis_testing/$analysis.jest.$dissimilarity.lt.tre.psummary" target="vis">tre.psummary</a>
    <li>Thetayc <a href="hypothesis_testing/$analysis.thetayc.$dissimilarity.lt.tre.parsimony" target="vis">tre.parsimony</a> & <a href="hypothesis_testing/$analysis.thetayc.$dissimilarity.lt.tre.psummary" target="vis">tre.psummary</a>
    <li>Braycurtis <a href="hypothesis_testing/$analysis.braycurtis.$dissimilarity.lt.tre.parsimony" target="vis">tre.parsimony</a> & <a href="hypothesis_testing/$analysis.braycurtis.$dissimilarity.lt.tre.psummary" target="vis">tre.psummary</a>
    </ul>
    </div>
    <div id="right_col">

    <ul id="lfile">
    <li>Jclass <a href="hypothesis_testing/$analysis.$dissimilarity.subsample.jclass.$dissimilarity.lt.tre.parsimony" target="vis">tre.parsimony</a> & <a href="hypothesis_testing/$analysis.$dissimilarity.subsample.jclass.$dissimilarity.lt.tre.psummary" target="vis">tre.psummary</a>
    <li>Jest <a href="hypothesis_testing/$analysis.$dissimilarity.subsample.jest.$dissimilarity.lt.tre.parsimony" target="vis">tre.parsimony</a> & <a href="hypothesis_testing/$analysis.$dissimilarity.subsample.jest.$dissimilarity.lt.tre.psummary" target="vis">tre.psummary</a>
    <li>Thetayc <a href="hypothesis_testing/$analysis.$dissimilarity.subsample.thetayc.$dissimilarity.lt.tre.parsimony" target="vis">tre.parsimony</a> & <a href="hypothesis_testing/$analysis.$dissimilarity.subsample.thetayc.$dissimilarity.lt.tre.psummary" target="vis">tre.psummary</a>
    <li>Braycurtis <a href="hypothesis_testing/$analysis.$dissimilarity.subsample.braycurtis.$dissimilarity.lt.tre.parsimony" target="vis">tre.parsimony</a> & <a href="hypothesis_testing/$analysis.$dissimilarity.subsample.braycurtis.$dissimilarity.lt.tre.psummary" target="vis">tre.psummary</a>
    </ul>
    </div>
    </div>
    <p  class="thick">» Unifrac.weighted &nbsp;&nbsp; <input TYPE="button" class="btn" id="wunifrac" value="HELP" onClick="help('wunifracModal', 'wunifrac', 'closeWunifrac');"></p>

     <div id="wunifracModal" class="modal">
    <!-- Modal content -->
    <div class="modal-content">
    <span class="close" id="closeWunifrac">&times;</span>
    <p>Mothur implements the weighted UniFrac algorithm.The UniFrac methods are generic tests that describes whether two or more communities have the same structure. The significance of the test statistic can only indicate the probability that the communities have the same structure by chance. The value does not indicate a level of similarity.</p>
    </div>
    </div>




    <div id="wrap">
    <div id="left_col">

    <ul id="lfile">
    <li>Jclass <a href="hypothesis_testing/$analysis.jclass.$dissimilarity.lt.tre1.weighted" target="vis">tre1.weighted</a> & <a href="hypothesis_testing/$analysis.jclass.$dissimilarity.lt.trewsummary" target="vis">trewsummary</a> <a href="https://www.mothur.org/wiki/Unifrac.weighted" target="vis"><img src="images/help.jpg" height="20" class="help">
    <li>Jest <a href="hypothesis_testing/$analysis.jest.$dissimilarity.lt.tre1.weighted" target="vis">tre1.weighted</a> & <a href="hypothesis_testing/$analysis.jest.$dissimilarity.lt.trewsummary" target="vis">trewsummary</a>
    <li>Thetayc <a href="hypothesis_testing/$analysis.thetayc.$dissimilarity.lt.tre1.weighted" target="vis">tre1.weighted</a> & <a href="hypothesis_testing/$analysis.thetayc.$dissimilarity.lt.trewsummary" target="vis">trewsummary</a>
    <li>Braycurtis <a href="hypothesis_testing/$analysis.braycurtis.$dissimilarity.lt.tre1.weighted" target="vis">tre1.weighted</a> & <a href="hypothesis_testing/$analysis.braycurtis.$dissimilarity.lt.trewsummary" target="vis">trewsummary</a>
    </ul>
    </div>
    <div id="right_col">
    <ul id="lfile">
    <li>Jclass <a href="hypothesis_testing/$analysis.$dissimilarity.subsample.jclass.$dissimilarity.lt.tre1.weighted" target="vis">tre1.weighted</a> & <a href="hypothesis_testing/$analysis.$dissimilarity.subsample.jclass.$dissimilarity.lt.trewsummary" target="vis">trewsummary</a>
    <li>Jest <a href="hypothesis_testing/$analysis.$dissimilarity.subsample.jest.$dissimilarity.lt.tre1.weighted" target="vis">tre1.weighted</a> & <a href="hypothesis_testing/$analysis.$dissimilarity.subsample.jest.$dissimilarity.lt.trewsummary" target="vis">trewsummary</a>
    <li>Thetayc <a href="hypothesis_testing/$analysis.$dissimilarity.subsample.thetayc.$dissimilarity.lt.tre1.weighted" target="vis">tre1.weighted</a> & <a href="hypothesis_testing/$analysis.$dissimilarity.subsample.thetayc.$dissimilarity.lt.trewsummary" target="vis">trewsummary</a>
    <li>Braycurtis <a href="hypothesis_testing/$analysis.$dissimilarity.subsample.braycurtis.$dissimilarity.lt.tre1.weighted" target="vis">tre1.weighted</a> & <a href="hypothesis_testing/$analysis.$dissimilarity.subsample.braycurtis.$dissimilarity.lt.trewsummary" target="vis">trewsummary</a>
    </ul>
    </div>
    </div>

    <p  class="thick">» Unifrac.unweighted &nbsp;&nbsp; <input TYPE="button" class="btn" id="uunifrac" value="HELP" onClick="help('uunifracModal', 'uunifrac', 'closeUunifrac');"></p>

     <div id="uunifracModal" class="modal">
    <!-- Modal content -->
    <div class="modal-content">
    <span class="close" id="closeUunifrac">&times;</span>
    <p>Mothur implements the unweighted UniFrac algorithm.The UniFrac methods are generic tests that describes whether two or more communities have the same structure. The significance of the tes\t statistic can only indicate the probability that the communities have the same structure by chance. The value does not indicate a level of similarity.</p>
    </div>
    </div>



<div id="wrap">
    <div id="left_col">

    <ul id="lfile">
    <li>Jclass <a href="hypothesis_testing/$analysis.jclass.$dissimilarity.lt.tre1.unweighted" target="vis">tre1.weighted</a> & <a href="hypothesis_testing/$analysis.jclass.$dissimilarity.lt.uwsummary" target="vis">uwsummary</a><a href="https://www.mothur.org/wiki/Unifrac.unweighted" target="vis"><img src="images/help.jpg" height="20" class="help">
    <li>Jest <a href="hypothesis_testing/$analysis.jest.$dissimilarity.lt.tre1.unweighted" target="vis">tre1.weighted</a> & <a href="hypothesis_testing/$analysis.jest.$dissimilarity.lt.uwsummary" target="vis">uwsummary</a>
    <li>Thetayc <a href="hypothesis_testing/$analysis.thetayc.$dissimilarity.lt.tre1.unweighted" target="vis">tre1.weighted</a> & <a href="hypothesis_testing/$analysis.thetayc.$dissimilarity.lt.uwsummary" target="vis">uwsummary</a>
    <li>Braycurtis <a href="hypothesis_testing/$analysis.braycurtis.$dissimilarity.lt.tre1.unweighted" target="vis">tre1.weighted</a> & <a href="hypothesis_testing/$analysis.braycurtis.$dissimilarity.lt.uwsummary" target="vis">uwsummary</a>
    </ul>
    </div>
    <div id="right_col">
    <ul id="lfile">
    <li>Jclass <a href="hypothesis_testing/$analysis.$dissimilarity.subsample.jclass.$dissimilarity.lt.tre1.unweighted" target="vis">tre1.weighted</a> & <a href="hypothesis_testing/$analysis.$dissimilarity.subsample.jclass.$dissimilarity.lt.uwsummary" target="vis">uwsummary</a>
    <li>Jest <a href="hypothesis_testing/$analysis.$dissimilarity.subsample.jest.$dissimilarity.lt.tre1.unweighted" target="vis">tre1.weighted</a> & <a href="hypothesis_testing/$analysis.$dissimilarity.subsample.jest.$dissimilarity.lt.uwsummary" target="vis">uwsummary</a>
    <li>Thetayc <a href="hypothesis_testing/$analysis.$dissimilarity.subsample.thetayc.$dissimilarity.lt.tre1.unweighted" target="vis">tre1.weighted</a> & <a href="hypothesis_testing/$analysis.$dissimilarity.subsample.thetayc.$dissimilarity.lt.uwsummary" target="vis">uwsummary</a>
    <li>Braycurtis <a href="hypothesis_testing/$analysis.$dissimilarity.subsample.braycurtis.$dissimilarity.lt.tre1.unweighted" target="vis">tre1.weighted</a> & <a href="hypothesis_testing/$analysis.$dissimilarity.subsample.braycurtis.$dissimilarity.lt.uwsummary" target="vis">uwsummary</a>
    </ul>
    </div>
    </div>
    <p  class="thick">» AMOVA(Analysis of Molecular Variance) &nbsp;&nbsp; <input TYPE="button" class="btn" id="amova" value="HELP" onClick="help('amovaModal', 'amova', 'closeAmova');"> </p>
    <div id="amovaModal" class="modal">
    <!-- Modal content -->
    <div class="modal-content">
    <span class="close" id="closeAmova">&times;</span>
    <p>Analysis of molecular variance is a nonparametric analog of traditional analysis of variance. This method is widely used in population genetics to test the hypothesis that genetic diversity within two populations is not significantly different from that which would result from pooling the two populations (Excoffier et al., 1992; Anderson, 2001; Martin, 2002). Amova  tests whether the centers of the clouds representing a group are more separated than the variation among samples of the same treatment.</p>
    </div>
    </div>



    <div id="wrap">
    <div id="left_col">

    <ul id="lfile">
	<li>Jclass <a href="hypothesis_testing/$analysis.jclass.$dissimilarity.lt.amova" target="vis">amova testing result </a> <a href="https://mothur.org/wiki/Amova" target="vis"><img src="images/help.jpg" height="20" class="help">

	<li>Jest <a href="hypothesis_testing/$analysis.jest.$dissimilarity.lt.amova" target="vis">amova testing result </a>
	<li>Thetayc <a href="hypothesis_testing/$analysis.thetayc.$dissimilarity.lt.amova" target="vis">amova testing result </a>
	<li>Braycurtis <a href="hypothesis_testing/$analysis.braycurtis.$dissimilarity.lt.amova" target="vis">amova testing result </a>
    </ul>
    </div>
    <div id="right_col">
    <ul id="lfile">
    <li>Jclass <a href="hypothesis_testing/$analysis.$dissimilarity.subsample.jclass.$dissimilarity.lt.amova" target="vis">amova testing result </a>
    <li>Jest <a href="hypothesis_testing/$analysis.$dissimilarity.subsample.jest.$dissimilarity.lt.amova" target="vis">amova testing result </a>
    <li>Thetayc <a href="hypothesis_testing/$analysis.$dissimilarity.subsample.thetayc.$dissimilarity.lt.amova" target="vis">amova testing result </a>
    <li>Braycurtis <a href="hypothesis_testing/$analysis.$dissimilarity.subsample.braycurtis.$dissimilarity.lt.amova" target="vis">amova testing result </a>
    </ul>
    </div>
    </div>
    <p  class="thick">» HOMOVA (Homogenity of Molecular Variance) &nbsp;&nbsp; <input TYPE="button" class="btn" id="homova" value="HELP" onClick="help('homovaModal', 'homova', 'closeHomova');"></p>

     <div id="homovaModal" class="modal">
    <!-- Modal content -->
    <div class="modal-content">
    <span class="close" id="closeHomova">&times;</span>
    <p>Homogeneity of molecular variance (HOMOVA) is a nonparametric analog of Bartlett’s test for homo- geneity of variance, which has been used in population genetics to test the hypothesis that the genetic diversity within two or more populations is homogeneous (Stewart and Excoffier, 1996);</p>
    </div>
    </div>


    <div id="wrap">
    <div id="left_col">

    <ul id="lfile">
    <li>Jclass <a href="hypothesis_testing/$analysis.jclass.$dissimilarity.lt.homova" target="vis">homova testing result </a> <a href="https://mothur.org/wiki/Homova" target="vis"><img src="images/help.jpg" height="20" class="help">

	<li>Jest <a href="hypothesis_testing/$analysis.jest.$dissimilarity.lt.homova" target="vis">homova testing result </a>
	<li>Thetayc <a href="hypothesis_testing/$analysis.thetayc.$dissimilarity.lt.homova" target="vis">homova testing result </a>
	<li>Braycurtis <a href="hypothesis_testing/$analysis.braycurtis.$dissimilarity.lt.homova" target="vis">homova testing result </a>
      </ul>
    </div>
    <div id="right_col">
    <ul id="lfile">
    <li>Jclass <a href="hypothesis_testing/$analysis.$dissimilarity.subsample.jclass.$dissimilarity.lt.homova" target="vis">homova testing result </a>
    <li>Jest <a href="hypothesis_testing/$analysis.$dissimilarity.subsample.jest.$dissimilarity.lt.homova" target="vis">homova testing result </a>
    <li>Thetayc <a href="hypothesis_testing/$analysis.$dissimilarity.subsample.thetayc.$dissimilarity.lt.homova" target="vis">homova testing result </a>
    <li>Braycurtis <a href="hypothesis_testing/$analysis.$dissimilarity.subsample.braycurtis.$dissimilarity.lt.homova" target="vis">homova testing result </a>
    </ul>
    </div>


</div>
RESULTS_END
}
sub print_web_common_end{
my ($RESULTS) = @_;
    print $RESULTS <<RESULTS_END;
    </article>
	</section>
	<div id="sep"></div>
	<div id="foot" align="center">Problems? Questions? Suggestions? Please contact <a href="mailto:xdong\@ucalgary.ca">Xiaoli Dong</a> or <a HREF="mailto:mstrous\@ucalgary.ca">Marc Strous</a><br><a href="http://www.ucalgary.ca/ebg/">Energy Bioengineering Group</a> in <a href="http://geoscience.ucalgary.ca/">Department of Geoscience</a> at <a href="http://www.ucalgary.ca">Univeristy of Calgary</a></div>

	<div>

	</div>

	</html>
RESULTS_END

}
# The information is used when the results are packed into a subdirectory
# and the Web page of results is built.
sub save_job_info {
    my ($file, $analysis, $email) = @_;

    if(! -e $file){
	my($human_time, $machine_time) = set_time_strings();

	open JOB_INFO, ">$file"
	    or notify_internal_error("Can't open job information file $file to write: $!", "Job information file creation failed");
	print JOB_INFO "$analysis\n$email\n$human_time\n$machine_time\n";
    }
}
# Creates unique string based on time and analysis name to use
# as the name of the server directory to store uploaded files in.
sub set_time_strings {
    my @months = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov iDec);
    my @days = qw(Sun Mon Tue Wed Thu Fri Sat Sun);
    my ($sec, $min, $hour, $mday, $mon, $years_since_1900, $wday, $yday, $isdst) = localtime;
    my $year = $years_since_1900 + 1900;
    my $tzone = $isdst ? "MDT": "MST";
    my $true_mon = $mon + 1;
    my $human_time = sprintf("%s %s %d %02d:%02d:%02d %s %d", $days[$wday], $months[$mon], $mday, $hour, $min, $sec, $tzone, $year);
    my $machine_time = sprintf("%4d-%02d-%02d-%02d%02d", $year, $true_mon, $mday, $hour, $min);
    return ($human_time, $machine_time);
}
#----------------------------------------------------------------------

sub msg {
  my $t = localtime;
  my $line = "[".$t->hms."] @_\n\n";
  print LOG $line if openhandle(\*LOG);
  print STDERR $line unless $quiet;
}

sub status {
  my $t = localtime;
  my $line = "[".$t->hms."] @_\n\n";
  print STATUS $line if openhandle(\*STATUS);
  print STDERR $line unless $quiet;
}

#----------------------------------------------------------------------

sub err {
  $quiet=0;
  msg(@_);
  exit(2);
}

sub runcmds {

     my ($max_cmds, @cmds) = @_;

     my ($num_children, $pid);

     return unless @cmds; # get out of the function if there is nothing in @cmds

     for($num_children = 1; $num_children < $max_cmds && @cmds; $num_children++){
	 # initialize the number of child processes at 1, and increment it by one
	 #while it is less than $max_cmds
	 my $cmd = shift (@cmds);
	 if($pid = fork) {
	     # do nothing if parent
	 } elsif(defined $pid) { # $pid is zero here if defined
	     msg("Running:", $cmd);
	     system($cmd)==0 or err("Could not run command:", $cmd);

	     exit;
	 } else {
	     #weird fork error
	     die "Can't fork: $!\n";
	 }
     }

     while(@cmds) {
	 undef $pid;
       FORK: {
	   my $cmd = shift (@cmds);
	   if($pid = fork) {
	       # parent here
	       $num_children++;
	       wait;
	       $num_children--;
	       next;

	   } elsif(defined $pid) { # $pid is zero here if defined

	       msg("Running:", $cmd);
	       system($cmd)==0 or err("Could not run command:", $cmd);
	       exit;

	   } else {
	       #weird fork error
	       die "Can't fork: $!\n";
	   }
	 }
     }
     wait while $num_children--;
}
sub delfile {
    for my $file (@_) {
	msg("Deleting unwanted file:", $file);
	unlink $file or warn "Could not unlink $file: $!\n";
    }
}

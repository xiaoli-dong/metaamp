#!/usr/bin/env perl
use warnings;
use strict;
use FindBin;
use lib "$FindBin::Bin";
require "util.pl";

our $bin = "$FindBin::RealBin";
our $proj = $bin;
$proj =~ s/\/bin//;
our $dbdir = "$proj/database";
our $usearch = "$proj/programs/usearch64";
our $mothur = "$proj/programs/mothur/mothur";
our $template = "$dbdir/silva.nr.fasta";
our $taxonomy = "$dbdir/silva.nr.tax";
our %mapHash = ();
our %fname2sname = ();
our $fprimers = "";
our $rprimers = "";
our $cpus = 1;

sub getTaxonInfo{
    my($analysisname,$cpus,$label, $fasta) = @_;
    my $prefix = $label eq "asv" ? "asv" : "otu";
    my $cmd = "$mothur \"\#classify.seqs(fasta=$fasta,template=$template,taxonomy=$taxonomy,method=wang,ksize=8, iters=100, cutoff=80, processors=$cpus);rename.file(taxonomy=current,prefix=$analysisname.$prefix)\";";

    $cmd .= "$bin/getTaxonSummary.pl -i1 $analysisname.$prefix.taxonomy -i2 $analysisname.shared > $analysisname.taxonomy.summary.txt;";
    $cmd .= "$bin/getOTUTaxonTable.pl -i1 $analysisname.$prefix.taxonomy  -i2 $analysisname.shared -l $label > $analysisname.$prefix.taxonomy.summary.txt;";
    $cmd .= "$bin/getOTUTaxonTable.relabund.pl -i1 $analysisname.$prefix.taxonomy  -i2 $analysisname.shared -l $label > $analysisname.$prefix.taxonomy.summary.relabund.txt";
    runcmds($cpus, $cmd);


}
#get primer sequences
sub getPrimerInfo{
    my ($oligos) = @_;

    my $err = "";

    open(OLIGOS, $oligos) or die "Could not open $oligos to read, $!\n";
    while(<OLIGOS>){
	chomp;
	next if /^\#/;
	if(/forward\s+(\w+)/i){
	    my $fprimer = $1;
	    #valied primer
	    if($fprimer =~ /[abcdghmnrstuvwxy]/gmi){
		$fprimers .= ":$fprimer";
	    }
	    else{
		$err = "user supplied forward primer contains invalid character: $fprimer";
		last;
	    }
	}
	elsif(/reverse\s+(\w+)/i){
	    my $rprimer = $1;

            #valied primer
	    if($rprimer =~ /[abcdghmnrstuvwxy]/gmi){
		$rprimers .= ":$rprimer";
	    }
	    else{
		$err = "user supplied reverse primer contains invalid character: $rprimer";
		last;
	    }
	}
    }
    close(OLIGOS);

    $fprimers =~ s/^:// if $fprimers ne "";
    $rprimers =~ s/^:// if $rprimers ne "";

    return $err;
}

sub makeDesignFile{
    my ($mapHash) = @_;

#generate design file for the sample comparsion and multivariate analysis
    open(DESIGN, ">design.txt") or die "Could not open design.txt to write, $!\n";
    foreach my $samplename (keys %$mapHash){
	my $treatment = $mapHash->{$samplename}->{treatment};
	print DESIGN "$samplename\t$treatment\n";
    }
    close(DESIGN);
}
sub readMappingFile{
    my ($map, $seqtype, $seqformat, $mapHash) = @_;

    open(INPUT, ">raw_file.list.txt") or die "Could not open raw_file.list.txt to write, $!\n";
    open(MAPPING, "$map") or die "Could not open $map to read, $!\n";

#Samplename Treatment Strand Read1 Read2 runid
#Samplename Treatment Strand Read1 runid
#runid is optional. When using dada2,, each run the error profile can be different.
# with the runid supplied, each run will be process separatly. Otherwise, it will be
#process together

    my %fws = ();
    my %rws = ();
    my $err = "";
    my %samplenames = ();
    while(<MAPPING>){
	chomp;
	#skip the header of the mapping file
	next if /^#/;
	next if !/\S/;
	my @sample_prop = ();

	my @line = split(/\s+/, $_);

	my $name = $line[0];
	$name =~ tr/[ \-]/_/;
	my $treatment = $line[1];
	$treatment =~ tr/[ \-]/_/;
	push(@sample_prop, $treatment);


	my $strand = $line[2];

	if($name =~ /^\d/){
	    $err = "Sample name cannot start with number";
	    last;
	}
	if(exists $mapHash->{$name}){
	    $err = "Samplename in mapping file must be unique: $line[0] has duplicate entries";
	    last;
	}
	if($strand !~ /[+-]/){
	    $err = "Sequnce strand can only be positive strand (+) or negative strand (-)";
	    last;
	}
	$mapHash->{$name}->{treatment} = $treatment;
	$mapHash->{$name}->{strand} = $strand;


	if($seqtype eq "single"){
	    if(@line < 4){
		$err = "The mapping file for single-end data must have at least 4 colums";
		last;
	    }
	    my $read1 = $line[3];

	    if(exists $fws{$read1}){
		$err = "Read entries in mapping file must be unique: $read1 has duplicate entries in the mapping file";
		last;
	    }
	    if(! -e $read1){
		$err = "Read1 $read1 does not exist";
		last;
	    }
	    $fws{$read1} = 1;
	    $mapHash->{$name}->{read1} = $read1;
	    $fname2sname{$read1} = $name;

	    print INPUT $read1, "\n";

	    if(@line > 4){
		my $runid = $line[4];
		if($runid =~ /\S/){
		    $mapHash->{$name}->{runid} = $runid;
		}
	    }

	    #convert fasta to fastq
	    if($seqformat eq "fasta"){
		my ($read1_prefix) = $read1 =~ /(.*)\.[^.]+$/;
		my $read1_fastq = "$read1_prefix.fastq";
		my $qual1 = "$read1_prefix\.qual";
		if( -e $qual1){
		    my $cmd = "fasta2fastq.pl -fasta $read1 -qual $qual1 > $read1_fastq";
		    runcmds($cpus, $cmd);
		}
		else{
		     my $cmd = "fasta2fastq.pl -fasta $read1  > $read1_fastq";
		     runcmds($cpus, $cmd);
		}
		$mapHash{$name}{read1} = $read1_fastq;
		$fname2sname{$read1_fastq} = $name;
	    }

	}
	elsif($seqtype eq "paired"){

	    if(@line < 5){
		$err = "The mapping file for paired-end data must have at least 5 colums";
		last;
	    }

	    my $read1 = $line[3];
	    my $read2 = $line[4];

	    if(exists $fws{$read1}){
		$err = "Read entries in mapping file must be unique: $read1 has duplicate entries in the mapping file";
		last;
	    }
	    if(exists $rws{$read2}){
		$err = "Read entries in mapping file must be unique: $read2 has duplicate entries in the mapping file";
		last;
	    }
	    if(! -e $read1){
		$err = "Read1 $read1 does not exist";
		last;
	    }
	    if(! -e $read2){
		$err = "Read2 does not exist";
		last;
	    }

	    $fws{$read1} = 1;
	    $rws{$read2} = 1;
	    $mapHash->{$name}->{read1} = $read1;
	    $mapHash->{$name}->{read2} = $read2;
	    $fname2sname{$read1} = $name;
	    $fname2sname{$read2} = $name;
	    print INPUT $read1, "\n";
	    print INPUT $read2, "\n";
	    if(@line > 5){
		my $runid = $line[5];
		if($runid =~ /\S/){
		    $mapHash->{$name}->{runid} = $runid;
		}
	    }

	    #convert fasta to fastq
	    if($seqformat eq "fasta"){
		my ($read1_prefix) = $read1 =~ /(.*)\.[^.]+$/;
		my $read1_fastq = "$read1_prefix.fastq";
		my $qual1 = "$read1_prefix\.qual";
		if( -e $qual1){
		    my $cmd = "fasta2fastq.pl -fasta $read1 -qual $qual1 > $read1_fastq";
		    runcmds($cpus, $cmd);
		}
		else{
		    my $cmd = "fasta2fastq.pl -fasta $read1  > $read1_fastq";
		    runcmds($cpus, $cmd);
		}
		$mapHash{$name}{read1} = $read1_fastq;
		$fname2sname{$read1_fastq} = $name;
		my ($read2_prefix) = $read2 =~ /(.*)\.[^.]+$/;
		my $read2_fastq = "$read2_prefix.fastq";
		my $qual2 = "$read2_prefix\.qual";
		if( -e $qual2){
			my $cmd = "fasta2fastq.pl -fasta $read2 -qual $qual2 > $read2_fastq";
			runcmds($cpus, $cmd);
		}
		else{
		    my $cmd = "fasta2fastq.pl -fasta $read2 > $read2_fastq";
		    runcmds($cpus, $cmd);
		}
		$mapHash{$name}->{read2} = $read2_fastq;
		$fname2sname{$read2_fastq} = $name;
	    }



	}



    }
    close(INPUT);
    return $err;
}


sub trimPrimers_dada2{
    my($mapHash, $seqtype, $fprimers, $rprimers) = @_;
    my %runids = ();
    #perl ../../bin/trimPrimer.pl -seqtype paired -s test -r1 rep1_C1_S244_L001_R1_001.fastq -r2 rep1_C1_S244_L001_R2_001.fastq -fp CCTACGGGAGGCAGCAG -rp GACTACHVGGGTATCTAATCC
    open(INPUT, ">after_strip_primer_file.list.txt") or die "Could not open after_strip_primer_file.list.txt to write, $!\n";

    foreach my $samplename (keys %$mapHash){
	my $runid = "NORUNID";

	if(exists $mapHash->{$samplename}->{runid}){
	    $runid = $mapHash->{$samplename}->{runid};
	}
	$runids{$runid}++;
	if($seqtype eq "paired"){
	    my $read1 = $mapHash->{$samplename}->{read1};
	    my $read2 = $mapHash->{$samplename}->{read2};
	    my $cmd = "$bin/trimPrimer.pl -seqtype $seqtype -s $samplename -r1 $read1 -r2 $read2";
	    $cmd .= " -fp $fprimers" if $fprimers ne "";
	    $cmd .= " -rp $rprimers" if $rprimers ne "";
	    runcmds($cpus, $cmd);

	    $cmd = "$bin/addBarcodelabel.pl  $samplename.R1.trim.fastq $samplename > $samplename\_R1.NOPRIM.$runid.fastq;";
	    $cmd .= "$bin/addBarcodelabel.pl  $samplename.R2.trim.fastq $samplename > $samplename\_R2.NOPRIM.$runid.fastq;";
	    runcmds($cpus, $cmd);

	    $mapHash->{$samplename}->{read1} = "$samplename\_R1.NOPRIM.$runid.fastq";
	    $mapHash->{$samplename}->{read2} = "$samplename\_R2.NOPRIM.$runid.fastq";
	    print INPUT "$samplename\_R1.NOPRIM.$runid.fastq\n", "$samplename\_R2.NOPRIM.$runid.fastq\n";
	}
	elsif($seqtype eq "single"){
	    my $read1 = $mapHash->{$samplename}->{read1};
	    my $cmd = "$bin/trimPrimer.pl -seqtype $seqtype -s $samplename -r1 $read1";
	    $cmd .= " -fp $fprimers" if $fprimers ne "";
	    $cmd .= " -rp $rprimers" if $rprimers ne "";
	    runcmds($cpus, $cmd);

	    $cmd = "$bin/addBarcodelabel.pl  $samplename.trim.fastq $samplename > $samplename.NOPRIM.$runid.fastq;";
	    runcmds($cpus, $cmd);
	    $mapHash->{$samplename}->{read1} = "$samplename.NOPRIM.$runid.fastq";
	    print INPUT "$samplename.NOPRIM.$runid.fastq\n"
	}
    }
    close(INPUT);
    return \%runids;
}
sub trimPrimers_uparse{
    my($mapHash, $seqtype, $fprimers, $rprimers) = @_;

    #perl ../../bin/trimPrimer.pl -seqtype paired -s test -r1 rep1_C1_S244_L001_R1_001.fastq -r2 rep1_C1_S244_L001_R2_001.fastq -fp CCTACGGGAGGCAGCAG -rp GACTACHVGGGTATCTAATCC
    open(INPUT, ">after_strip_primer_file.list.txt") or die "Could not open after_strip_primer_file.list.txt to write, $!\n";
    foreach my $samplename (keys %$mapHash){
	if($seqtype eq "paired"){
	    my $read1 = $mapHash->{$samplename}->{read1};
	    my $read2 = $mapHash->{$samplename}->{read2};
	    my $cmd = "$bin/trimPrimer.pl -seqtype $seqtype -s $samplename -r1 $read1 -r2 $read2";
	    $cmd .= " -fp $fprimers" if $fprimers ne "";
	    $cmd .= " -rp $rprimers" if $rprimers ne "";
	    runcmds($cpus, $cmd);

	    $cmd = "$bin/addBarcodelabel.pl  $samplename.R1.trim.fastq $samplename > $samplename\_R1.NOPRIM.fastq;";
	    $cmd .= "$bin/addBarcodelabel.pl  $samplename.R2.trim.fastq $samplename > $samplename\_R2.NOPRIM.fastq;";
	    runcmds($cpus, $cmd);

	    $mapHash->{$samplename}->{read1} = "$samplename\_R1.NOPRIM.fastq";
	    $mapHash->{$samplename}->{read2} = "$samplename\_R2.NOPRIM.fastq";
	    print INPUT "$samplename\_R1.NOPRIM.fastq\n", "$samplename\_R2.NOPRIM.fastq\n";
	}
	elsif($seqtype eq "single"){
	    my $read1 = $mapHash->{$samplename}->{read1};
	    my $cmd = "$bin/trimPrimer.pl -seqtype $seqtype -s $samplename -r1 $read1";
	    $cmd .= " -fp $fprimers" if $fprimers ne "";
	    $cmd .= " -rp $rprimers" if $rprimers ne "";
	    runcmds($cpus, $cmd);

	    $cmd = "$bin/addBarcodelabel.pl  $samplename.trim.fastq $samplename > $samplename.NOPRIM.fastq;";
	    runcmds($cpus, $cmd);
	    $mapHash->{$samplename}->{read1} = "$samplename.NOPRIM.fastq";
	    print INPUT "$samplename.NOPRIM.fastq\n"
	}
    }
    close(INPUT);

}




# Sends an email to the user notifying the job is finished
sub old_send_job_finished_email {

    my $results_dir = "results";

    open(JOB_INFO, "JobInfo.txt")
	or die "Can't open job information file for reading: $!";

    chomp(my $analysis = <JOB_INFO>);
    chomp(my $email = <JOB_INFO>);
    chomp(my $human_time = <JOB_INFO>);
    chomp(my $machine_time = <JOB_INFO>);
    close JOB_INFO;

    open(MAIL, "|/usr/sbin/sendmail -t");
    print MAIL "To: $email\n";
    #print MAIL "Bcc: xdong\@ucalgary.ca\n";
    #print MAIL "From: metaAmp\@ebg.ucalgary.ca\n";
    print MAIL "From: xdong\@ucalgary.ca\n";
    print MAIL "Subject: metaAmp job done: $analysis\n\n";
    print MAIL "Dear metaAmp user,\n\nYour metaAmp job with the analysis name: $analysis has been finished. Please visit the link below to view or download your analysis results:\n\nhttp://ebg.ucalgary.ca/metaamp/tmp/$machine_time-$analysis/$results_dir/index.html\n\nThank you for using metaAmp.\n";
    close(MAIL);
}
# The information is used when the results are packed into a subdirectory
# and the Web page of results is built.
sub save_job_info {
    my ($file, $analysis, $email, $param) = @_;
    my($human_time, $machine_time) = set_time_strings();
    open JOB_INFO, ">$file" or die "Can't open job information file $file to write: $!", "Job information file creation failed";
    print JOB_INFO "Job name:$analysis\nEmail:$email\nSubmit human time:$human_time\nSubmit machine time:$machine_time\nCommand:$param\n";
}
# Creates unique string based on time and analysis name to use
# as the name of the server directory to store uploaded files in.
sub set_time_strings {
    my @months = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
    my @days = qw(Sun Mon Tue Wed Thu Fri Sat Sun);
    my ($sec, $min, $hour, $mday, $mon, $years_since_1900, $wday, $yday, $isdst) = localtime;
    my $year = $years_since_1900 + 1900;
    my $tzone = $isdst ? "MDT": "MST";
    my $true_mon = $mon + 1;
    my $human_time = sprintf("%s %s %d %02d:%02d:%02d %s %d", $days[$wday], $months[$mon], $mday, $hour, $min, $sec, $tzone, $year);
    my $machine_time = sprintf("%4d-%02d-%02d-%02d%02d", $year, $true_mon, $mday, $hour, $min);
    return ($human_time, $machine_time);
}


1;

__END__

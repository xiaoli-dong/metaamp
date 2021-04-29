#!/gpfs/ebg_data/lib/perl/bin/perl

# The environment variable defined here are propagated to all Perl scripts in
# $ENV{METAAMPHOME}/bin directory through the "env perl" command when invoking
# the scripts.
BEGIN {
    # Web document root directory
    $ENV{METAAMPHOME} = "/export/web/metaamp";
    $ENV{PERL5LIB} = "/export/home/xdong/perl5/lib/perl5";
    #get rid of "TERM environment variable not set." message
    $ENV{TERM}="xterm";
    $ENV{PATH} = "$ENV{METAAMPHOME}/bin:$ENV{METAAMPHOME}/programs:$ENV{METAAMPHOME}/programs/mothur:/gpfs/ebg_data/programs/cutadapt/bin:$ENV{PATH}";

}
use lib "/gpfs/home/xdong/perl5/lib/perl5" ;
use strict;
use warnings;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
use File::Basename;
use File::Copy;
use File::Path qw(remove_tree);
use Cwd;
use Archive::Extract;
use File::Path qw(make_path remove_tree);
use Cwd qw(getcwd);
use Encode;

# need either this or to explicitly flush stdout, etc.
$| = 1;


# In order to prevent the server being overloaded by huge file uploads,
# we limit the allowable size of the upload to 15000MB
$CGI::POST_MAX = 1048576 * 15000;
my $query = new CGI;
my $survive = 0;

# Send HTTP header to show subsequent messages as HTML
print $query->header();

my $analysis = $query->param("analysisname");
my $email = $query->param("email");
my $fprimers = $query->param("fprimers");
my $rprimers = $query->param("rprimers");
my $seqtype = $query->param("seqtype");
my $seqformat = $query->param("seqformat");
my $minoverlen = $query->param("minoverlen");
my $maxdiffs = $query->param("maxdiffs");
my $maxee = $query->param("maxee");
my $trunclen = $query->param("trunclen");
#my $action = $query->url_param("action");
my $pdiffs = $query->param("pdiffs");
my $scutoff = $query->param("similarity");

# Get file upload parameters
my $seqArchiveFile = $query->param("seqArchiveFile");
my $mappingFile = $query->param("mappingFile");
my $url = "http://ebg.ucalgary.ca/metaamp";
my $template_dir = "$ENV{METAAMPHOME}/template";
# Catch CGI error when the post size exceeds POST_MAX and other unexpected errors happen
catch_cgi_error();
my @seqfiles;
my @fastq_bases;
my $num_seqfiles;

my $human_time;
my $machine_time;

umask 022;

########## Single archive file uploaded ##########

make_and_move_to_upload_dir();
my $archive = convert_client_filename($seqArchiveFile, "archive");
upload_files(["seqArchiveFile"], $archive);
extract_archive($archive, $seqformat);

my $upload_mapping_file = convert_client_filename($mappingFile, "mapping");
upload_files(["mappingFile"], $upload_mapping_file);
#verify the file exist
open(MAP, $upload_mapping_file) or die "Could not open $upload_mapping_file to read, $!\n";

while(<MAP>){
    chomp;
    next if (/^#/);
    tr/\r/\n/;
    next if (!/\S+/);
    my @line = split(/\s+/, $_);

    if(! -e $line[3]) {
	print_error_web("File is not existing","Please make sure the read1 name $line[3] in mapping file has the same name as the sequence file you included in the upload archived sequence file", "Cannot fine read1 file: $line[3] in the uploaded archived sequence file");
    }
    if(@line == 5 && (! -e $line[4])) {
	print_error_web("File is not existing","Please make sure the read2 name $line[4] in mapping file has the same name as the sequence file you included in the upload archived sequence file", "Cannot find read2 file: $line[4] in the uploaded archived sequence file");

    }

}

createOligoFile($fprimers, $rprimers);


my $cmd = "metaamp.pl -map $upload_mapping_file -an $analysis -seqtype $seqtype -seqformat $seqformat -s $scutoff -minoverlen $minoverlen -maxdiffs $maxdiffs -trunclen $trunclen -maxee $maxee -oligos oligos.txt -pdiffs $pdiffs -email $email  >& debug.txt";
my $metaamp_param = $cmd;

save_job_info("JobInfo.txt");

my $outDir = "results";

make_path $outDir;
chmod 0777, $outDir;
my $outLink = "$analysis\_$outDir";
symlink $outDir, $outLink;
copy "../../html/dummyResult.html", "$outDir/index.html";
print_job_received_html();


###############fork and return the browser to user
if ( !defined(my $pid = fork())) {
    print STDERR "fork error\n";
    exit(1);
}
elsif ($pid == 0) {
# child
    close(STDOUT);close(STDIN);close(STDERR);
    #exec('./fork-long-process-test-process');     # lengthy processing

    #You'll need to bitshift the return value by 8 (or divide by 256) to get the return value of the program called:
    my $code = system($cmd) >> 8;
    if($code != 0){
	if($code == 3){notify_internal_error("None of reads has passed quality control stage of the MetaAmp pipeline\nFor the detail reason please check quality control status file at:\n$url/tmp/$machine_time-$analysis/MetaAmp.qcstatus.txt\n")}
	elsif($code == 4){notify_internal_error("Mapping file is not properly formated, check help page for how to make a mapping file\n")}
    }
}
else {
# parent
    #print "forked child \$pid= $pid\n";
    exit 0;
}
###################end##########################


# Creates unique string based on time and analysis name to use
# as the name of the server directory to store uploaded files in.
sub set_time_strings {
    my @months = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov iDec);
    my @days = qw(Sun Mon Tue Wed Thu Fri Sat Sun);
    my ($sec, $min, $hour, $mday, $mon, $years_since_1900, $wday, $yday, $isdst) = localtime;
    my $year = $years_since_1900 + 1900;
    my $tzone = $isdst ? "MDT": "MST";
    my $true_mon = $mon + 1;
    $human_time = sprintf("%s %s %d %02d:%02d:%02d %s %d", $days[$wday], $months[$mon], $mday, $hour, $min, $sec, $tzone, $year);
    $machine_time = sprintf("%4d-%02d-%02d-%02d%02d", $year, $true_mon, $mday, $hour, $min);
}



# Uploads files input in the form by copying them into the
# upload directory created.
sub upload_files {
    my ($ref_form_input_names, @files) = @_;
    my @form_input_names = @{$ref_form_input_names};
    for my $i (0..$#files) {
	my $upload_filehandle = $query->upload($form_input_names[$i]);
	open UPLOAD, ">$files[$i]"
	    or notify_internal_error("Can't open file $files[$i] to write: $!", "File upload failed");
	binmode UPLOAD;
	while (<$upload_filehandle>) {
	    print UPLOAD;
	}
	close UPLOAD;
    }
}

# Prints an error message HTML, telling the user about the problem.
# Usually asks the user to retry.
sub print_error_web {
    my ($heading, $instruction, @messages) = @_;
    print "<html>\n<head>\n";
    print "<link rel=stylesheet type=\"text/css\" href=\"../css/metaamp.css\"/>\n";
    print "<title>metaamp problem</title>\n</head>\n";
    print "<body>\n<h1>$heading</h1>\n";
    for my $i (0..$#messages) {
	print "<p><font color=\"red\">$messages[$i]</font></p>\n";
    }
    print "<p>$instruction</p>\n";
    print "<!-- form><input type=\"button\" value=\"Retry\" style=\"height:2.5em; width:6em\" onClick=\"history.go(-1);return true;\"></form -->\n";
    print "</body>\n</html>\n";
    exit;
}

# Make a directory and move to it to prepare for writing files in it.
sub make_and_move_to_upload_dir {
    set_time_strings();
    my $dir = "$ENV{METAAMPHOME}/tmp/$machine_time-$analysis";
    mkdir $dir
	or notify_internal_error("Can't create upload directory $dir: $!", "Upload directory $dir creation failed");
    chmod 0777, $dir;
    chdir $dir;
}



# Prints out job submission confirmation HTML, listing all job.
# information and files received
sub print_job_received_html {

    use File::Slurper 'read_text';
    my $link_str = read_text("$template_dir/result_link.html");
    $link_str =~ s/METAAMPTYPE/OTU/g;
    $link_str =~ s/METAAMP_ANALYSISNAME/$analysis/g;
    $link_str =~ s/METAAMP_MACHINE_TIME/$machine_time/g;
    $link_str =~ s/METAAMP_USER_EMAIL/$email/g;
    $link_str =~ s/METAAMP_SUBMIT_TIME/$human_time/g;
    my $result_url = "$url/tmp/$machine_time\-$analysis/results/index.html";
    $link_str =~ s/METAAMP_RESULT_URL/$result_url/g;
    my $debug_url = "$url/tmp/$machine_time\-$analysis/debug.txt";
       $link_str =~ s/METAAMP_DEBUG_URL/$debug_url/g;
    $link_str =~ s/METAAMPFASTQXOUNT/$num_seqfiles/g;
    #binmode(OUT, ":utf8");
    binmode STDOUT, ":utf8";
    print $link_str;
    #close(OUT);

}

# Catches CGI error.
# Currently the only known one is the POST size exeeding the limit.
sub catch_cgi_error {
    my $error = $query->cgi_error;
    if ($error) {
	print_error_web("metaamp job request has not been processed", "Please send email to the metaamp Team giving the error message below:<br>$query->strong($error)", "metaamp server has encountered a problem.");
    }
}

# Saves necessary job information into a text file.
# The information is used when the results are packed into a subdirectory
# and the Web page of results is built.
sub save_job_info {
    my ($file) = @_;
    open JOB_INFO, ">$file"
	or notify_internal_error("Can't open job information file $file to write: $!", "Job information file creation failed");
    #print JOB_INFO "$analysis\n$email\n$human_time\n$machine_time\n";
    print JOB_INFO "Job name:$analysis\nEmail:$email\nSubmit human time:$human_time\nSubmit machine time:$machine_time\nCommand:$metaamp_param\n";
}

sub convert_client_filename {
    my ($client_filename, $file_type) = @_;
    my ($base, $path, $ext) = fileparse($client_filename, '\.[^.]+$');
    my $converted_filename = safe_filename("$base$ext");
    if ($converted_filename) {
	return $converted_filename;
    }
    else {
	print_error_web("metaamp has a problem with the submitted $file_type file", "Please ensure your $file_type file name contains only alphanumeric characters, \".\", \"-\", \"_\", and no space.", "Submitted $file_type file \"$base$ext\" contains unsafe characters.");
    }
}

# Checks the archive file and extract all files from it if possible.
# Checks the extracted files for safe file names and puts them into
# fasta and quality filename/basename arrays, based on file extensions.
sub extract_archive {
    my ($archive, $seqformat) = @_;
    my $ae = Archive::Extract->new(archive => $archive);
    if ($ae) {
	#($ae->is_tar || $ae->is_tgz || $ae->is_gz || $ae->is_zip || $ae->is_bz2 || $ae->is_tbz || $ae->is_lzma || $ae->is_xz || ae->is_txz)
	$ae->extract;
	for my $file (@{$ae->files}) {
	    #print "Extracted: $file<br>";
	    # Strip leading / if it exists (packed with absolute path name)
	    $file =~ s/^\///;
	    # If this entry is a directory, change permission (for later local access) and skip
	    if (-d $file) {
		chmod 0755, $file;
		next;
	    }
	    # Skip if this is a bogus entry (e.g. those starting with __MACOSX/)
	    unless (-e $file) {
		next;
	    }
	    chmod 0644, $file;
	    my ($base, $path, $ext) = fileparse($file, '\.[^.]+$');
	    # Skip if this file starts with "." (e.g. those added by Macs)
	    if ($base eq "" || $base =~ /^\./) {
		next;
	    }
	    #path=mock_data/, base=rep2_C3_S246_L001_R2_001.fastq, ext=.gz
	    #print "path=$path, base=$base, ext=$ext<br>";
	    my $file_safe = safe_filename("$base$ext");
	    if ($file_safe) {
		# If file was extracted into a subdirectory or file name needs to be changed
		if ($path ne "./" || $file_safe ne "$base$ext") {
		    move $file, $file_safe;
		}
	    }
	    else {
		print_error_web("metaamp has a problem with the submitted data files", "Please ensure your data file names contain only alphanumeric characters, \".\", \"-\", \"_\", and no space.", "Data file \"$base$ext\" in the submitted archive contains unsafe characters.");
	    }

	    #if ($seqformat eq "fastq" && $ext =~ /^\.fastq/) {
	    if ($seqformat eq "fastq" && $file =~ /\.fastq/) {
		push @seqfiles, $file_safe;
		push @fastq_bases, safe_filename($base);
	    }
	    elsif($seqformat eq "fasta" && $file =~ /\.fna|\.fasta|\.qual/) {
                push @seqfiles, $file_safe;
                push @fastq_bases, safe_filename($base);
            }
	    else {
		if($seqformat eq "fasta"){
		    print_error_web("Uploaded file extensions are not fasta, fna, qual",
				    "You selected $seqformat, Please ensure your fasta file extension as \".fasta, fasta.gz, fna, fna.gz, qual, qual.gz\"",
				    "Extracted file $file has unknown extension \"$ext\".");
		}
		elsif($seqformat eq "fastq"){
                    print_error_web("Uploaded files are not having fastq extenson",
				    "You selected $seqformat, Please ensure your fastq file extension as \".fastq, .fastq.gz\"",
				    "Extracted file $file has unknown extension \"$ext\".");
                }
	    }
	}
	$num_seqfiles = @seqfiles;

	# Remove any subdirectory created during archive extraction to
	# avoid error from cd-hit-est when it attempts to open a file
	# with an identical name
	for my $file (glob "*") {
	    remove_tree($file) if -d $file;
	}
    }
    else {
	print_error_web("metaamp has a problem with the submitted archive file", "Please ensure your archive is in one of these formats:<ul><li>tar: Standard tar file (.tar)</li><li>tgz: Gzip compressed tar file (.tgz or .tar.gz)</li><li>gz: Gzip compressed file (.gz)</li><li>Z: Lempel-Ziv compressed file (.Z)</li><li>zip: Zip compressed file (.zip, .jar or .par)</li><li>bz2: Bzip2 compressed file (.bz2)</li><li>tbz: Bzip2 compressed tar file (.tbz or .tar.bz2)</li><li>lzma: Lzma compressed file (.lzma)</li><li>xz: Xz compressed file (.xz)</li><li>txz: Xz compressed tar file (.txz or .tar.xz)</li></ul>", "metaamp can't extract the archive file $archive.");
	exit;
    }
}

# Removes potentially unsafe characters from a file name and returns
# the converted one. If the conversion fails somehow (usually by
# becoming an empty string), returns false.
sub safe_filename {
    my ($file) = @_;
    my $safe_filename_characters = "a-zA-Z0-9_.-";
    # Space is bad, and mothur breaks file names around "-" (it shouldn't)
    #$file =~ tr/[ \-]/_/;
    #$file =~ s/[^$safe_filename_characters]//g;
    if ($file =~ /^([$safe_filename_characters]+)$/) {
	return $1;
    }
    else {
	return 0;
    }
}

# When an error internal to the server happens, records the message
# in a log and sends email to the person in charge of metaamp web.
# Exits this CGI unless a flag is given to continue.
sub notify_internal_error {
    my ($internal_message, $external_message, $survive) = @_;
    my $cgi_log = "cgi.log";
    if (open CGI_LOG, ">>$cgi_log") {
	print CGI_LOG "$internal_message\n";
	close CGI_LOG;
    }
    if (open MAIL, "|/usr/sbin/sendmail -t") {
	print MAIL "To: $email\n";
	#print MAIL "Bcc: xdong\@ucalgary.ca\n";
	print MAIL "From: metaamp\@ebg.ucalgary.ca\n";
	print MAIL "Subject: metaamp job $analysis problem\n\n";
	print MAIL "Dear metaamp user,\n\nmetaamp job $analysis submitted by $email at $human_time had a problem:\n\n$internal_message";
	close MAIL;
    }
    unless ($survive) {
	print_error_web("metaamp server has a problem", "metaamp Team has been notified. Please retry when you receive a messsage that this problem has been fixed.", $external_message);
	exit;
    }
}

sub createOligoFile{
    my ($forward,$reverse) = @_;
    open OLIGOS, ">oligos.txt"
	or notify_internal_error("Can't open oligos.txt to write: $!", "createOligoFile failed");
    chomp($forward);
    chomp($reverse);
    print STDERR
	my @fa = split(/\n/, $forward);

    foreach (@fa){
	next if !/\w/;
	s/\s+//g;
	next if /[^acgturymkwsbdhvnACGTURYMKWSBDHVN]/;
	print OLIGOS "forward\t$_\n";
    }
    my @ra = split(/\n/, $reverse);
    foreach (@ra){
	next if !/\w/;
	s/\s+//g;
	next if /[^acgturymkwsbdhvnACGTURYMKWSBDHVN]/;
	print OLIGOS "reverse\t$_\n";
    }
    close(OLIGOS);
}


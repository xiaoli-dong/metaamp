#!/usr/bin/env perl

use strict;
use warnings;
use Time::Piece;
use Scalar::Util qw(openhandle);

my $quiet = 0;
my $logfile = "MetaAmp.log.txt";
open LOG, '>', $logfile or err("Can't open logfile");

my $qc_status_file = "MetaAmp.qcstatus.txt";
open STATUS, '>', $qc_status_file or err("Can't open $qc_status_file");

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


sub reverse_complement_IUPAC {
        my ($dna) = @_;
	
        # reverse the DNA sequence
        my $revcomp = reverse($dna);

        # complement the reversed DNA sequence
        $revcomp =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
        return $revcomp;
}





1;

__END__

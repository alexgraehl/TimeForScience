#!/usr/bin/perl -w

use strict; use warnings; use diagnostics;
use File::Temp;
use Getopt::Long;

my $CMD_SUFFIX = ".cmd.tmp.sh";

sub agw_make_temp() {
	my $timeInSecondsSince1970 = time();
	my ($temp_filehandle, $temp_filename) = File::Temp::tempfile(TEMPLATE=>"cluster.${timeInSecondsSince1970}.XXXXX", DIR=>"", SUFFIX=>"$CMD_SUFFIX", UNLINK=>0); return($temp_filename);
}

if (scalar(@ARGV) == 0) { die "[ERROR]: Cannot execute with NO commands... try some example like 'qplz echo \"Hello\"'...\n"; }

my @args = @ARGV;

my $submission_looks_like_a_pbs_script = 0;

# Short = 30 minutes
# Long = 14 days

# $ qsub -q General -W group_list="genqueue" [/path/to/script]
# $ qsub -q Bio -W group_list="bioqueue" [/path/to/script]

#my $filename  = agw_make_temp();
#my $jobname   = $filename; $jobname =~ s/${CMD_SUFFIX}$//; # filename WITHOUT the annoying suffix
#
#print("Filename is: $filename\nJobname is: $jobname\n");

my $ncpus      = 1;
my $memMB     = 150;

my $memString  = $memMB . "mb";
my $timeString = "00:29:00";

#my $pbsheader = qq{}
#  . qq{#!/usr/bin/bash} . "\n"
#  . qq{#PBS -N $jobname} . "\n"
#  . qq{#PBS -l ncpus=${ncpus}} . "\n"
#  . qq{#PBS -l mem=${memMB}mb} . "\n"
#  . qq{#PBS -l walltime=${timeInMin}:00} . "\n";
#
#open(my $fh, '>', $filename) or die "Could not open file '$filename' $!";
#print $fh $pbsheader . "\n";
#print $fh join(" ", @args) . "\n";
#close($fh);

my $qsub_exe="qsub";
my $QGRP="Bio";
my $QQUE="bioqueue";
my $cmd="";

my $cmdargs = join(" ", @args);

my $qsub_common = qq{$qsub_exe -q "$QGRP" -W group_list="$QQUE" -l ncpus=$ncpus -l mem=$memString -l walltime=$timeString};

if ($submission_looks_like_a_pbs_script) {
	# looks like someone passed in a script...
	$cmd = qq{$qsub_common $cmdargs};
} else {

	$cmd = qq{echo '$cmdargs' | $qsub_common};
}

print("Executing: $cmd\n");

my $exCode = system($cmd);

print STDERR qq{[QUEUED] Submitted this job to the queue: $cmd\n};
if ($exCode == 0) {
	print STDERR "[QUEUED] You can type 'qstat' (basic) or 'qstats' (fancy) to check on the status of your job. --Alex\n";
} else {
	die qq{[ERROR] Curses! Something went wrong, and the queue command returned the error code '$exCode'. Unclear what this means, but something is probably wrong with your input command.\n};
}

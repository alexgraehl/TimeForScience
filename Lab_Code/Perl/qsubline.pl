#!/usr/bin/perl -w

use strict; use warnings; use diagnostics;
use File::Temp;
use Getopt::Long;

my $CMD_SUFFIX = ".cmd.tmp.sh";

sub agw_make_temp() {
	my $timeInSecondsSince1970 = time();
	my ($temp_filehandle, $temp_filename) = File::Temp::tempfile(TEMPLATE=>"cluster.${timeInSecondsSince1970}.XXXXX", DIR=>"", SUFFIX=>"$CMD_SUFFIX", UNLINK=>0); return($temp_filename);
}

if (scalar(@ARGV) == 0) { die "Cannot execute with NO commands...\n"; }

my @args = @ARGV;

# Short = 30 minutes
# Long = 14 days

# $ qsub -q General -W group_list="genqueue" [/path/to/script]
# $ qsub -q Bio -W group_list="bioqueue" [/path/to/script]

my $filename  = agw_make_temp();
my $jobname   = $filename; $jobname =~ s/${CMD_SUFFIX}$//; # filename WITHOUT the annoying suffix

print("Filename is: $filename\nJobname is: $jobname\n");

my $ncpu      = 1;
my $memMB     = 150;
my $timeInMin = 29;

my $pbsheader = qq{}
  . qq{#!/usr/bin/bash} . "\n"
  . qq{#PBS -N $jobname} . "\n"
  . qq{#PBS -l ncpus=${ncpu}} . "\n"
  . qq{#PBS -l mem=${memMB}mb} . "\n"
  . qq{#PBS -l walltime=${timeInMin}:00} . "\n";

open(my $fh, '>', $filename) or die "Could not open file '$filename' $!";
print $fh $pbsheader . "\n";
print $fh join(" ", @args) . "\n";
close($fh);

my $QGRP="Bio";
my $QQUE="bioqueue";
my $cmd=qq{qsub -q "$QGRP" -W group_list="$QQUE" $filename};

print("Executing: $cmd\n");
system($cmd);
print("[DONE] submitting that job to the queue... now you should wait for it to finish running!\n");

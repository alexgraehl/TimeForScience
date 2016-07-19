#!/usr/bin/perl
use strict;  use warnings;  use diagnostics;
use POSIX;
use List::Util qw(max min);
use Getopt::Long;
use Carp; # backtrace on errors. Has the "confess" function. Use this instead of "die" if you want useful information!   

#use File::Basename;

$| = 1; # Always flush text output IMMEDIATELY to the console, don't wait to buffer terminal output! Setting this to zero can cause STDERR and STDOUT to be interleaved in weird ways.

#no warnings 'numeric';
#use Scalar::Util;
#print Scalar::Util::looks_like_number($string), "\n";

my $UNIX_BIOGRP = "bioqueue";
my $UNIX_GENGRP = "genqueue";
my %QSETTINGS = ( "unix_gname_to_gid" => {"$UNIX_BIOGRP"   => "35098"
					  , "$UNIX_GENGRP" => "35099" }
		  , "dest"             => {"$UNIX_BIOGRP"   => "-q Bio"
					   , "$UNIX_GENGRP" => "-q General" }
		  , "grouplist"        => {"$UNIX_BIOGRP"   => "-W group_list=bioqueue"
					   , "$UNIX_GENGRP" => "-W group_list=genqueue" }
		  , "max_ncpus"        => {"$UNIX_BIOGRP"   => 16
					   , "$UNIX_GENGRP" => 16 }
		  , "max_mem"       => {"$UNIX_BIOGRP"   => 256 # GB of ram!
					   , "$UNIX_GENGRP" => 64 }
		  , "default_ncpus"    => {"$UNIX_BIOGRP"   => 1
					   , "$UNIX_GENGRP" => 1 }
		  , "default_mem"      => {"$UNIX_BIOGRP"   => 4
					   , "$UNIX_GENGRP" => 4 }
		  , "default_walltime"    => {"$UNIX_BIOGRP"   => "00:29:00"
					   , "$UNIX_GENGRP" => "00:29:00" }
	 );


sub tryToLoadModule($) {
	my $x = eval("require $_[0]");
	if ((defined($@) && $@)) {
		warn "We FAILED to load module $_[0]. Skipping this module, but continuing with the program.";
		return 0;	# FAILURE
	} else {
		$_[0]->import();
		return 1;	# SUCCESS
	}
}

my $SHOULD_USE_COLORS = tryToLoadModule("Term::ANSIColor");
if ($SHOULD_USE_COLORS) {
	use Term::ANSIColor;
}

sub warnPrint($) { chomp($_[0]); warn(safeColor("[WARNING]: " . $_[0] . "", "yellow on_black")); } # regarding "warn": if the string ends with a newline it WON'T print the line number!

sub safeColor($;$) {		# one required and one optional argument
	## Returns colored text, but only if $SHOULD_USE_COLORS is set.
	## Allows you to totally disable colored printing by just changing $SHOULD_USE_COLORS to 0 at the top of this file
	# Colorstring is OPTIONAL, and can be something like "red on_blue" or "red" or "magenta on_green"
	# Example usage:
	#    *    print STDERR safeColor("This warning message is red on yellow", "red on_yellow");
	my ($message, $color) = @_;
	return (($SHOULD_USE_COLORS && defined($color)) ? (Term::ANSIColor::colored($message, $color) . Term::ANSIColor::color("reset")) : $message);
}

sub printColorStderr($;$) {
	# prints color to STDERR *UNLESS* it is re-directed to a file, in which case NO COLOR IS PRINTED.
	my ($msg, $col) = @_; # Only prints in color if STDERR is to a terminal, NOT if it is redirected to an output file!
	if (! -t STDERR) {
		$col = undef;
	}			# no coloration if this isn't to a terminal
	print STDERR safeColor($msg, $col);
}

sub printColorStdout($;$) {
	# prints color to STDOUT *UNLESS* it is re-directed to a file, in which case NO COLOR IS PRINTED.
	my ($msg, $col) = @_; # Only prints in color if STDOUT is to a terminal, NOT if it is redirected to an output file!
	if (! -t STDOUT) {
		$col = undef;
	}			# no coloration if this isn't to a terminal
	print STDOUT safeColor($msg, $col);
}

sub dryNotify(;$) {		# one optional argument
	my ($msg) = @_;
	$msg = (defined($msg)) ? $msg : "This was only a dry run, so we skipped executing a command.";
	print STDERR safeColor("[DRY RUN]: $msg\n", "black on_yellow");
}

sub notify($) {			# one required argument
	my ($msg) = @_;
	warn safeColor("[DRY RUN]: $msg\n", "cyan on_blue");
}

sub main();
sub quitWithUsageError($) { print("[ERROR]: " . $_[0] . "\n"); printUsageAndQuit(); print($_[0] . "\n"); }
sub printUsageAndQuit() { printUsage(); exit(1); }

sub printUsage() {
	print STDOUT <DATA>;
	exit(0);
}

sub getOurQueueGroup() {
	# == See if we are in the bio group
	my @gids = POSIX::getgroups(); # getgroups() is from the POSIX module
	my @BIOQUEUE_GROUP_IDS = ( $QSETTINGS{"unix_gname_to_gid"}{"$UNIX_BIOGRP"} );
	for my $priority (@BIOQUEUE_GROUP_IDS) {
		if (grep(/^$priority$/, @gids)) { return $UNIX_BIOGRP; } # apparently the user belongs to a privileged group--let them use more CPUs, etc.
	}
	return $UNIX_GENGRP; # otherwise...
}

# ==1==
sub main() {			# Main program
	my ($delim) = "\t";
	my ($decimalPlaces) = 4; # How many decimal places to print, by default
	$Getopt::Long::passthrough = 1; # ignore arguments we don't recognize in GetOptions, and put them in @ARGV

	my ($pbs_ncpus, $pbs_mem, $pbs_walltime) = (undef, undef, undef);

	GetOptions("help|?|man" => sub { printUsageAndQuit(); }
		   , "delim|d=s" => \$delim
		   , "ncpus|c=i" => \$pbs_ncpus
		   , "mem|m=i" => \$pbs_mem
		   , "walltime|t=s" => \$pbs_walltime
		  ) or printUsageAndQuit();

	my $grp = getOurQueueGroup();

	if (defined($pbs_mem)) {
		($pbs_mem =~ m/^\d+(gb|g|)$/i) or quitWithUsageError("(Bad value to --mem / -m): Your memory request (in GIGABYTES) was invalid. You need to speciy an integer number of GB (e.g. '10' or '10gb' or '10g'. You specified this value: $pbs_mem");
		$pbs_mem =~ s/[A-Za-z]//g; # remove any letters from it, now $pbs_mem is purely numeric
		($pbs_mem =~ m/^\d+$/) or confess("Programming error: somehow failed to make pbs_mem numeric! Offending variable was: $pbs_mem");; # remove any letters from it, now $pbs_mem is purely numeric
		($pbs_mem <= $QSETTINGS{max_mem}{$grp}) or quitWithUsageError("(Bad value to --mem / -m): You requested too much RAM! The maximum you can specify at this specific queue/user/group combination is " . $QSETTINGS{max_mem){$grp} . " (in gigabytes), but your request was for this number: $pbs_mem");
	} else {
		$pbs_ncpus = $QSETTINGS{default_mem}{$grp};
	}

	if (defined($pbs_ncpus)) {
		($pbs_ncpus =~ m/[1-9]\d+/) or quitWithUsageError("(Bad value to --ncpus / -c): You need to specify a non-zero integer number of CPU cores to use. You specified this value: $pbs_ncpus");
		($pbs_ncpus <= $QSETTINGS{max_ncpus){$grp}) or quitWithUsageError("(Bad value to --ncpus / -c): You specified TOO MANY cpus. The maximum you can specify at this specific queue/user/group combination is " . $QSETTINGS{max_ncpus){$grp} . ", but your request was for this number: $pbs_ncpus");
	} else {
		$pbs_ncpus = $QSETTINGS{default_ncpus}{$grp};
	}

	if (defined($pbs_walltime)) {
		($pbs_walltime =~ m/\d+:\d\d:\d\d/) or quitWithUsageError("(Bad value to --walltime / -t): You need to specify a valid walltime in this format: 11:22:33 (hours, minutes, seconds). You specified this value: $pbs_walltime");
	} else {
		$pbs_walltime = $QSETTINGS{default_walltime}{$grp};
	}


	my $numUnprocessedArgs = scalar(@ARGV);
	if ($numUnprocessedArgs == 0) {
		quitWithUsageError("Cannot execute with NO commands... try some example like 'qplz echo \"Hello\"'...\n");
	}

	my $filename1 = undef;
	foreach (@ARGV) { # these were arguments that were not understood by GetOptions
		if (!defined($filename1)) {
			$filename1 = $_;
		} else {
			print STDERR "Unprocessed argument: $_\n";
		}
	}

	print STDERR "===============================\n";
	print STDOUT safeColor("test color\n", "blue on_yellow");
	printColorStdout("hey", "red");
	printColorStdout("what", "red");
	printColorStderr("stderr what", "blue on_white");
	printColorStdout("is", "red on_blue");
	printColorStdout("this", "red on_green");
	printColorStderr("stderr what", "red on_white");
	print STDERR "hellow\n";
	print STDERR "Test color!\n";
	print STDERR "===============================\n";

	my $qsub_exe = "qsub";

	#sub agw_make_temp() {
	#	my $timeInSecondsSince1970 = time();
	#	my ($temp_filehandle, $temp_filename) = File::Temp::tempfile(TEMPLATE=>"cluster.${timeInSecondsSince1970}.XXXXX", DIR=>"", SUFFIX=>"$CMD_SUFFIX", UNLINK=>0); return($temp_filename);
	#}

	my $submission_looks_like_a_pbs_script = 0;

	# $ qsub -q General -W group_list="genqueue" [/path/to/script]
	# $ qsub -q Bio -W group_list="bioqueue" [/path/to/script]

	#my $filename  = agw_make_temp();
	#my $jobname   = $filename; $jobname =~ s/${CMD_SUFFIX}$//; # filename WITHOUT the annoying suffix
	#
	#print("Filename is: $filename\nJobname is: $jobname\n");

	my $timeString = "00:29:00";

	my $qdest      = $QSETTINGS{dest}{$grp};
	my $qgrouplist = $QSETTINGS{grouplist}{$grp};
	my $cmd="";

	my $cmdargs = join(" ", @ARGV);

	my $qsub_common = qq{$qsub_exe -q "$qdest" -W group_list="$qgrouplist" -l ncpus=${pbs_ncpus} -l mem=\"${pbs_mem}gb\" -l walltime=\"${pbs_walltime}\"};

	#my $mem      = "-l mem="      . (defined($hr->{pbs_mem})      ? $hr->{pbs_mem}      : "8gb");
	#my $walltime = "-l walltime=" . (defined($hr->{pbs_walltime}) ? $hr->{pbs_walltime} : "299:59:59");
	#my $ncpus    = "-l ncpus="    . (defined($hr->{pbs_ncpus})    ? $hr->{pbs_ncpus}    : "1");
	#my $nodes    = "-l nodes=1:ppn=1"; # <-- ????????? unclear if this is useful
	my $stderr   = "";	#"-e /dev/null"
	my $stdout   = "";	#"-o /dev/null"



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

	print STDERR "Ok, now you should run 'qstats' and look for your output in these STDERR / STDOUT files...\n";

	print STDERR safeColor("===============================\n", "green");
} # end main()


main();


END {
	# Runs after everything else.
	# Makes sure that the terminal text is back to its normal color.
	if ($SHOULD_USE_COLORS) {
		print STDERR Term::ANSIColor::color("reset"); print STDOUT Term::ANSIColor::color("reset");
	}
}

exit(0);
# ====

__DATA__

MYPROGRAM.pl  [OPTIONS]

  by YOUR NAME, THE_YEAR

  THIS PROGRAM DOES SOMETHING. YOU GIVE IT ONE THING AND THIS OTHER THING,
  AND THE OUTPUT IS THIS FINAL THING.

  See the examples below for more information.

 CAVEATS:

MAYBE IT TAKES 30 MINUTES TO RUN.

 OPTIONS:

--delim = DELIMITER   (Default: tab)
  Sets the input delimiter to DELIMITER.

 EXAMPLES:

MYPROGRAM.pl --help
  Displays this help

  MYPROGRAM.pl  --works=yes --bugs=4  -q
  Does nothing. -q indicates "quiet" operation.


KNOWN BUGS:

None known.

TO DO:

Add ???.

  --------------




#!/usr/bin/perl
use strict;  use warnings;  #use diagnostics;
use POSIX;
use List::Util qw(max min);
use Getopt::Long;
use Carp; # backtrace on errors. Has the "confess" function. Use this instead of "die" if you want useful information!   

#use File::Basename;

$| = 1; # Always flush text output IMMEDIATELY to the console, don't wait to buffer terminal output! Setting this to zero can cause STDERR and STDOUT to be interleaved in weird ways.

#no warnings 'numeric';
#use Scalar::Util;
#print Scalar::Util::looks_like_number($string), "\n";

my $QSUB_EXE = "qsub";
my $UNIX_BIOGRP = "bioqueue";
my $UNIX_GENGRP = "genqueue";
my %QSETTINGS = ( "unix_gname_to_gid" => {"$UNIX_BIOGRP"   => "35098"
					  , "$UNIX_GENGRP" => "35099" }
		  , "dest"             => {"$UNIX_BIOGRP"   => "Bio"
					   , "$UNIX_GENGRP" => "General" }
		  , "grouplist"        => {"$UNIX_BIOGRP"   => "bioqueue" # <-- could theoretically be different from the UNIX_BIOGRP, so don't change this from plain text
					   , "$UNIX_GENGRP" => "genqueue" } # <-- could theoretically be different from the UNIX_GENGRP
		  , "max_ncpus"        => {"$UNIX_BIOGRP"   => 16
					   , "$UNIX_GENGRP" => 16 }
		  , "max_mem"          => {"$UNIX_BIOGRP"   => 256 # GB of ram!
					   , "$UNIX_GENGRP" => 64 }
		  , "max_walltime_hours" => {"$UNIX_BIOGRP"   => 335
					   , "$UNIX_GENGRP" => 335 }
		  , "default_ncpus"    => {"$UNIX_BIOGRP"   => 1
					   , "$UNIX_GENGRP" => 1 }
		  , "default_mem"      => {"$UNIX_BIOGRP"   => 4
					   , "$UNIX_GENGRP" => 4 }
		  , "default_walltime" => {"$UNIX_BIOGRP"   => "00:29:00"
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
	my ($msg) = @_; chomp($msg);
	$msg = (defined($msg)) ? $msg : "This was only a dry run, so we skipped executing a command.";
	print STDERR safeColor("[DRY RUN]: $msg\n", "black on_yellow");
}

sub notify($) {			# one required argument
	my ($msg) = @_; chomp($msg);
	warn safeColor("[DRY RUN]: $msg\n", "cyan on_blue");
}


sub printCool($) {			# one required argument
	my ($msg) = @_; chomp($msg);
	print STDERR safeColor("[PROGRESS REPORT]: $msg\n", "white on_blue");
}

sub main();
sub quitWithUsageError($) { print("[ERROR]: " . $_[0] . "\n"); printUsageAndQuit(); print($_[0] . "\n"); }

my $GLOBAL_WARN_STRING = "";
sub ourWarn($) { my $s = $_[0]; chomp($s); print("[WARNING]: $s\n"); $GLOBAL_WARN_STRING .= "$s\n"; }

sub printUsageAndQuit() { printUsage(); exit(1); }

sub printUsage() {
	print STDOUT <DATA>;
	exit(0);
}




sub fileIsProbablySomeScript($) { # detect if a filename seems to be an ok-to-submit PBS script
	my ($filename) = @_; # filename
	if (-e $filename) {
		# it's a file that EXISTS at this path, then maybe it's a script we can run??
		return 1;
	}
	#if ($maybeScriptToSubmit =~ m/[.](sh|pl|py|R|rb)$/i) {
		# it has a common script extension like ".sh" or ".pl"
	#return 1;
	#}
	return 0;
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
sub main() { # Main program
	my ($delim) = "\t";
	my ($decimalPlaces) = 4; # How many decimal places to print, by default
	$Getopt::Long::passthrough = 1; # ignore arguments we don't recognize in GetOptions, and put them in @ARGV

	my ($pbs_ncpus, $pbs_mem, $pbs_walltime) = (undef, undef, undef);
	my ($pbs_submit_file) = undef;


	GetOptions("help|?|man" => sub { printUsageAndQuit(); }
		   , "delim|d=s" => \$delim
		   , "ncpus|c=i" => \$pbs_ncpus
		   , "mem|m=i" => \$pbs_mem
		   , "walltime|t=s" => \$pbs_walltime
		   , "f=s"          => \$pbs_submit_file
		  ) or printUsageAndQuit();

	my $grp = getOurQueueGroup();

	if (defined($pbs_mem)) {
		($pbs_mem =~ m/^\d+(gb|g|)$/i) or quitWithUsageError("(Bad value to --mem / -m): Your memory request (in GIGABYTES) was invalid. You need to speciy an integer number of GB (e.g. '10' or '10gb' or '10g'. You specified this value: $pbs_mem");
		$pbs_mem =~ s/[A-Za-z]//g; # <-- remove any letters from it, now $pbs_mem is purely numeric (i.e., "24gb => 24")
		($pbs_mem =~ m/^\d+$/) or confess("Programming error: somehow failed to make pbs_mem numeric! Offending variable was: $pbs_mem"); # remove any letters from it, now $pbs_mem is purely numeric
		($pbs_mem <= $QSETTINGS{max_mem}{$grp}) or quitWithUsageError("(Bad value to --mem / -m): You requested too much RAM! The maximum you can specify at this specific queue/user/group combination is " . $QSETTINGS{max_mem}{$grp} . " (in gigabytes), but your request was for this number: $pbs_mem");
	} else {
		$pbs_mem = $QSETTINGS{default_mem}{$grp};
	}

	if (defined($pbs_ncpus)) {
		($pbs_ncpus =~ m/[1-9]\d*/) or quitWithUsageError("(Bad value to --ncpus / -c): You need to specify a non-zero integer number of CPU cores to use. You specified this value: $pbs_ncpus");
		($pbs_ncpus <= $QSETTINGS{max_ncpus}{$grp}) or quitWithUsageError("(Bad value to --ncpus / -c): You specified TOO MANY cpus. The maximum you can specify at this specific queue/user/group combination is " . $QSETTINGS{max_ncpus}{$grp} . ", but your request was for this number: $pbs_ncpus");
	} else {
		$pbs_ncpus = $QSETTINGS{default_ncpus}{$grp};
	}

	if (defined($pbs_walltime)) {
		($pbs_walltime =~ m/\d+:\d\d:\d\d/) or quitWithUsageError("(Bad value to --walltime / -t): You need to specify a valid walltime in this format: 11:22:33 (hours, minutes, seconds). You specified this value: $pbs_walltime");
	} else {
		$pbs_walltime = $QSETTINGS{default_walltime}{$grp};
	}


	my ($pbs_wall_hr, $pbs_wall_min, $pbs_wall_sec) = split(':', $pbs_walltime);
	
	($pbs_wall_hr <= $QSETTINGS{max_walltime_hours}{$grp}) or quitWithUsageError("You requested TOO MUCH walltime! The maximum is: $QSETTINGS{max_walltime_hours}{$grp}:00:00 . Try again with a smaller value!");

	my $numUnprocessedArgs = scalar(@ARGV);
	if (0 == $numUnprocessedArgs) {
		quitWithUsageError("Cannot execute with NO commands... try some example like 'qplz echo \"Hello\"'...\n");
	}

	if (1 == $numUnprocessedArgs and fileIsProbablySomeScript($ARGV[0])) {
		# See if the unprocessed argument is a FILENAME (a script to submit)
		$pbs_submit_file = $ARGV[0];
	}



	my $qdest      = $QSETTINGS{dest}{$grp};
	my $qgrouplist = $QSETTINGS{grouplist}{$grp};

	my $stderr   = "";	#"-e /dev/null"
	my $stdout   = "";	#"-o /dev/null"

	my $qsub_common = qq{$QSUB_EXE }
	  . qq{ -q "$qdest" }
	  . qq{ -W group_list="$qgrouplist" }
	  . qq{ -l ncpus=${pbs_ncpus} }
	  . qq{ -l mem="${pbs_mem}gb" }
	  . qq{ -l walltime="${pbs_walltime}"};

	my $cmd = undef;
	if (defined($pbs_submit_file)) {
		(-e $pbs_submit_file) or quitWithUsageError("It looked like you submitted a script file directly on the command line (we though this was a script to submit to PBS: \"$pbs_submit_file\"), but somehow it seems like that file did not exist. Weird!");
		if ($pbs_submit_file !~ m/[.](sh|pl|py|R|rb)$/) {
			ourWarn(qq{You submitted a file to qsub (specifically, "$pbs_submit_file") that did not have a common script ending... just be aware of this!});
		}
		$cmd = qq{$qsub_common $pbs_submit_file};
	} else {
		# ok, looks like we are just submitting a QUICK job right on the command line
		my $cmdArgs = join(" ", @ARGV); # mash the command line arguments together
		$cmd = qq{echo '$cmdArgs' | $qsub_common};
	}

	printCool(qq{Submitting this job to the queue:\n    YOUR JOB -->   $cmd\n});
	my $exitText = `$cmd`;
	my $exitCode = $?; # <-- would be the result of system($cmd)
	if ($exitCode == 0) {
		printCool("You can type 'qstat' (basic) or 'qstats' (fancy) to check on the status of your job. --Alex\n");
	} else {
		die qq{[ERROR] Curses! Something went wrong, and the queue command returned the error code '$exitCode'. Unclear what this means, but something is probably wrong with your input command.\n};
	}

	#my $filename1 = undef;
	#foreach (@ARGV) { # these were arguments that were not understood by GetOptions
	#	if (!defined($filename1)) {
	#		$filename1 = $_;
	#	} else {
	#		print STDERR "Unprocessed argument: $_\n";
	#	}
	#}

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

	#sub agw_make_temp() {
	#	my $timeInSecondsSince1970 = time();
	#	my ($temp_filehandle, $temp_filename) = File::Temp::tempfile(TEMPLATE=>"cluster.${timeInSecondsSince1970}.XXXXX", DIR=>"", SUFFIX=>"$CMD_SUFFIX", UNLINK=>0); return($temp_filename);
	#}

	# $ qsub -q General -W group_list="genqueue" [/path/to/script]
	# $ qsub -q Bio -W group_list="bioqueue" [/path/to/script]

	#my $filename  = agw_make_temp();
	#my $jobname   = $filename; $jobname =~ s/${CMD_SUFFIX}$//; # filename WITHOUT the annoying suffix
	#
	#print("Filename is: $filename\nJobname is: $jobname\n");
	

	print STDERR "Your job will be allowed to run for this long: ${pbs_wall_hr} hours and ${pbs_wall_min} minutes\n";
	print STDERR "Your job has been allocated the following:\n";
	print STDERR "                CPU CORES: ${pbs_ncpus}\n";
	print STDERR "                      RAM: ${pbs_mem}gb\n";
	print STDERR "                     TIME: ${pbs_walltime}\n";

	print STDERR "Now you can type these commands to check your job:\n";
	print STDERR "             qstats   (print out a color list of jobs that are running\n";
	print STDERR "             qstat    (print out monochrome list of the above jobs\n";
	print STDERR "             qstat -u \$USER (print out YOUR jobs only\n";
	#print STDERR "Ok, now you should run 'qstats' and look for your output in these STDERR / STDOUT files...\n";

	#print STDERR $GLOBAL_WARN_STRING . "\n";
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




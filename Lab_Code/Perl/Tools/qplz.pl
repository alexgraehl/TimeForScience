#!/usr/bin/perl
use strict;  use warnings;  #use diagnostics;
use POSIX;
use List::Util qw(max min);
use Getopt::Long;
use File::Basename;
use Carp; # backtrace on errors. Has the "confess" function. Use this instead of "die" if you want useful information!   

#use File::Basename;

$| = 1; # Always flush text output IMMEDIATELY to the console, don't wait to buffer terminal output! Setting this to zero can cause STDERR and STDOUT to be interleaved in weird ways.

#no warnings 'numeric';
#use Scalar::Util;
#print Scalar::Util::looks_like_number($string), "\n";

my $JOBNAME_MAX_LEN = 32; # cut it down to some reasonable size

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
		  , "default_walltime" => {"$UNIX_BIOGRP"   => "23:59:59"
					   , "$UNIX_GENGRP" => "23:59:59" }
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

sub printBadNews($) {			# one required argument
	my ($msg) = @_; chomp($msg);
	print STDERR safeColor("[ERROR]: $msg\n", "white on_red");
}
sub printBadNewsAndDie($) { printBadNews($_[0]); die $_[0]; }

sub main();

sub quitWithUsageError($) { printBadNews($_[0]); printUsage(); printBadNews($_[0]); exit(1); }

my $GLOBAL_WARN_STRING = "";
sub ourWarn($) { my $s = $_[0]; chomp($s); print("[WARNING]: $s\n"); $GLOBAL_WARN_STRING .= "$s\n"; }

sub printUsageAndQuit() { printUsage(); exit(1); }

sub printUsage() {
	print STDOUT <DATA>;
}


sub regarg($) {
	my ($filename) = @_;
	open my $fff, '<', $filename or die "whoops, no file somehow";
	while (my $line = <$fff>) {
		if ($line =~ /REGEX/) {
			
		}
	}
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
	GetOptions("h|help|?|man" => sub { printUsageAndQuit(); }
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

	if (not defined($pbs_walltime)) {
		$pbs_walltime = $QSETTINGS{default_walltime}{$grp};
	} elsif ($pbs_walltime =~ m/^\d+$/) { # if it's literally JUST a single number (the number of hours)
		$pbs_walltime = "${pbs_walltime}:00:00"; # I guess it's JUST the number of hours
	}
	
	($pbs_walltime =~ m/^(\d+):(\d\d):(\d\d)$/) or quitWithUsageError("(Bad value to --walltime / -t): You need to specify a valid walltime in this format: 11:22:33 (hours, minutes, seconds). You specified this value: $pbs_walltime");
	my ($pbs_wall_hr, $pbs_wall_min, $pbs_wall_sec) = ($1, $2, $3); # <-- results from the match expression above
	($pbs_wall_hr <= $QSETTINGS{max_walltime_hours}{$grp}) or quitWithUsageError("You requested TOO MUCH walltime! The maximum is: $QSETTINGS{max_walltime_hours}{$grp}:00:00 . Try again with a smaller value!");

	my $numUnprocessedArgs = scalar(@ARGV);
	if (defined($pbs_submit_file) and $numUnprocessedArgs > 0) {
		quitWithUsageError("Weird, you specified a pbs_submit_file on the command line with the '-f' option, yet there is still an unrecognized argument at the end (which normally would be a script file, maybe)). We have some unrecognized arguments on the command line. See why these are here, as qplz does not understand them: " . join(" ", @ARGV)); # ok, we should NOT have any unprocessed arguments if the user defined a script file
	}

	if (!defined($pbs_submit_file) and 0 == $numUnprocessedArgs) {
		quitWithUsageError("Cannot submit a job with no job script or command! Make sure you actually submitted a script filename (or arguments on the command line, for example:  qplz echo \"Hello\")...\n");
	}



	if (1 == $numUnprocessedArgs and fileIsProbablySomeScript($ARGV[0])) {
		# See if the unprocessed argument is a FILENAME (a script to submit)
		$pbs_submit_file = $ARGV[0];
	} else {
		# looks like the user submitted a 'quick command' on the command line, like 'qplz.pl pwd'
	}

	my $qdest      = $QSETTINGS{dest}{$grp};
	my $qgrouplist = $QSETTINGS{grouplist}{$grp};
	my $stderr   = "";	#"-e /dev/null"
	my $stdout   = "";	#"-o /dev/null"
	my $jobname  = undef;
	if (defined($pbs_submit_file)) {
		(-e $pbs_submit_file) or quitWithUsageError("It looked like you submitted a script file directly on the command line (we though this was a script to submit to PBS: \"$pbs_submit_file\"), but somehow it seems like that file did not exist. Weird!");
		if ($pbs_submit_file !~ m/[.](sh|pl|py|R|rb)$/) {
			ourWarn(qq{You submitted a file to qsub (specifically, "$pbs_submit_file") that did not have a common script ending... just be aware of this!});
		}
		$jobname = File::Basename::basename($pbs_submit_file); # Job name is the SUBMITTED SCRIPT name (but not the full path)
	} else {
		$jobname = join("_", @ARGV);
	}
	$jobname =~ s/[\W\s]/_/g; # replace any "weird" non-word characters with underscores
	substr($jobname, $JOBNAME_MAX_LEN); # Job name is the first N characters of the command line arguments.

	my $qsub_common = qq{$QSUB_EXE }
	  . qq{ -V }
	  . qq{ -N "$jobname" }
	  . qq{ -q "$qdest" }
	  . qq{ -W group_list="$qgrouplist" }
	  . qq{ -l ncpus=${pbs_ncpus} }
	  . qq{ -l mem="${pbs_mem}gb" }
	  . qq{ -l walltime="${pbs_walltime}"};

	my $pwd = `pwd`; chomp($pwd);

	my $cmd = undef;
	if (defined($pbs_submit_file)) {
		$cmd = qq{$qsub_common $pbs_submit_file};
	} else { # ok, looks like we are just submitting a QUICK job right on the command line
		my $cmdArgs = join(" ", @ARGV); # mash the command line arguments together
		$cmd = qq{echo 'cd "$pwd" && $cmdArgs' | $qsub_common};
	}
	printCool(qq{Submitting this actual command to the queue:\n    ACTUAL JOB SUBMISSION COMMAND -->   $cmd\n});
	my $exitText = `$cmd`;
	my $exitCode = $?; # <-- would be the result of system($cmd)
	if ($exitCode == 0) {
		printCool("You can type 'qstat' (basic) or 'qstats' (fancy) to check on the status of your job. --Alex\n");
	} else {
		printBadNewsAndDie("[ERROR] Curses! Something went wrong, and the queue command returned the error code '$exitCode'. Unclear what this means, but something is probably wrong with your input command.\n");
	}

	#my $filename1 = undef;
	#foreach (@ARGV) { # these were arguments that were not understood by GetOptions
	#	if (!defined($filename1)) {
	#		$filename1 = $_;
	#	} else {
	#		print STDERR "Unprocessed argument: $_\n";
	#	}
	#}

	# print STDERR "===============================\n";
	# print STDOUT safeColor("test color\n", "blue on_yellow");
	# printColorStdout("hey", "red");
	# printColorStdout("what", "red");
	# printColorStderr("stderr what", "blue on_white");
	# printColorStdout("is", "red on_blue");
	# printColorStdout("this", "red on_green");
	# printColorStderr("stderr what", "red on_white");
	# print STDERR "hellow\n";
	# print STDERR "Test color!\n";
	# print STDERR "===============================\n";

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

	print STDERR "Your job has been allocated the following:\n";
	print STDERR "                CPU CORES: ${pbs_ncpus}\n";
	print STDERR "                      RAM: ${pbs_mem}gb\n";
	print STDERR "                     TIME: ${pbs_walltime}\n";

	print STDERR "Now you can type these commands to check your job:\n";
	print STDERR "To check your job 1:  qstats   (print out a color list of jobs that are running\n";
	print STDERR "To check your job 2:  qstat    (print out monochrome list of the above jobs\n";
	print STDERR "To check your job 3:  qstat -u \$USER (print out YOUR jobs only\n";
	print STDERR "If you want to cancel your job (maybe you just realized that it needs more time / RAM):\n";
	print STDERR "    To delete a job: 1. Find the 'Job id' number with 'qstat' (leftmost column)\n";
	print STDERR "    To delete a job: 2. Then use 'qdel ####' (that same number) to cancel it\n";
	print STDERR "    To delete all your jobs ever (dangerous!): qselect -u \$USER | xargs qdel  <-- deletes all your jobs\n";
	#print STDERR "Ok, now you should run 'qstats' and look for your output in these STDERR / STDOUT files...\n";

	printColorStderr("Your job will be allowed to use ${pbs_mem} GB of RAM and run for ${pbs_wall_hr} hours and ${pbs_wall_min} minutes before it is cancelled.\n", "white on_red");

	if ($pbs_wall_hr <= 0) {
		printColorStderr("WARNING: Note that your job will be cancelled if it takes longer than ${pbs_wall_min} minutes to run!\n", "white on_red");
	}

	# maybe make and then delete a fake file here just to refresh the queue

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

qplz.pl [OPTIONS]  <script_or_command>
  by Alex Williams, 2016
  Tested with PBS Pro v 13.

"QUEUE PLEASE" (qplz) is a frontend to 'qsub' that makes it easier to submit queue jobs to PBS Pro 13.

Examples: qplz.pl "pwd"   or  qplz.pl my_script.sh

OPTIONS:

--walltime=HH:MM:SS or "-t HH:MM:SS"
 Default: 00:30:00 = 30 minutes
  Lets your job run for at least this long. It gets auto-killed if this time is exceeded.
  Example of a job that wants to run for 123 hours and 45 minutes:  -t 123:45:00

--mem=INTEGER or "-m INTEGER"
 Default: 4 (4 gigabytes)  (Example ways to request 4 GB of RAM: "-m 4" or "-m 4gb" or "-m 4g")
  Lets your job use this much memory, in gigabytes. Job will be auto-killed if it tries to use more.

--ncpus=INTEGER or "-c INTEGER"
 Default: 1 (1 core)
  Lets your job use this many cores.
  A job can *try* to use more, but all of its threads will be placed onto this many cores.

  (Note: If hyperthreading is enabled (which it typically is NOT on our server), then this number
  would be the number of *hyperthreaded* cores instead of physical cores.)

-f FILENAME (or just put the filename at the end of the line)
  Submit the commands in this filename. Example:   qplz.pl -f myscript.sh  -t 1:00:00
  Note that this is normally equivalent to just putting the filename at the end of the line, too, like:
     * qplz.pl myscript.sh

TO DO: add '-o' (STDOUT) and '-e' (STDERR) options

EXAMPLES:

qplz.pl ls
  Lists your home directory

qplz.pl -t 1:00:00 pwd
or
qplz.pl -t 1 pwd
  Prints the current directory, allowing the script one hour to do so.

qplz.pl -t 24 -m 8 -c 2 myscript.pl
  Runs 'myscript.pl' for up to 24 hours and 8 GB of RAM, with 2 CPU cores.

CAVEATS:

 Warning: if you try to run a 'quick' command (e.g. qplz.pl pwd), yet there is ALSO a file in the directory
  named 'pwd', we will assume 'pwd' is a script and submit that instead of running the standard UNIX 'pwd'.
 This is likely to occur if you're trying to submit a program with arguments, for example "bowtie arg1" and you
  are in the same directory as the actual 'bowtie' exectuable. "qplz.pl" will assume 'bowtie' is a PBS-able script and will
  FAIL, because it then doesn't understand what 'arg1' is supposed to be.
  This can probably be fixed with more smart detection of files or maybe a 'quick job' flag.

 -o and -e (STDERR / STDOUT redirection) is not working yet.
  By default, PBS will write to your job submission directory with a bunch of stuff like "STDIN.o1234".

 Does not yet support other fancy PBS options like joining STDERR and STDOUT.
  --------------

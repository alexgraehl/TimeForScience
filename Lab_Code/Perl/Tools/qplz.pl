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

my $JOBNAME_MAX_LEN = 128; # cut it down to some reasonable size

my @RECOGNIZED_PBS_OPTIONS = ("walltime", "mem", "ncpus"); # The ones we handle in this script

my $PBS_DIRECTIVE_PREFIX = "#PBS"; # <-- the thing at the top of 'your_submitted_script.sh' that means a PBS directive is coming. Usually "#PBS"
my $GLOBAL_WARN_STRING = ""; # record any POSSIBLE problems so that we can print them out at the very end where the user can see them in a convenient summarized format
my $QSUB_EXE = "qsub";
my $QSTAT_EXE = "qstat";

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
		  , "pbs_defaults" => { "ncpus" => {"$UNIX_BIOGRP"   => 1
					   , "$UNIX_GENGRP" => 1 }
					, "mem" => {"$UNIX_BIOGRP"   => 4
						    , "$UNIX_GENGRP" => 4 }
					, "walltime" => {"$UNIX_BIOGRP"   => "23:59:59"
							 , "$UNIX_GENGRP" => "23:59:59" }
				 }
	 );


# another thing to check for: qstat -f :   "job held, too many failed attempts to run"



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

sub printCool($) {			# one required argument
	my ($msg) = @_; chomp($msg);
	print STDERR safeColor("[PROGRESS REPORT]: $msg\n", "white on_blue");
}

sub printProgressWait($) {			# one required argument
	my ($msg) = @_; chomp($msg);
	print STDERR safeColor("[PROGRESS REPORT]: $msg\n", "green on_black");
}

sub printJobTechnicalDetails($) {
	my ($msg) = @_; chomp($msg);
	print STDERR safeColor("[NOTE]: $msg\n", "green on_black");
}

sub printNote($) {			# one required argument
	my ($msg) = @_; chomp($msg);
	print STDERR safeColor("[NOTE]: $msg\n", "white on_blue");
}

sub printBadNews($) {			# one required argument
	my ($msg) = @_; chomp($msg);
	print STDERR safeColor("[ERROR]: $msg\n", "white on_red");
}
sub printBadNewsAndDie($) { printBadNews($_[0]); die $_[0]; }

sub main();

sub quitWithComplaintAboutScript($$$) { my ($scriptFile, $lineNum, $msg) = @_; printBadNews(qq{On line $lineNum of the script file <$scriptFile>, we ran into this problem: $msg}); exit(1); }

sub quitWithComplaintAboutScriptOption($$$) {
	my ($offendingOption, $scriptFile, $lineNum) = @_;
	printBadNews(qq{On line $lineNum of the script file <$scriptFile>:});
	printBadNews(qq{Although you have specified '$offendingOption' on the command line...});
	printBadNews(qq{...that same setting ALSO appears in the script (<$scriptFile>) on line $lineNum.});
	printBadNews(qq{TO FIX THIS: Either remove '$PBS_DIRECTIVE_PREFIX -l $offendingOption'});
	printBadNews(qq{             from line $lineNum in your script...});
	printBadNews(qq{         ... OR omit '--$offendingOption' on the command line...});
	printBadNews(qq{         ... OR run qplz with '--override' to overrule the script.});
	exit(1);
}

sub quitWithUsageError($) { printBadNews($_[0]); printUsage(); printBadNews($_[0]); exit(1); }


sub debugPrint($) {			# one required argument
	my ($msg) = @_; chomp($msg);
	print STDERR safeColor("[DEBUG] $msg\n", "yellow on_red");
}

sub ourWarn($) { my $s = $_[0]; chomp($s); print("[WARNING]: $s\n"); $GLOBAL_WARN_STRING .= "$s\n"; }

sub printUsageAndQuit() { printUsage(); exit(1); }
sub printUsage() { print STDOUT <DATA>; }




sub getAllowedTimeFromQstatOutputText($) {
	# input: the STDOUT results from 'qstat -f THIS_JOB_NAME'
	# output: the hours, minutes, and seconds that this job was allocatd according to 'Resource_List.walltime = SOMETHING'
	my ($qstatOutputText) = @_;
	chomp($qstatOutputText);
	my @qstatLines = split(/\n/, $qstatOutputText);
	my @allowedTimeStrs = grep { /Resource_List.walltime\s*=\s*/i } @qstatLines;
	if (scalar(@allowedTimeStrs) == 1) {
		chomp($allowedTimeStrs[0]);
		if ($allowedTimeStrs[0] =~ m/Resource_List.walltime\s*=\s*(\d+):(\d+):(\d+)/i) {
			return ($1, $2, $3); # <-- results from the match expression above
		} else {
			printBadNews("Weird... we were unable to parse the 'walltime' string from qstat...");
		}
	} else {
		printBadNews("Weird... we were unable to parse the 'qstat' output for some reason...");
	}
	return(undef, undef, undef); # <-- failure
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
	my ($decimalPlaces) = 4; # How many decimal places to print, by default
	$Getopt::Long::passthrough = 1; # ignore arguments we don't recognize in GetOptions, and put them in @ARGV

	my %copt = ();  # Command OPTions. Hash of final values
	foreach my $opt (@RECOGNIZED_PBS_OPTIONS) { $copt{$opt} = undef; }

	my ($pbs_submit_file) = undef;
	my ($shouldOverridePbsDirectives) = 0;
	GetOptions("h|help|?|man"     => sub { printUsageAndQuit(); }
		   , "ncpus|ncpu|c=i" => \$copt{ncpus}
		   , "mem|m=i"        => \$copt{mem}
		   , "walltime|t=s"   => \$copt{walltime}
		   , "f=s"            => \$pbs_submit_file
		   , "override!"      => \$shouldOverridePbsDirectives
		  ) or printUsageAndQuit();

	my $grp = getOurQueueGroup();

	my %optIsAlreadyInFile = (); # remember which $PBS directives were defined in the script
	foreach my $opt (@RECOGNIZED_PBS_OPTIONS) { $optIsAlreadyInFile{opt} = 0; }
	my $numUnprocessedArgs = scalar(@ARGV);

        my $isRunningAScript = (1 == $numUnprocessedArgs and fileIsProbablySomeScript($ARGV[0]));
	if ($isRunningAScript) {
		$pbs_submit_file = $ARGV[0]; # See if the unprocessed argument is a FILENAME (a script to submit)
		# Let's check the script for #PBS directives... if there ARE any, then do not let us override them!
		open my $fff, '<', $pbs_submit_file or die "Cannot read <$pbs_submit_file>...";
		my $lineNum = 0;
		while (my $line = <$fff>) {
			$lineNum++;
			next if (not $line =~ m/^$PBS_DIRECTIVE_PREFIX/i); # skip any non-pbs lines
			(not $line =~ m/\bncpu\b/i) or die "Hey, it looks like you specified 'ncpu' instead of ncpus (plural) in the PBS directive in your script named <$pbs_submit_file>. Fix this!";
			foreach my $opt (@RECOGNIZED_PBS_OPTIONS) {
				if ($line =~ m/^$PBS_DIRECTIVE_PREFIX.*\b$opt\b/i) {
					$optIsAlreadyInFile{$opt} = 1; # remember that this was defined in the script
					if (defined($copt{$opt}) and not $shouldOverridePbsDirectives) { quitWithComplaintAboutScriptOption($opt, $pbs_submit_file, $lineNum); } # Here's a problem, you both defined the option on the command line AND ALSO defined it in the script! We don't know which one to use.
				}
			}
		}
	} else {
		# looks like the user submitted a 'quick command' on the command line, like 'qplz.pl pwd'. So no need to handle overrides.
	}

	foreach my $opt (@RECOGNIZED_PBS_OPTIONS) {
		# Set some defaults... unless they are already defined in the script
		if ( not defined($copt{$opt}) and not $optIsAlreadyInFile{$opt}) {
			$copt{$opt} = $QSETTINGS{pbs_defaults}{$opt}{$grp}; # Set it to the defaults!
		} elsif (defined($copt{$opt}) and     $optIsAlreadyInFile{$opt}) {
			# possible problem --collision! Defined both on command line AND in the file!
			$shouldOverridePbsDirectives or confess "We should have ALREADY exited here with the error checking above.";
			printNote(qq{Since you specified '--override', we are overriding the script's PBS directive for '$opt' and will be using the one you specified on the command line instead (specifically, '${opt}=$copt{$opt}').});
		}
	}

	if (defined($copt{mem})) {
		($copt{mem} =~ m/^\d+(gb|g|)$/i) or quitWithUsageError("(Bad value to --mem / -m): Your memory request (in GIGABYTES) was invalid. You need to speciy an integer number of GB (e.g. '10' or '10gb' or '10g'. You specified this value: $copt{mem}");
		$copt{mem} =~ s/[A-Za-z]//g; # <-- remove any letters from it, now $copt{mem} is purely numeric (i.e., "24gb => 24")
		($copt{mem} =~ m/^\d+$/) or confess("Programming error: somehow failed to make copt{mem} numeric! Offending variable was: $copt{mem}"); # remove any letters from it, now $copt{mem} is purely numeric
		($copt{mem} <= $QSETTINGS{max_mem}{$grp}) or quitWithUsageError("(Bad value to --mem / -m): You requested too much RAM! The maximum you can specify at this specific queue/user/group combination is " . $QSETTINGS{max_mem}{$grp} . " (in gigabytes), but your request was for this number: $copt{mem}");
	}

	if (defined($copt{ncpus})) {
		($copt{ncpus} =~ m/[1-9]\d*/) or quitWithUsageError("(Bad value to --ncpus / -c): You need to specify a non-zero integer number of CPU cores to use. You specified this value: $copt{ncpus}");
		($copt{ncpus} <= $QSETTINGS{max_ncpus}{$grp}) or quitWithUsageError("(Bad value to --ncpus / -c): You specified TOO MANY cpus. The maximum you can specify at this specific queue/user/group combination is " . $QSETTINGS{max_ncpus}{$grp} . ", but your request was for this number: $copt{ncpus}");
	}

	my ($pbs_wall_hr, $pbs_wall_min, $pbs_wall_sec) = (undef, undef, undef); # <-- results from the match expression above
	if (defined($copt{walltime})) {
		if ($copt{walltime} =~ m/^\d+$/) { # if it's literally JUST a single number (the number of hours)
			$copt{walltime} = "$copt{walltime}:00:00"; # I guess it's JUST the number of hours
		}
		($copt{walltime} =~ m/^(\d+):(\d\d):(\d\d)$/) or quitWithUsageError("(Bad value to --walltime / -t): You need to specify a valid walltime in this format: 11:22:33 (hours, minutes, seconds). You specified this value: $copt{walltime}");
		($pbs_wall_hr, $pbs_wall_min, $pbs_wall_sec) = ($1, $2, $3); # <-- results from the match expression above
		($pbs_wall_hr <= $QSETTINGS{max_walltime_hours}{$grp}) or quitWithUsageError("You requested TOO MUCH walltime! The maximum is: $QSETTINGS{max_walltime_hours}{$grp}:00:00 . Try again with a smaller value!");
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

	if (length($jobname) > $JOBNAME_MAX_LEN) { $jobname = substr($jobname, 0, $JOBNAME_MAX_LEN); } # Trim very long job names
	my $qsub_common = qq{$QSUB_EXE }
	  . qq{ -V }
	  . qq{ -N "$jobname" }
	  . qq{ -q "$qdest" }
	  . qq{ -W group_list="$qgrouplist" };
	if (defined($copt{ncpus}))    { $qsub_common .= qq{ -l ncpus=$copt{ncpus} }; }
	if (defined($copt{mem}))      { $qsub_common .= qq{ -l mem="$copt{mem}gb" }; }
	if (defined($copt{walltime})) { $qsub_common .= qq{ -l walltime=$copt{walltime} }; }

	my $pwd = `pwd`; chomp($pwd);

	my $cmd = undef;
	if (defined($pbs_submit_file)) {
		$cmd = qq{$qsub_common $pbs_submit_file};
	} else { # ok, looks like we are just submitting a QUICK job right on the command line
		my $cmdArgs = join(" ", @ARGV); # mash the command line arguments together
		$cmd = qq{echo 'cd "$pwd" && $cmdArgs' | $qsub_common};
	}
	printCool(qq{Submitting this actual command to the queue:\n    ACTUAL JOB SUBMISSION COMMAND -->   $cmd\n});
	my $exitText = `$cmd`; chomp($exitText);
	my $exitCode = $?; # <-- the exit code (i.e., did `$cmd` run properly---same as the result of system($cmd))
	my $jobID = $exitText;
	if (0 == $exitCode) {
		printCool("You can type 'qstat' (basic) or 'qstats' (fancy) to check on the status of your job. --Alex\n");
	} else {
		printBadNewsAndDie("[ERROR] Curses! Something went wrong, and the queue command returned the error code '$exitCode'. Unclear what this means, but something is probably wrong with your input command.\n");
	}

	$jobID =~ m/^(\d+)/ or die "Weird... unexpected exit text ('$exitText'). We thought it would start with a numeric-only job number (like '1234.my-cluster').";
	my $jobNumber = $1; # grab the job number from the $exitText
	my $expectedStdout = "${pwd}/${jobname}.o${jobNumber}"; # expected filename
	my $expectedStderr = "${pwd}/${jobname}.e${jobNumber}"; # expected filename

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

	printJobTechnicalDetails("Your job has been allocated the following:\n");
	defined($copt{ncpus})    and printJobTechnicalDetails("                CPU CORES: $copt{ncpus}\n");
	defined($copt{mem})      and printJobTechnicalDetails("                  MAX RAM: $copt{mem}gb\n");
	defined($copt{walltime}) and printJobTechnicalDetails("                 MAX TIME: $copt{walltime}\n");
	print STDERR "Be aware that your job will be instantly cancelled if it exceeds the MAX RAM or MAX TIME specified above.";
	print STDERR "Now you can type these commands to check your job:\n";
	print STDERR "To check your job 1:  qstats   (print out a color list of jobs that are running\n";
	print STDERR "To check your job 2:  qstat    (print out monochrome list of the above jobs\n";
	print STDERR "To check your job 3:  qstat -u \$USER (print out YOUR jobs only\n";
	print STDERR "If you want to cancel your job (maybe you just realized that it needs more time / RAM):\n";
	print STDERR "    To delete a job: 1. Find the 'Job id' number with 'qstat' (leftmost column)\n";
	print STDERR "    To delete a job: 2. Then use 'qdel ####' (that same number) to cancel it\n";
	print STDERR "    To delete all your jobs (dangerous!): qselect -u \$USER | xargs qdel  <-- deletes all your jobs\n";
	#print STDERR "Ok, now you should run 'qstats' and look for your output in these STDERR / STDOUT files...\n";

	#printJobTechnicalDetails("Your job will be allowed to use $copt{mem} GB of RAM and run for ${pbs_wall_hr} hours and ${pbs_wall_min} minutes before it is cancelled.\n");
	my $secondsToMonitor = 1800;
	my $shouldMonitorForever = 1;
	my $REFRESH_FILE_INTERVAL = 5; # N seconds between filesystem refreshes
	my $REFRESH_QSTAT_INTERVAL = 3; # N seconds between qstat calls

	my $walltimeInSec = undef;

	my $qstatClaimsJobIsDone = 0;

	sleep(1);
	for (my $sec = 1; ($sec <= $secondsToMonitor or $shouldMonitorForever); $sec++) {
		if (0 == $sec % $REFRESH_FILE_INTERVAL) {
			my $refreshCmd = ' FAKEFILE=$(mktemp ./tmp.qplz.XXXXXX.tmp) '. ' && ' . ' /bin/rm $FAKEFILE ';
			system($refreshCmd); # Make and then delete a fake file in order to refresh the filesystem. This lets us catch immediate problems that lead to output files being generated, without having to wait 60 seconds for the filesystem to update.
		}



		foreach my $STDHUH ("STDERR", "STDOUT") {
			my $expectedFilename = ($STDHUH eq "STDERR") ? $expectedStderr : $expectedStdout;
			if (not -e $expectedFilename) {
				#printNote("The $STDHUH file (which is expected to be named <$expectedFilename>) seems not to be available yet. Check for it later.");
			} else {
				my $numLinesToShow = 15;
				my $txt = `tail -n $numLinesToShow "$expectedFilename"`; chomp($txt);
				my $looksLikeError = 0;
				if ($txt =~ m/(not found|not exist|cannot access)/i) {
					$looksLikeError = 1; printBadNews(qq{It looks like either a command or a file was NOT FOUND in your qsub submission! Check the logs for details.\nHOW TO FIX THIS: ***Probably*** you need to specify a full path (like /path/to/my/file.txt) instead of just 'file.txt'. Or, if the problem was an executable that was not found, maybe you did not set your \$PATH variable to include all the special lab-specific tools?});
				}
				if ($txt =~ m/\b(usage:)/i) {
					$looksLikeError = 1; printBadNews(qq{Looks like this program doesn't run properly with the current commands. Better double check the 'usage'.});
				}
				if ($txt =~ m/\b(Segmentation\s+fault|segfault)/i) {
					$looksLikeError = 1; printBadNews(qq{Looks like this program doesn't run properly for some reason!});
				}
				$txt =~ s/^/[Most recent $numLinesToShow lines of $STDHUH] /;
				if ($looksLikeError) {
					printBadNews($txt);
					printBadNews("*"x80);
					printBadNews("      Your job probably FAILED TO RUN!                              ");
					printBadNews("      You should check these files for details:                     ");
					printBadNews("                 STDOUT: $expectedStdout");
					printBadNews("                 STDERR: $expectedStderr");
					printBadNews("*"x80);
					exit(1); # deadly!
				}
				printNote($txt); # Looks like there was no major error, that's good!
			}
		}


		my $shouldUpdateQstat = (!defined($walltimeInSec) or (0 == $sec % $REFRESH_QSTAT_INTERVAL));
		if ($shouldUpdateQstat) {
			my $qstatCmd  = "$QSTAT_EXE -f $jobID";
			debugPrint($qstatCmd);
			my $qstatText = `$qstatCmd 2>&1`; chomp($qstatText);
			my $qstatCode = $?;

			my $JOB_DONE_OK_EXIT_CODE = 35;
			my $JOB_DONE_OTHER_EXIT_CODE = 8960;
			if ($qstatCode =~ m/^($JOB_DONE_OK_EXIT_CODE|$JOB_DONE_OTHER_EXIT_CODE)$/) { # ok exit codes: 35 (job finished) and 8960 (job finished)
				# looks like the job is probably done??
				printCool(qq{Looks like your job has FINISHED! Qstat text is <<$qstatText>>});
				last;
			} elsif (0 == $qstatCode) {
				# probably still running
			} else {
				# got a weird qstat result
				printBadNews(qq{Weird--we didn't recognize the qstat exit code "$qstatCode", which came along with this message: $qstatText"});
				last;
			}

			print $qstatCode . " is the code \n";
			($pbs_wall_hr, $pbs_wall_min, $pbs_wall_sec) = getAllowedTimeFromQstatOutputText($qstatText);
		}
		
		
		my $timeLeftString = "";
		if (defined($pbs_wall_hr)) {
			$walltimeInSec = 3600*$pbs_wall_hr + 60*$pbs_wall_min + $pbs_wall_sec;
			my $remainInSec = $walltimeInSec - $sec;
			my $hhr = floor($remainInSec / 3600);
			my $mmr = floor(($remainInSec / 60) / 60);
			my $ssr = $remainInSec % 60;
			printProgressWait("[$sec] [$hhr:$mmr:$ssr remaining]");
		} else {
			printProgressWait("[$sec] ...");
		}
		sleep(1);
	}

	printNote("Output text may be in these files. Check them as follows:");
	printNote("         STDOUT:   cat $expectedStdout");
	printNote("         STDERR:   cat $expectedStderr");
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

--override : (default: do not override PBS directives in the script file)
  If you have a '#PBS' directive in your script file, like "$PBS_DIRECTIVE_PREFIX -l mem=4gb" and you also specify
  '--mem=8gb' on the command line, you must also pass in '--override' so 'qplz' will know to use your
  command line option instead of the PBS directive in the script.

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

#!/usr/bin/perl -w

#@COMMENT@ mdverify.pl is a script to easily verify a bunch of files with md5 sums. Frequency-of-use rating: 2/10.

# April 2016 by Alex Williams.

use strict; use warnings; use Getopt::Long;
use File::Basename;
use Carp; # for "confess"
use Digest::MD5;
#use Digest::MD5 qw( md5_hex );

#use File::Slurp; # <-- for reading an entire file into memory.
# If the system doesn't have slurp, then:
#  sudo perl -MCPAN -e shell
#  install File::Slurp

use Term::ANSIColor;
$| = 1;  # Flush output to STDOUT immediately.

my $verbose = 0;
my $isDebugging = 0;
my $useColor    = 1;
my $onlyPrintFailures = 0; # by default, print ALL file results
my $shouldRenameBad = 0; # should we rename bad files on the filesystem?
my $shouldDeleteBad = 0; # should we DELETE bad files on the filesystem? May be more dangerous!

sub quitWithUsageError($) { print($_[0] . "\n"); printUsageAndQuit(); print($_[0] . "\n"); }
sub printUsageAndQuit() { printUsageAndContinue(); exit(1); }
sub printUsageAndContinue() {    print STDOUT <DATA>; }
sub debugPrint($) {   ($isDebugging) and print STDERR $_[0]; }
sub verboseWarnPrint($) { if ($verbose) { print STDERR Term::ANSIColor::colored($_[0] . "\n", "yellow on_blue"); } }
sub verboseUpdatePrint($) { if ($verbose) { print STDERR Term::ANSIColor::colored($_[0] . "\n", "black on_green"); } }

my $MD5_OK_COLOR       = "green";
my $MD5_FAILURE_COLOR  = "red";
my $FILE_MISSING_COLOR = "yellow on_black";
my $VERSION = "1.00";

my $G_NOK      = 0;
my $G_NBAD     = 0;
my $G_NMISSING = 0;

sub tryToLoadModule($) {
	my $x = eval("require $_[0]");
	if ((defined($@) && $@)) {
		warn "We FAILED to load module $_[0]. Skipping this module, but continuing with the program.";
		return 0; # FAILURE
	} else {
		$_[0]->import();
		return 1; # SUCCESS
	}
}

my $CAN_USE_COLORS = tryToLoadModule("Term::ANSIColor");
if ($CAN_USE_COLORS) { use Term::ANSIColor; }

sub warnPrint($) { chomp($_[0]); warn(safeColor("[WARNING]: " . $_[0] . "", "yellow on_black")); } # regarding "warn": if the string ends with a newline it WON'T print the line number!

sub safeColor($;$) { # one required and one optional argument
    ## Returns colored text, but only if $CAN_USE_COLORS is set.
    ## Allows you to totally disable colored printing by just changing $CAN_USE_COLORS to 0 at the top of this file
    # Colorstring is OPTIONAL, and can be something like "red on_blue" or "red" or "magenta on_green"
    # Example usage:
    #    *    print STDERR safeColor("This warning message is red on yellow", "red on_yellow");
    my ($message, $color) = @_;
    return (($CAN_USE_COLORS && defined($color)) ? (Term::ANSIColor::colored($message, $color) . Term::ANSIColor::color("reset")) : $message);
}

sub printColorStderr($;$) {
    # prints color to STDERR *UNLESS* it is re-directed to a file, in which case NO COLOR IS PRINTED.
    my ($msg, $col) = @_; # Only prints in color if STDERR is to a terminal, NOT if it is redirected to an output file!
    if (! -t STDERR) { $col = undef; } # no coloration if this isn't to a terminal
    print STDERR safeColor($msg, $col);
}

sub printColorStdout($;$) {
    # prints color to STDOUT *UNLESS* it is re-directed to a file, in which case NO COLOR IS PRINTED.
    my ($msg, $col) = @_; # Only prints in color if STDOUT is to a terminal, NOT if it is redirected to an output file!
    if (! -t STDOUT) { $col = undef; } # no coloration if this isn't to a terminal
    print STDOUT safeColor($msg, $col);
}

sub isMd5($) { # tells you if a string is an md5 sum (32 chars, hex digits only
	return ($_[0]=~m/^[A-Fa-f0-9]{32}$/);
}

sub md5_of_file($) {
	debugPrint("Calculating md5 for the file $_[0]...\n");
	open my $fh, '<', $_[0] or confess "Failure to read file $_[0]...";
	my $md5 = Digest::MD5->new->addfile($fh)->hexdigest;
	close $fh;
	return($md5);
}





$Getopt::Long::passthrough = 1; # ignore arguments we don't recognize in GetOptions, and put them in @ARGV
GetOptions("help|?|man" => sub { printUsageAndQuit(); }
	   , "q|onlybad" => sub { $onlyPrintFailures = 1; }
	   , "nc" => sub { $useColor = 0; } # 'no color'
	   , "c!" => $useColor # default is "yes, use color"
	   , "rename!" => \$shouldRenameBad
	   , "del|delete!" => \$shouldDeleteBad
	   , "vv" => sub { $isDebugging = 1; $verbose = 1; }
	   , "v|version" => sub { print "VERSION $VERSION\n"; exit(0); }
    ) or printUsageAndQuit();

my $numUnprocessedArgs = scalar(@ARGV);
($numUnprocessedArgs >= 1) or quitWithUsageError("[Error] in arguments! You must send at least ONE path (or one filename and '-' for STDIN) to this program.");

my @md5files = @ARGV;

#                   "................................\t................................\t""
my $G_FILE_STATUS = "Status__\tActual_MD5______________________\tExpected_MD5____________________\tFile_path\n";

sub appendStatus($$$) {
	my ($actualMd5, $expectedMd5, $path) = @_;
	my ($statusStr, $shortStatus, $outColor);
	if (!defined($actualMd5)) {
		#             .2345678.2345678.2345678.2345678
		$actualMd5 = "<FILE_MISSING>                  "; # replace 'undefined' with this (32 chars)
		$statusStr    = "[MISSING] ";
		$shortStatus = "missing";
		$G_NMISSING++;
		$outColor = $FILE_MISSING_COLOR;
	} elsif ($actualMd5 ne $expectedMd5) {
		$statusStr = "[MD5_FAIL]";
		$shortStatus = "bad";
		$G_NBAD++;
		$outColor = $MD5_FAILURE_COLOR;
	} else {
		# I guess it was oK!
		$statusStr   = "[OK]      "; # 10 characters
		$shortStatus = "ok";
		$outColor    = $MD5_OK_COLOR;
		$G_NOK++;
	}
	my $line = $statusStr . "\t" . $actualMd5 . "\t" . $expectedMd5 . "\t" . $path . "\n"; # global variable
	$G_FILE_STATUS .= $line;

	if (!$useColor) { $outColor = undef; }

	if ($onlyPrintFailures and $shortStatus eq 'ok') {
		# don't print if the user only wants failures
	} else {
		printColorStdout($line, $outColor);
	}
	#debugPrint($line);
}

foreach my $md5file (@md5files) {
	my $dir = File::Basename::dirname($md5file);
	my $base= File::Basename::basename($md5file);
	(($md5file eq '-') or (-f $md5file and -r $md5file)) or quitWithUsageError("[ERROR]: md5 file '$md5file' ($base) at path ($dir) did not appear to exist!");
	debugPrint("Handling md5 file <$md5file> ($base) at path ($dir)...\n");

	open(my $fh, "<", "$md5file") or confess "[ERROR] Cannot open the (supposedly existing) MD5 checksum file named <$md5file>--check to make sure it really exists!";
	my $lineNum=0;
	while (my $x = <$fh>) {
		$lineNum++;
		chomp($x);
		my $md = undef;
		my $fn = undef;
		my @a = split(/[ ][ ]/, $x);
		if (scalar(@a) != 2) {
			# Ok, this wasn't a normal "two spaces is the delimiter" md5sum output file. Maybe we can still figure it out, though.
			# Sometimes the delimiter is " = ", and in fact, it looks like:   MD5 (filename here) = 12345...
			if ($x =~ m/^MD5 [(]([^)]+)[)]\s+=\s+([0-9a-zA-Z]{32})$/i) {
				@a = ($1, $2); # ok, the regexp matched!
			}
		}
		(scalar(@a) == 2) or confess "[ERROR] Cannot properly parse the (supposedly) MD5 checksum file named <$md5file>---we expect it to have two columns of data, with either TWO SPACES as the delimiter OR possibly ' = ' (equals sign surrounded by a single space--this is apparently what the BSD version of 'md5' prints) as the delimiter!";
		if (isMd5($a[0])) {
			($md, $fn) = ($a[0], $a[1]);
		} elsif (isMd5($a[1])) {
			($md, $fn) = ($a[1], $a[0]); # switch the order
		} else {
			confess "Failure to find an md5 sum on line $lineNum of $md5file...";
		}

		my $filePathOurSystem = undef;
		if (-e $fn) {
			$filePathOurSystem = $fn; # cool, the file exists
		} else {
			my $filebase       = File::Basename::basename($fn); # hmm, maybe only the BASENAME of this file exists, from the current path? (it could be some path like /other_server/whatever/myfile.txt, where only "myfile.txt" is correct
			$filePathOurSystem = $dir . "/" . $filebase;
		}

		my $observedMd5;
		if (not -e $filePathOurSystem) {
			$observedMd5 = undef; # can't find the file on the filesystem...
		} else {
			$observedMd5 = md5_of_file($filePathOurSystem);

			if ($observedMd5 ne $md) {
				# Bad file detected!
				if ($shouldRenameBad) {
					debugPrint("Would be renaming the file now...\n");
				}
				if ($shouldDeleteBad) {
					debugPrint("Would be deleting the file now...\n");
				}
			}

		}
		appendStatus($observedMd5, $md, $filePathOurSystem);
		#debugPrint("Observed md5 = $observedMd5 . Expected md5 = $md . For file $filePathOurSystem.\n");
	}
	close($fh);
	# ff is a full path
}

#print $G_FILE_STATUS;

my $totalFiles = ($G_NOK + $G_NMISSING + $G_NBAD);
printColorStderr(" TOTAL FILES CHECKED:\t$totalFiles\n", undef);
printColorStderr("       [CHECKSUM OK]:\t$G_NOK\n", undef);
printColorStderr("     [CHECKSUM FAIL]:\t$G_NBAD\n", ($G_NBAD > 0 ? $MD5_FAILURE_COLOR : undef));
printColorStderr("     [MISSING FILES]:\t$G_NMISSING\n", ($G_NMISSING > 0 ? $FILE_MISSING_COLOR : undef));

exit(0); # looks like we were successful

################# END MAIN #############################

__DATA__
mdverify.pl [OPTIONS] <list_of_input_md5_files>

by Alex Williams, 2016

This program calculates 'md5' checksums and verifies them against an input file of expected checksums.

 OPTIONS:
-nc:  No color (default: yes, color please)

--onlybad (or -q): Quiet mode. Do not report the files where the MD5 sum matched OK.

--rename (Not implemented yet -- will rename bad files)
--delete (Not implemented yet -- will delete bad files)

--v: verbose
--vv (debugging)

EXAMPLES:
  mdverify.pl --color **/*md5*  (check a bunch of subdirectories)

  mdverify.pl --onlybad md5.txt (do not report the OK files)



KNOWN BUGS:
(None yet)

TO DO:

(Nothing yet)

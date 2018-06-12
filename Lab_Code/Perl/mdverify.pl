#!/usr/bin/env perl

#@COMMENT@ mdverify.pl is a script to easily verify a bunch of files with md5 checksums. It understands the difference between 'md5' and 'md5sum' (Mac vs Linux) and can handle md5 files with either the filename first OR the checksum first, unlike normal md5(sum). Frequency-of-use rating: 2/10.

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


my $EXIT_CODE_OK                       = 0;
my $EXIT_CODE_MISSING_BUT_NO_BAD_FILES = 1;
my $EXIT_CODE_BAD                      = 2;
my $EXIT_CODE_MALFORMED_LINES          = 3; # malformed lines in the checksum file
my $EXIT_CODE_USAGE                    = 31;


my $verbose = 1;
my $isDebugging = 0;
my $useColor    = 1;
my $onlyPrintFailures = 0; # by default, print ALL file results
my $shouldRenameBad = 0; # should we rename bad files on the filesystem?
my $shouldDeleteBad = 0; # should we DELETE bad files on the filesystem? May be more dangerous!

sub quitWithUsageError($)   { print($_[0] . "\n"); printUsageAndContinue(); print($_[0] . "\n"); exit($EXIT_CODE_USAGE); }
sub printUsageAndFail()     { printUsageAndContinue(); exit($EXIT_CODE_USAGE); }
sub printUsageAndContinue() { print STDOUT <DATA>; }
sub debugPrint($)           { ($isDebugging) and print STDERR $_[0]; }
sub verboseWarnPrint($)     { if ($verbose) { print STDERR Term::ANSIColor::colored($_[0] . "\n", "yellow on_red"); } }
sub verboseUpdatePrint($)   { if ($verbose) { print STDERR Term::ANSIColor::colored($_[0] . "\n", "black on_green"); } }

my $RENAME_SUFFIX      = "VERIFICATION_FAILED";
my $MAX_RENAMED_FILES  = 999;
my $MD5_OK_COLOR       = "green";
my $MD5_FAILURE_COLOR  = "red";
my $FILE_MISSING_COLOR = "yellow on_black";
my $VERSION            = "1.00";

my %checkedHash = (); # files that have been checked already

my $G_NOK      = 0;
my $G_NBAD     = 0;
my $G_NMISSING = 0;
my $G_NMALFORMED_LINES = 0;

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
GetOptions("help|?|man" => sub { printUsageAndContinue(); exit($EXIT_CODE_OK); }
	   , "q|onlybad" => sub { $onlyPrintFailures = 1; }
	   , "nc|bw" => sub { $useColor = 0; } # 'no color'
	   , "c!" => $useColor # default is "yes, use color"
	   , "r" => sub { die "-r is not a valid command line option!"; }
	   , "rename!" => \$shouldRenameBad
	   , "delete!" => \$shouldDeleteBad
	   , "vv" => sub { $isDebugging = 1; $verbose = 1; }
	   , "v|version" => sub { print "$VERSION\n"; exit(0); }
    ) or printUsageAndFail();

my $numUnprocessedArgs = scalar(@ARGV);
($numUnprocessedArgs >= 1) or quitWithUsageError("[Error] in arguments! You must send at least ONE md5 filepath to this program. It cannot accept STDIN reading currently.");

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
			# or maybe we can try just ANY spaces... sometimes it's just a tab
			if ($x =~ m/^(\S+)[\s]+(\S{32})$/i) {
				@a = ($1, $2); # ok, the regexp matched (it was just spaces/tabs, no equal sign)
			}
		}
		(scalar(@a) == 2) or confess "[ERROR] Cannot properly parse the (supposedly) MD5 checksum file named <$md5file>---we expect it to have two columns of data, with either TWO SPACES as the delimiter OR possibly ' = ' (equals sign surrounded by a single space) as the delimiter, or failing that, just whitespace as a delimiter. Whitespace in filenames will probably cause this script to FAIL.";
		if    (isMd5($a[0])) {  ($md, $fn) = ($a[0], $a[1]); }
		elsif (isMd5($a[1])) {  ($md, $fn) = ($a[1], $a[0]); } # switch the order
		else {
			verboseWarnPrint("[ERROR on checksum file line $lineNum] Could not find a proper 32-character md5 checksum on line $lineNum of the checksum file '$md5file'! This probably indicates a malformed input file.");
			$G_NMALFORMED_LINES++; # this was a bad line for sure!
			next;
		}
		
		my $filebase       = File::Basename::basename($fn); # hmm, maybe only the BASENAME of this file exists, from the current path? (it could be some path like /other_server/whatever/myfile.txt, where only "myfile.txt" is correct
		my $filePathOurSystem = undef;

		if (-e $fn) {   $filePathOurSystem = $fn; }# cool, the file exists
		else {        $filePathOurSystem = $dir . "/" . $filebase; }
		
		my $observedMd5;
		if (not -e $filePathOurSystem) {
			appendStatus(undef, $md, $filePathOurSystem); # <-- increments 'G_NBAD' and more
		} elsif (exists($checkedHash{$filePathOurSystem}) && defined($checkedHash{$filePathOurSystem})) {
			verboseWarnPrint("[SKIPPING] a duplicate checksum entry for the file '$filePathOurSystem', which we have already checked.");
		} else {
			$checkedHash{$filePathOurSystem} = 1; # remember that we checked this file, in case there are duplicate entries (which occasionally happens!)
			$observedMd5 = md5_of_file($filePathOurSystem);
			if ($observedMd5 ne $md) {
				# Bad file detected!
				if ($shouldRenameBad) {
					debugPrint("Would be renaming the file now...\n");
					my $renamePath = undef;
					for (my $i = 1; $i <= $MAX_RENAMED_FILES; $i++) {
						$renamePath = $filePathOurSystem . "." . sprintf("%04d", $i) . "." . $RENAME_SUFFIX;
						if (not -e $renamePath) { # great, no existing file at this path!
							last;
						}
					}
					if (-e $renamePath) {
						verboseWarnPrint("*** [ERROR] Can't rename a file without overwriting another one (renaming bad files due to '--rename' option) Hit the limit of renaming files! Tried to rename <$filePathOurSystem>, but there were already at least $MAX_RENAMED_FILES files with the rename suffix!");
					} else {
						printColorStderr("* Since --rename was specified, now renaming the MD5-mismatching file <$filebase> to " . File::Basename::basename($renamePath) . "...", $MD5_FAILURE_COLOR);
						system(("mv", "$filePathOurSystem", "$renamePath"));
					}
				}
				if ($shouldDeleteBad) {
					# not currently implemented
					printColorStderr("* Since --delete was specified, now DELETING the MD5-mismatching file <$filePathOurSystem>.", $MD5_FAILURE_COLOR);
					unlink($filePathOurSystem) or verboseWarnPrint("*** [ERROR]: Failed to delete <$filePathOurSystem>");
				}
			}
			appendStatus($observedMd5, $md, $filePathOurSystem); # <-- increments 'G_NBAD' and more
		}
		
		#debugPrint("Observed md5 = $observedMd5 . Expected md5 = $md . For file $filePathOurSystem.\n");
	}
	close($fh);
	# ff is a full path
}

#print $G_FILE_STATUS;

my $totalFiles = ($G_NOK + $G_NMISSING + $G_NBAD + $G_NMALFORMED_LINES);
printColorStderr("   TOTAL FILES CHECKED:\t$totalFiles\n", undef);
printColorStderr("         [CHECKSUM OK]:\t$G_NOK\n", undef);
printColorStderr("       [CHECKSUM FAIL]:\t$G_NBAD\n", ($G_NBAD > 0 ? $MD5_FAILURE_COLOR : undef)); # prints in RED color if there are > 0 bad files!
printColorStderr("       [MISSING FILES]:\t$G_NMISSING\n", ($G_NMISSING > 0 ? $FILE_MISSING_COLOR : undef)); # prints in color if there are > 0 bad files!
printColorStderr(" [MALFORMED CHECKSUMS]:\t$G_NMALFORMED_LINES\n", ($G_NMALFORMED_LINES > 0 ? $MD5_FAILURE_COLOR : undef)); # prints in color if there are > 0 malformed lines

($G_NBAD > 0)             and exit($EXIT_CODE_BAD);  # bad news! there were BAD files. Exit code 2.
($G_NMALFORMED_LINES > 0) and exit($EXIT_CODE_MALFORMED_LINES); # we had some MALFORMED input lines.
($G_NMISSING > 0)         and exit($EXIT_CODE_MISSING_BUT_NO_BAD_FILES); # we had some MISSING files, but no checksum failures. This could mean a lot of things.
exit($EXIT_CODE_OK); # Apparently everything was fine

################# END MAIN #############################

__DATA__
mdverify.pl [OPTIONS] md5file.txt [optionally, additional md5_checksum.txt files]
by Alex Williams, 2016

mdverify.pl calculates 'md5' checksums and verifies them against input file(s) of expected checksums.

You supply one or more "md5.txt" checksum files, and it tries to find the matching files and see if they have the
correct md5sum. See the EXAMPLES section below. There are a few different checksum file formats (e.g., checksum
first / filename first). mdverify tries to find the checksum, and it can handle both tabs and spaces.

OPTIONS:
-q or --onlybad: Quiet mode. Only report BAD and MISSING files, not OK ones.

--rename: Instead of just checking files, if a file has a non-matching MD5 sum, it is renamed
          to have the suffix "VERIFICATION_FAILED". We try to avoid overwriting files in the case of multiple
          failed checksums, so the first one would be e.g. failed_file.0001.VERIFICATION_FAILED. A second
          bad copy of the same file would be named failed_file.0002.VERIFICATION_FAILED , etc, up through a
          maximum of $MAX_RENAMED_FILES.

--delete: DELETES any files with bad checksums. May be dangerous!

--nc|bw:  No color (default: colored text). "bw" is short for black & white.

--vv (debugging)

--version or -v or -V: Print the version number and exit.

EXAMPLES:
  mdverify.pl **/*md5*  (check a bunch of subdirectories)

  mdverify.pl --bw **/*md5*  (check a bunch of subdirectories, no color text)

  mdverify.pl --onlybad md5.txt (do not report the OK files)

  mdverify.pl -q  Download[12]/*/*md5_checksum.txt (similar to the above example, but for multiple files)

  mdverify.pl --rename  */**/md5.txt (rename files where the MD5 sum does not match)

KNOWN BUGS:
  (None yet)

----------


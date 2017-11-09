#!/usr/bin/perl -w

#@COMMENT@ join.pl is a modified version of UNIX join. It can handle un-sorted input and deal with case-insensitive joins. It behaves in a manner that is more similar to what you would expect from a database join. If you want to join multiple files at once, see "join_multi.pl". Frequency-of-use rating: 10/10.

# New version of join.pl by Alex Williams. (This isn't related to the previous UCSC code at all... and probably produces different results! Note that both versions will occasionally produce different results from UNIX join, even on properly sorted input!) UNIX join maybe does the cartesian product sometimes? Anyway, it's probably not what you want.

# Nov 10, 2015: handles Mac '\r'-only input files better. Previously had a bug and output additional "no match" entries no matter what. Whoops, this broken in certain cases. Now it just errors out no matter what if it sees a '\r'
# Now 24, 2015: --multi=intersect added. Now works for multi-file joining (3 or more files).
# Nov 8, 2017: refactored the code to account for some oddities. Removed the special case code for N=2, so now
#              there's only one thing to maintain. Slightly slows down n=2 joins, though. Alas.
use strict;  use warnings;  use diagnostics;
use POSIX      qw(ceil floor);
use List::Util qw(max min);
use Getopt::Long;
use Carp; # backtrace on errors. Has the "confess" function. Use this instead of "die" if you want useful information!
$| = 1; # Always flush text output IMMEDIATELY to the console, don't wait to buffer terminal output! Setting this to zero can cause STDERR and STDOUT to be interleaved in weird ways.

sub tryToLoadModule($) {
    my $x = eval("require $_[0]");
    if ((defined($@) && $@)) { warn "We FAILED to load module $_[0]. Skipping this module, but continuing with the program."; return 0; } # FAILURE
	else { $_[0]->import(); return 1; } # SUCCESS
}

my $SHOULD_USE_COLORS = tryToLoadModule("Term::ANSIColor");
if ($SHOULD_USE_COLORS) { use Term::ANSIColor; }   #print $SHOULD_USE_COLORS . "<-----\n\n\n";

sub main(); # function prototype

my $isDebugging = 0;
my $verbose     = 1; # use 'q' (quiet) to suppress this


sub quitWithUsageError($) { print STDOUT ($_[0] . "\n"); printUsageAndContinue(); verboseWarnPrint($_[0] . "\n"); exit(1); }
sub printUsageAndQuit()   { printUsageAndContinue(); exit(1); }
sub printUsageAndContinue() {    print STDOUT <DATA>; }
sub debugPrint($)         { ($isDebugging) and print STDERR $_[0]; }
sub verboseWarnPrint($)   { if ($verbose) { print STDERR Term::ANSIColor::colored($_[0] . "\n", "yellow on_black"); } } # on_black"); } }
sub verboseUpdatePrint($) { if ($verbose) { print STDERR Term::ANSIColor::colored($_[0] . "\n", "black on_green"); } }

my $KEY_COLUMN_HEADER_STRING_TEXT = "KEY"; # the key is always named KEY no matter what
my $DEFAULT_DELIM = "\t";
my $SPLIT_WITH_TRAILING_DELIMS = -1; # You MUST specify this constant value of -1, or else split will by default NOT split consecutive trailing delimiters! This is super important and belongs in EVERY SINGLE split call.
my $MAX_DUPE_KEYS_TO_REPORT          = 10;
my $MAX_WEIRD_LINE_LENGTHS_TO_REPORT = 10;
my $MAX_BLANK_LINES_TO_REPORT        = 10;
my $MAX_WHITESPACE_KEYS_TO_REPORT    = 50;
my $MAX_MISSING_KEY_COLS_TO_REPORT   = 10;

my $SETDIFF_STR  = "set_subtract";
my $UNION_STR     = "union";
my $INTERSECT_STR = "intersect";

# Global I guess
my $numDupeKeysMultiJoin   = 0;
my $numBlankLines          = 0;
my $numWeirdLengths        = 0;
my $numMissingKeyCols      = 0;
my $numWhitespaceKeys      = 0;
my $numWeirdLineLengths    = 0;

my ($keyCol1,$keyCol2,$keyBoth) = (undef, undef, undef); # indexed from ONE rather than 0!
my ($delim1,$delim2)  = ($DEFAULT_DELIM, $DEFAULT_DELIM); # input deilmiter for files 1 and 2
my ($delimBoth, $outDelim) = (undef, undef); # input delimiter
my $stringWhenNoMatch = undef;
my $allowEmptyKey     = 0; # whether we allow a TOTALLY EMPTY value to be a key (default: no)
my $numHeaderLines    = 0;
my $multiJoinType     = undef;
my $includeFilenameInHeader = 0;
my $shouldSetSubtract = 0;
my $shouldUnixSortKeys    = 0;
my $shouldIgnoreCase  = 0; # by default, case-sensitive
my $sumDuplicates     = 0; # if multiple keys are seen in ONE FILE, try to SUM them.



sub get_key_accounting_for_case_sensitivity($$\%) {
	my ($ignore_case, $key, $mappingHashPtr) = @_;
	return (($ignore_case) ? $$mappingHashPtr{$key} : $key);
}

sub guess_multi_join_type_from_parameters($$$$) {
	my ($join_type, $should_set_subtract, $num_input_files, $na_string) = @_;

	if ($should_set_subtract || (defined($join_type) && $join_type =~ m/^${SETDIFF_STR}/i)) {
		(2 == $num_input_files)  or quitWithUsageError("You can only specify set difference (-v) if there are EXACTLY TWO input files. The first file is the 'gold standard' file, and the keys in the second one are used to subtract out lines from the 'gold' file.");
		!defined($join_type)     or quitWithUsageError("You cannot specify the '-v' (negate / set subtract) operator AND ALSO specify a join type. -v ALREADY implies 'set subtraction' (first file minus keys in second file.)");
		!defined($na_string)     or quitWithUsageError("You cannot specify BOTH set subtraction (the '-v' operator) AND ALSO a missing-value string (e.g. '-o NA' or '-ob'). Try removing the '-o' argument.");
		$join_type = $SETDIFF_STR;
	}
	
	if (!defined($join_type)) {
		($num_input_files <= 2) or quitWithUsageError("Since you are handling 3 or more files, you must specify a join type, either '--multi=union' or '--multi=intersect'. Re-run the command with --multi=i or --multi=u");
		$join_type = defined($na_string) ? $UNION_STR : $INTERSECT_STR;
	}
	if (defined($join_type)) {
		if ($join_type =~ m/^(${UNION_STR}|${INTERSECT_STR}|${SETDIFF_STR})$/i) { } # do nothing---it was already set correctly
		elsif ($join_type =~ m/^u/i)    { $join_type = $UNION_STR; }
		elsif ($join_type =~ m/^i/i) { $join_type = $INTERSECT_STR; }
		else { quitWithUsageError("Unrecognized --multi value: must be either 'union' or 'intersect' (or 'u' or 'i'). Your invalid value was: $join_type"); }
	}
	if (defined($join_type) and ($join_type eq $INTERSECT_STR) and defined($na_string)) {
		quitWithUsageError("You cannot specify a 'string when there is no match' value when you are doing an INTERSECTION. This value would never be used. Remove the '-o' or '--ob' arguments!");
	}
	return $join_type;
}

sub maybeWarn($$$) {
	my ($weirdCount, $maxWeirdCount, $message) = @_;
	if ($weirdCount == $maxWeirdCount) { $message .= " (suppressing further warnings about this issue)"; }
	($weirdCount <= $maxWeirdCount ) and verboseWarnPrint("[WARNING] $message");
}

sub closeSmartFilehandle($) { my($handle)=@_; if ($handle ne *STDIN) { close $handle; } }# Don't close STDIN, but close anything else!
sub openSmartAndGetFilehandle($) {
    # returns a FILEHANDLE. Can be standard in, if the 'filename' was specified as '-'
    # Transparently handles ".gz" and ".bz2" files.
    # This is the MARCH 6, 2013 version of this function.
    my ($filename) = @_;
    if ($filename eq '-') { return(*STDIN); } # <-- RETURN!!!
    my $reader;
    if    ($filename =~ /[.](gz|gzip)$/i)  { $reader = "gzip  -d --stdout $filename |"; }    # Un-gzip a file and send it to STDOUT.
    elsif ($filename =~ /[.]bz2$/i) { $reader = "bzip2 -d -c       $filename |"; }    # Un-bz2 a file and send it to STDOUT
    elsif ($filename =~ /[.]xz$/i)  { $reader = "xz -d -c       $filename |"; }    # Un-xz a file and send it to STDOUT
    elsif ($filename =~ /[.]zip$/i) { $reader = "unzip -p $filename |"; } # Un-regular-zip a file and send it to STDOUT with "-p": which is DIFFERENT from -c (-c is NOT what you want here). See 'man unzip'
    else                            { $reader = "$filename"; }  # Default: just read a file normally
    my $fh;
    open($fh, "$reader") or die("Couldn't read from <$filename>: $!");
    return $fh;
}

sub assert_line_has_unix_line_endings($$$) {
	my ($filename, $lineToCheck, $lineNum) = @_;
	($lineToCheck !~ m/\r/) or die "ERROR: Exiting! The file <$filename> appears to have either WINDOWS-STYLE line endings or MAC-STYLE line endings ( with an '\\r' character) (as seen on line $lineNum).\nThis behavior is often seen when files are saved by Excel. You will need to manually convert the line endings from Mac / Win format to UNIX. Search online for a way to do this. We cannot handle this character automatically at this point in time, and are QUITTING.\n";
	return 1;
}

sub is_line_array_too_weird_to_use(\@$$$) {
	my ($lineArrPtr, $lnum, $filename, $keyColIndexedFromOne) = @_;
	if (0 == scalar(@$lineArrPtr)) {
		maybeWarn($numBlankLines++, $MAX_BLANK_LINES_TO_REPORT, "on line $lnum in file <$filename> was blank. Skipping it.");
		return 1; # Disqualify this line, since it is TOTALLY EMPTY
	} elsif (scalar(@$lineArrPtr) < $keyColIndexedFromOne) { # The line was SHORTER than the demanded key column! $keycol is indexed from 1, not 0;
		maybeWarn($numMissingKeyCols++, $MAX_MISSING_KEY_COLS_TO_REPORT, "line $lnum was missing the key column in file <$filename>. (Key column: $keyColIndexedFromOne, Columns on line: " . scalar(@$lineArrPtr) . ").");
		# Disqualify this line, since it has no key!
		#$key = "";  # totally blank line maybe? Or at least, no key.
		#@valArr = ();
		return 1; # too weird
	} else {
		return 0; # not too weird
	}
}

sub handleMultiJoin($$$$$$) {
	my ($k1, $k2, $filenameArrPtr, $mergeType, $d1, $d2) = @_;
	
	my %datHash               = (); # this stores all the lines, and gets very large. Key = filename, value = hash with a second key = line key, value = array of data on that line
	my %longestLineInFileHash = (); # key = filename, value = how long the longest line in that file is
	my %keysHash              = (); # key = the keys seen in ALL files, value = order that this key was added (0 is the first key, 1 is the second, etc...)
	my %ocKeyHash             = (); # Maps keys back to their ORIGINAL CASE versions, in case we were doing case-insensitive (case insensitive) key output.
	my @keysInOrder           = ();
	
	($mergeType =~ m/^(${UNION_STR}|${INTERSECT_STR}|${SETDIFF_STR})$/i) or die "Programming error: multi-intersection type must be '$INTERSECT_STR' or '$UNION_STR'! But it was this: $mergeType";
	my $filenameHeaderDelim = "::"; # example:  "Filename1::headerCol1   Filename2::headerCol2"
	my $na = (defined($stringWhenNoMatch)) ? $stringWhenNoMatch : ""; # use a blank value if (global) $stringWhenNoMatch is not defined
	my %headHash = ();
	my $is_intersection     = ($mergeType eq ${INTERSECT_STR});
	my $is_union            = ($mergeType eq ${UNION_STR});
	my $is_setdiff         = ($mergeType eq ${SETDIFF_STR});
	my %numFilesWithKeyHash = (); # Counts the number of files that we saw this line in. Only used if this is an intersection and not a union
	my $numFilesOpened      = 0;
	my $numFilesExpected    = scalar(@$filenameArrPtr);
	my $totalLinesReadAcrossAllFiles = 0;
	foreach my $filename (@$filenameArrPtr) {
		my $numItemsOnFirstLine = undef;
		$numFilesOpened++; # ok, remember that we read a file
		%{$datHash{$filename}} = (); # new hash value is an ARRAY for this
		if (0==${numHeaderLines} && $includeFilenameInHeader) {
			@{$headHash{$filename}} = ( () ); # zero header lines, but in this case we'll initialize the header anyway
		} else {
			@{$headHash{$filename}} = (()x$numHeaderLines); # it's one array element per line
		}
		$longestLineInFileHash{$filename} = 0; # longest line is length 0 to start...
		#print STDERR "Handling file named $filename ...\n";
		($filename ne '-') or die "If you are multi-joining, you CANNOT read input from STDIN. Sorry.";
		(-e $filename) or die "Cannot read input file $filename.";
		my $fh = openSmartAndGetFilehandle($filename);
		my $lnum = 0;

		my $delim  = ($numFilesExpected == 2 && $numFilesOpened == 2 && defined($d2)) ? $d2 : $d1; # If there are EXACTLY TWO input files, then we allow delim1 and delim2 to be separately specified
		my $keycol = ($numFilesExpected == 2 && $numFilesOpened == 2 && defined($k2)) ? $k2 : $k1; # If there are EXACTLY TWO input files, then we allow key1 and key2 to be separately specified
		(0!=${keycol} && defined($keycol)) or die "Keycol can't be 0 or undefined -- it's indexed with ONE as the first element.";
		(defined($delim))                  or die "delim can't be undefined.";
		foreach my $line (<$fh>) {
			$lnum++;
			$totalLinesReadAcrossAllFiles++;
			($totalLinesReadAcrossAllFiles % 25000 == 0) and verboseUpdatePrint("[STATUS UPDATE] Read $totalLinesReadAcrossAllFiles lines across $numFilesOpened files...");
			if ($lnum == 0) { assert_line_has_unix_line_endings($filename, $line, $lnum); }
			chomp($line);
			my @vals = split($delim, $line, $SPLIT_WITH_TRAILING_DELIMS); # split up the line
			if (is_line_array_too_weird_to_use(@vals, $lnum, $filename, $keycol)) {
				next; # line is too weird for us, maybe it's missing a key or something.
			}
			if (!defined($numItemsOnFirstLine)) { $numItemsOnFirstLine = scalar(@vals); }
			($numItemsOnFirstLine == scalar(@vals)) or maybeWarn($numWeirdLineLengths++, $MAX_WEIRD_LINE_LENGTHS_TO_REPORT, "The number of elements on each line of this file was NOT CONSTANT. The first line had $numItemsOnFirstLine columns, but line number $lnum in file <$filename> had " . scalar(@vals) . " instead. Continuing anyway.");

			my $key = $vals[($keycol-1)];  # get the correct key, since it might not be the first column, I guess!

			if ($shouldIgnoreCase) {
				$ocKeyHash{uc($key)} = $key; # save the original key
				$key = uc($key); # upper-case this key for any future work
			}
			
			splice(@vals, ($keycol-1), 1); # <-- Delete the key column from the array! Splice MODIFIES the array---it REMOVES the key column, keep everything else!
			# ==========================================
			# Check for a certain (non-fatal) warning issues.
			($key !~ /\s/) or maybeWarn($numWhitespaceKeys++, $MAX_WHITESPACE_KEYS_TO_REPORT, "on line $lnum in file <$filename>, the key <$key> had *whitespace* in it. This is often unintentional, but is not necessarily a problem!");
			# ==========================================
			
			# Note: line length does NOT include the key as an element. So it can be zero!
			if ($longestLineInFileHash{$filename} < scalar(@vals)) { $longestLineInFileHash{$filename} = scalar(@vals); }
			my $this_is_a_header_line = ($lnum <= $numHeaderLines);
			if ($this_is_a_header_line) {
				@{${$headHash{$filename}}[$lnum-1]} = ($includeFilenameInHeader) ? map{"${filename}${filenameHeaderDelim}$_"}@vals : @vals; # It is ok if @valsy has zero elements.
			} else {
				if (($lnum == 1) && ($numHeaderLines == 0) && ($includeFilenameInHeader)) {
					# Just for the very first line, if we DO NOT have a header line, we'll make one anyway if the user specified '--fnh'
					# Include the filename in the 'header', even though there isn't a header line per se---we create a new one.
					@{${$headHash{$filename}}[$lnum-1]} = map{"${filename}"}@vals; # <-- initialize a header line to something like  KEY   file1.txt   file1.txt   file2.txt  file3.txt etc...
				}
				
				if ((!$allowEmptyKey) && ("" eq $key)) {
					verboseWarnPrint("[WARNING] Skipping a blank key on line $lnum of file <$filename>!");
				} elsif (exists($datHash{$filename}{$key})) {
					# Already saw this key!
					maybeWarn($numDupeKeysMultiJoin++, $MAX_DUPE_KEYS_TO_REPORT, "on line $lnum, we saw a duplicate of key <$key> (in file <$filename>). We are only keeping the FIRST instance of this key.");
				} else {
					#print "File [$filename] added key <<$key>> ...\n";
					@{$datHash{$filename}{$key}} = @vals; # save this key with the value being the rest of the array
					if (!exists($keysHash{$key})) {
						$keysHash{$key} = 1; #${numKeysAdded}++;
						push(@keysInOrder, $key);
					} # Apparently this is a TOTALLY new key never before seen in any file, so save the order that this key was seen in!
				}
			}
		}
		closeSmartFilehandle($fh);
	}

	my $nHeaderLinesAccountingForFNH = ($includeFilenameInHeader and $numHeaderLines == 0) ? 1 : $numHeaderLines;  # basically, if include-filename-in-header is here, BUT we don't have any other header lines, then this value should be at least one!
	my $firstFile = $$filenameArrPtr[0];

	my @keysToCheck;
	if ($is_setdiff) {
		@keysToCheck = ();
		foreach my $candidate_key (@keysInOrder) { # if it's a set difference, we only try to print the keys in the FIRST file
			my $key_is_in_first_file = exists($datHash{$firstFile}{$candidate_key});
			if ($key_is_in_first_file) { push(@keysToCheck, $candidate_key); }
		}
	} else {
		@keysToCheck = @keysInOrder;
	}

	if ($shouldUnixSortKeys) { @keysToCheck = sort(@keysToCheck); }

	my @filesToCheckForPrinting = ($is_setdiff) ? ($firstFile) : @$filenameArrPtr;
	for (my $headerIndex = 0; $headerIndex < $nHeaderLinesAccountingForFNH; $headerIndex++) {
		my @head = (${KEY_COLUMN_HEADER_STRING_TEXT}); # array with one element to start
		for my $filename (@filesToCheckForPrinting) {
			my $thisHeadArrPtr = \@{${$headHash{$filename}}[$headerIndex]};
			my $numElemsToPad = $longestLineInFileHash{$filename} - scalar(@$thisHeadArrPtr);
			push(@head, @$thisHeadArrPtr, ($na)x$numElemsToPad);
		}
		print(join($outDelim, @head)."\n");
	}

	foreach my $k (@keysToCheck) {
		# ========== See if the key is in EVERY SINGLE file, but only if we are computing an intersection! ======
		my $key_is_in_every_file = 0;
		if ($is_intersection || $is_setdiff) { # gotta compute this for intersection AND ALSO set diff (counterintuitively, since for set diff, we DISQUALIFY any intersecting keys)
			my $numFilesWithKey = 0;
			foreach my $filename (@$filenameArrPtr) {
				if (exists($datHash{$filename}{$k})) { $numFilesWithKey++; } # see if this filename/key combination exists!
			}
			$key_is_in_every_file = ($numFilesWithKey == $numFilesOpened);
			#print "Skipping key <$k>: it was only in $numFilesWithKey files of $numFilesOpened, so we are skipping it for the intersection.\n";
		}
		#print "U:    $is_union        I: $is_intersection      D: $is_setdiff\n";
		my $key_has_succeeded = (($is_union) || ($is_intersection && $key_is_in_every_file) || ($is_setdiff && !$key_is_in_every_file));
		#print "$k ------> every file? $key_is_in_every_file\n";
		#print "$k ------> enough files? $key_has_succeeded\n";
		if ($key_has_succeeded) {
			my $printableKey = get_key_accounting_for_case_sensitivity($shouldIgnoreCase, $k, %ocKeyHash);
			my @outLine     = ($printableKey); # output line. starts with just the key and nothing else
			foreach my $filename (@filesToCheckForPrinting) {
				if (exists($datHash{$filename}{$k})) {
					my $numElemsToPad = $longestLineInFileHash{$filename} - scalar(@{$datHash{$filename}{$k}});
					push(@outLine, @{$datHash{$filename}{$k}}, ($na)x$numElemsToPad);
				} else {
					# This specific file didn't have this key!
					push(@outLine, ($na)x$longestLineInFileHash{$filename});
				}
			}
			print join($outDelim, @outLine), "\n"; # <-- somehow this results in "uninitialized value" warning sometimes...
		}
	}
} # end of handleMultiJoin(...)
# ========================== MAIN PROGRAM HERE

sub main() { # Main program
	$Getopt::Long::passthrough = 1; # ignore arguments we don't recognize in GetOptions, and put them in @ARGV
	GetOptions("help|?|man" => sub { printUsageAndQuit(); }
			   , "k=i" => \$keyBoth
			   , "f=i" => sub { quitWithUsageError("-f is not an option to this script! If you want to specify field separators, use -d (delim). If you want to specify keys, use -1 and -2 to pick the key columns.\n"); }
			   , "q" => sub { $verbose = 0; } # q = quiet
			   , "1=s" => \$keyCol1
			   , "2=s" => \$keyCol2
			   , "d1=s" => \$delim1
			   , "d2=s" => \$delim2
			   , "t|d|delim=s" => \$delimBoth # -t is the regular UNIX join name for this option
			   , "o=s"  => \$stringWhenNoMatch # -o "0.00" would be "print 0.00 for missing values"
			   , "ob!"  => sub { $stringWhenNoMatch = ''; } # shortcut for just a blank when there's no match. Default is to OMIT lines with no match.
			   , "do=s" => \$outDelim
			   , "v|neg!" => \$shouldSetSubtract
			   , "sort"  => sub { quitWithUsageError("--sort is not an argument, you want to instead type '--unix-sort' or '--usort' to emphasize that this is NOT REGULAR ALPHABETICAL sorting order!"); }
			   , "unix-sort|usort!"  => \$shouldUnixSortKeys
			   #	   , "sum|sum-duplicates!" => \$sumDuplicates # for numeric values, we can also SUM the duplicate values on a per-file basis instead of just overwriting them
			   , "h|header|headers=i" => \$numHeaderLines
			   , "fnh|filename-in-header" => \$includeFilenameInHeader # only for multi-join
			   , "eok|allow-empty-key!" => \$allowEmptyKey # basically, do we skip blank lines?
			   , "i|ignore-case!" => \$shouldIgnoreCase # Ignore case in the keys
			   , "multi=s" => \$multiJoinType
			   , "union!"     => sub { $multiJoinType = $UNION_STR; }
			   , "intersect!" => sub { $multiJoinType = $INTERSECT_STR; }
			   , "debug!" => \$isDebugging
			  ) or printUsageAndQuit();

	my $numUnprocessedArgs = scalar(@ARGV);
	($numUnprocessedArgs >= 2) or quitWithUsageError("[Error] in arguments! You must send at least TWO filenames (or one filename and '-' for STDIN) to this program. Example: join.pl FILE1.txt FILE2.txt > OUTPUT.txt");

	my @files = @ARGV;			#my ($file1,$file2) = ($files[0], $files[1]);
	my $numFilesToJoin = scalar(@files);
	$multiJoinType = guess_multi_join_type_from_parameters($multiJoinType, $shouldSetSubtract, $numFilesToJoin, $stringWhenNoMatch);

	if (defined($keyBoth)) {
		(!defined($keyCol1) && !defined($keyCol2)) or quitWithUsageError("You cannot specify both -k (key) AND ALSO -1 (key for file 1) or -2 (key for file 2). -k sets both -1 and -2. Pick one or the other!");
		$keyCol1 = $keyBoth; $keyCol2 = $keyBoth;
	} else {
		$keyCol1 = (defined($keyCol1)) ? $keyCol1 : 1; # default value is 1
		$keyCol2 = (defined($keyCol2)) ? $keyCol2 : 1; # default value is 1
	}
	($keyCol1 =~ m/[0-9]+/ && $keyCol1 >= 0) or quitWithUsageError("[ERROR]: Key1 (-1 argument) to join.pl (which was specified as '$keyCol1') CANNOT BE ZERO or less than zero! These indices are numbered from ONE and not zero!");
	($keyCol2 =~ m/[0-9]+/ && $keyCol2 >= 0) or quitWithUsageError("[ERROR]: Key2 (-2 argument) to join.pl (which you specified as '$keyCol2') CANNOT BE ZERO or less than zero! These indices are numbered from ONE and not zero!");

	## ================ SET SOME DEFAULT VALUES ============================
	if (defined($delimBoth)) { ($delim1, $delim2) = ($delimBoth, $delimBoth); } # If "delimBoth" was specified, then set both of the input delimiters accordingly.

	if (!defined($outDelim)) { # Figure out what the output delimiter should be, if it wasn't explicitly specified.
		if ($delim1 eq $delim2) { $outDelim = $delim1; } # default: set the output delim to whatever the input delim was
		else { $outDelim = $DEFAULT_DELIM; } # otherwise, set it to the default delimiter
	}
	## ================ DONE SETTING SOME DEFAULT VALUES ====================

	## ================ SANITY-CHECK A BUNCH OF VARIABLES ==================
	foreach my $ff (@files) {
		((-f $ff) && (-r $ff)) or die "[ERROR]: join.pl cannot join these two files, because the file '$ff' did not exist or could not be read!"; # Specified files must be either - (for stdin) or a real, valid, existing filename)
	}

	if (defined($stringWhenNoMatch)) { # replace any "\t" with actual tabs! Lets the user specify a delimiter as -d '\t'
		$stringWhenNoMatch =~ s/[\\][t]/\t/g; # replace a backslash-then-t with a tab
		$stringWhenNoMatch =~ s/[\\][n]/\n/g; # replace a backslash-then-n with a newline
		$stringWhenNoMatch =~ s/[\\][r]/\r/g; # replace a backslash-then-r with a CR return. Why would you do this, I wonder!
	}
	
	handleMultiJoin($keyCol1, $keyCol2, \@files, $multiJoinType, $delim1, $delim2);
	($numDupeKeysMultiJoin > 0) and verboseWarnPrint("$numDupeKeysMultiJoin duplicate keys were skipped in the multi-joining.");
}

main();
END {
	($SHOULD_USE_COLORS) and print STDERR Term::ANSIColor::color("reset"); # Restores terminal text to its normal color.
	# DO NOT PRINT TO STDOUT --> print STDOUT Term::ANSIColor::color("reset"); <-- this prints weird characters to your results!
}
exit(0);

################# END MAIN #############################

__DATA__
syntax: join.pl [OPTIONS] LOOKUP_FILE  HUGE_DICTIONARY_FILE

This script takes two (or more) plain text files as input and
produces a new table that is a join of FILE1 and FILE2 (and perhaps additional files).

In other words, it does a database join on the inputs.

By default, files are assumed to be tab-delimited with the keys in the first column.

Unlike the UNIX "join", join.pl does NOT require sorted keys.

CAVEATS:

* If the DICTIONARY_FILE contains several lines with the same key,
  only the *LAST* key read will actually ever be used. Thus, join.pl
  does NOT exactly duplicate the function of
  GNU "join"--in particular it does not output the cross-product
  in multiple-key situations.

* join.pl reads the ENTIRE contents of each file into memory!
  It may be unsuitable for joining very large (> 1000 MB) files.

OPTIONS:
-1 COLNUM: Find the key in this column from FILE1. Default 1. Indexed from 1.
         Join.pl only supports ONE key field, unlike unix join.
-2 COLNUM: Same as above, but for finding the key in file 2. Default: 1.
-o FILLER: Print keys even if they do not have a match in FILE2: use FILLER
          as the filler text to print for missing values not present in FILE2.
          Example:  join -o "0.0" FILE1.txt FILE2.txt > OUTPUT_with_zeros.txt
--ob: Set a blank value for the FILLER above. Same as -o ''.
-h INT: Number of header lines to print verbatim (without joining) from FILE1. (default: 0)
          If you want filenames in the header, also specify the "--fnh" option.
--fnh or --filename-in-header: Print the filename in the header. Very useful!
           Example:    KEY    file1::col1.txt   file2::col1.txt     file2::col2.txt
--multi=union or --union:  Run a multi-file join, report the union of keys seen in ANY file.
--multi=intersect or --intersect: Run a multi-file join, report only keys that match in ALL files.
            For set difference, see the '-v' option (which only works with exactly two files)
--unix-sort or --usort: Print the keys in UNIX-sort order.
    Default behavior is to print the keys in the order they were first encountered.
-v or --neg: Set subtraction: print keys from FILE1 that are not in FILE2.
            Cannot specify both this AND ALSO '-ob' or '-o'.
-d DELIM or --delim=DELIM or -t DELIM:
            Set the input delimiters for both FILE1 and FILE2 to DELIM (default: tab)
            Equivalent to setting both --d1 and --d2.
--d1=DELIM: Set the input delimiter for FILE1 to DELIM (default: tab).
--d2=DELIM: Set the input delimiter for FILE2 to DELIM (default: tab).
--do=DELIM: Set the OUTPUT delimiter to DELIM. (default: same as input delim)
            You can set the output delimiter explicitly to a tab with --do='\t'.
--allow-empty-key or --eok: (Default: do not allow it). Allows a blank ("") key to be joined.
            Note: even when the empty key is NOT allowed for matching, if we do an outer join or
            negation, we will still print items from FILE1 where the key was blank.
-i or --ignore-case: (Case-insensitive join)
-q or --quiet : No verbose output. May hide some useful warning messages!

EXAMPLE USAGE SCENARIOS:

join.pl -1 1 -2 2 file1--key_in_first_col.txt  file2--key_in_second_col.compressed.gz > join.output.txt
  Print lines from file1 that are ALSO in file2, and append the data from file2 in the output.
  This is the most standard-plain-vanilla join, and should be similar to the unix "join" results.

join.pl -o "NO_MATCH_HERE_I_SEE" -1 4 -2 1 file_with_key_in_fourth_col.compressed.bz2  file2.gz > join.with.unmatched.rows.txt
  Also print the un-matching lines from the first file (lines with no match will say "NO_MATCH_HERE_I_SEE")

join.pl --multi=union   file1.txt   file2.txt   file3.txt
  Join MORE than two files. Has different behaviors from joining exactly two files. Some options do not work
  with multi-file joining (example: --neg does not work).
  Be aware of --fnh (filename in header), which can help a lot when multi-joining files with otherwise-identical keys!

cat myfile.txt | join.pl - b.txt | less -S
  Read from STDIN (use a '-' instead of a filename!), and pipe into the program "less".

Detailed example:

If you have the following two files:

File1 (tab-delimited)       <-- FIRST FILE
    AAA   avar    aard
    ZZZ   zebra                     <-- Note the differing number of items per line!
    BBB   beta    bead    been          This *usually* indicates a problem in your input data.
    MMM   man     most                  Alex wrote "table-no-ragged.py" to pad any blank spaces,
                                        so you can use that if your file has "ragged" edges.

File2 (tab-delimited)       <-- SECOND FILE
    AAA   111
    BBB   222
    CCC   333

Then the result of a regular join.pl invocation...
   join.pl File1 File2
...is:
    AAA   avar    aard    111             <-- Note: 4 columns in all
    BBB   beta    bead    been    222     <-- Note: 5 columns in all

The result of an outer join (i.e., use File1 as the "master" file)...
   join.pl -o "NONE!" File1 File2
Results in:
    AAA    avar     aard      111
    ZZZ    zebra    NONE!
    BBB    beta     bead      been      222
    MMM    man      mouse     NONE!

Note that if you switch the order of File1 and File2 for an outer join...
...you get a different result (the first file specifies the keys)!

'--fnh' is "filename in header," which is a convenient way of labeling all of the
columns in the (common) case where many files have identical column headers.

join.pl --fnh -o "NOPE!" File2.txt File1.txt
  Results in:
    KEY    File1.txt File2.txt
    AAA    111       avar      aard
    BBB    222       beta      bead       been
    CCC    333       NOPE!

----------
Note that UNIX join will behave differently from
join.pl on the input FILE1 and FILE2 given below:

FILE1: (tab-delimited)
Alpha   1   2
Alpha   3   4

FILE2: (tab-delimited)
Alpha   a   first
Alpha   b   middle
Alpha   c   last

----------


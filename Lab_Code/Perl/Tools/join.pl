#!/usr/bin/perl -w

#@COMMENT@ join.pl is a modified version of UNIX join. It can handle un-sorted input and deal with case-insensitive joins. It behaves in a manner that is more similar to what you would expect from a database join. If you want to join multiple files at once, see "join_multi.pl". Frequency-of-use rating: 10/10.

# New version of join.pl by Alex Williams. (This isn't related to the previous UCSC code at all... and probably produces different results! Note that both versions will occasionally produce different results from UNIX join, even on properly sorted input!) UNIX join maybe does the cartesian product sometimes? Anyway, it's probably not what you want.

use strict; use warnings; use Getopt::Long;
use Term::ANSIColor;
$| = 1;  # Flush output to STDOUT immediately.

my $isDebugging = 0;

my $verbose = 1; # use 'q' (quiet) to suppress this

sub quitWithUsageError($) { print($_[0] . "\n"); printUsageAndQuit(); print($_[0] . "\n"); }
sub printUsageAndQuit() { printUsageAndContinue(); exit(1); }
sub printUsageAndContinue() {    print STDOUT <DATA>; }
sub debugPrint($) {   ($isDebugging) and print STDERR $_[0]; }
sub verboseWarnPrint($) {
    if ($verbose) { print STDERR Term::ANSIColor::colored($_[0] . "\n", "yellow on_blue"); }
}
my $keyCol1 = 1; # indexed from ONE rather than 0!
my $keyCol2 = 1; # indexed from ONE rather than 0!

my $DEFAULT_DELIM_IF_NOTHING_ELSE_IS_SPECIFIED = "\t";
my $SPLIT_WITH_TRAILING_DELIMS = -1; # You MUST specify this constant value of -1, or else split will by default NOT split consecutive trailing delimiters! This is super important and belongs in EVERY SINGLE split call.

my $MAX_DUPE_KEYS_TO_REPORT = 10;
my $MAX_WEIRD_LINE_LENGTHS_TO_REPORT = 10;

my $delim1 = $DEFAULT_DELIM_IF_NOTHING_ELSE_IS_SPECIFIED; # input deilmiter for file 1
my $delim2 = $DEFAULT_DELIM_IF_NOTHING_ELSE_IS_SPECIFIED; # input delimiter for file 2
my $delimBoth   = undef; # input delimiter
my $outputDelim = undef;

my $shouldNegate = 0; # whether we should NEGATE the output
my $shouldIgnoreCase = 0; # by default, case-sensitive

my $stringWhenNoMatch = undef;

my $allowEmptyKey = 0; # whether we allow a TOTALLY EMPTY value to be a key (default: no)

my $shouldKeepKeyInOriginalPosition = 0; # whether we should KEEP the key in whatever column it was found in, instead of moving it to the front of the line.

my $numHeaderLines = 0;

$Getopt::Long::passthrough = 1; # ignore arguments we don't recognize in GetOptions, and put them in @ARGV
GetOptions("help|?|man" => sub { printUsageAndQuit(); }
	   , "q" => sub { $verbose = 0; }
	   , "1=s" => \$keyCol1
	   , "2=s" => \$keyCol2
	   , "d1=s" => \$delim1
	   , "d2=s" => \$delim2
	   , "t|d|delim=s" => \$delimBoth # -t is the regular UNIX join name for this option
	   , "o=s"  => \$stringWhenNoMatch
	   , "ob!"  => sub { $stringWhenNoMatch = ''; } # shortcut for just a blank when there's no match. Default is to OMIT lines with no match.
	   , "do=s" => \$outputDelim
	   , "neg!" => \$shouldNegate
	   , "h|header|headers=i" => \$numHeaderLines
	   , "nrk|no-reorder-key!"  => \$shouldKeepKeyInOriginalPosition
	   , "eok|allow-empty-key!" => \$allowEmptyKey
	   , "i|ignore-case!" => \$shouldIgnoreCase
	   , "debug!" => \$isDebugging
    ) or printUsageAndQuit();

my $numUnprocessedArgs = scalar(@ARGV);
($numUnprocessedArgs == 2) or quitWithUsageError("[Error] in arguments! You must send exactly TWO filenames (or one filename and '-' for STDIN) to this program. Example: join.pl FILE1.txt FILE2.txt > OUTPUT.txt");

my $file1 = $ARGV[0]; # first unprocessed argument
my $file2 = $ARGV[1]; # second unprocessed argument

## ================ SET SOME DEFAULT VALUES ============================
if (defined($delimBoth)) { # If "delimBoth" was specified, then set both of the input delimiters accordingly.
    $delim1 = $delimBoth; $delim2 = $delimBoth;
}

if (!defined($outputDelim)) { # Figure out what the output delimiter should be, if it wasn't explicitly specified.
    if (defined($delimBoth)) { $outputDelim = $delimBoth; } # default: set the output delim to whatever the input delim was
    elsif ($delim1 eq $delim2) { $outputDelim = $delim1; } # or we can set it to the manually-specified delimiters, if they are the SAME only
    else { $outputDelim = $DEFAULT_DELIM_IF_NOTHING_ELSE_IS_SPECIFIED; } # otherwise, set it to the default delimiter
}
## ================ DONE SETTING SOME DEFAULT VALUES ====================

## ================ SANITY-CHECK A BUNCH OF VARIABLES ==================
foreach my $ff ($file1, $file2) {
    (($ff eq '-') or (-f $ff and -r $ff)) or die "[ERROR]: join.pl cannot join these two files, because the file '$ff' did not exist or could not be read!"; # Specified files must be either - (for stdin) or a real, valid, existing filename)
}

(!$shouldNegate or !defined($stringWhenNoMatch)) or quitWithUsageError("[Error] in arguments! Cannot specify both --neg AND -o or --ob, because it doesn't make sense to both '--neg (negate)' the join AND ALSO specify '-o' or '--ob' -- the outer join specifies that we should print lines REGARDLESS of match, whereas the --neg specifies that we should ONLY print lines with no match. You cannot specifiy both of these options at the same time.");
($keyCol1 != 0) or quitWithUsageError("[ERROR]: Key1 (-1 argument) to join.pl CANNOT BE ZERO! These indices are numbered from ONE and not zero!");
($keyCol2 != 0) or quitWithUsageError("[ERROR]: Key2 (-2 argument) to join.pl CANNOT BE ZERO! These indices are numbered from ONE and not zero!");

if (defined($stringWhenNoMatch)) { # replace any "\t" with actual tabs! No idea why it doesn't work on the command line otherwise
    $stringWhenNoMatch =~ s/[\\][t]/\t/g; # replace a SINGLE backslash-then-t with a tab
    $stringWhenNoMatch =~ s/[\\][n]/\n/g; # replace a SINGLE backslash-then-n with a newline
    $stringWhenNoMatch =~ s/[\\][r]/\n/g; # replace a SINGLE backslash-then-r with a CR return
}

## ================ DONE SANITY-CHECKING A BUNCH OF VARIABLES ==================

sub openSmartAndGetFilehandle($) {
    # returns a FILEHANDLE. Can be standard in, if the 'filename' was specified as '-'
    # Transparently handles ".gz" and ".bz2" files.
    # This is the MARCH 6, 2013 version of this function.
    my ($filename) = @_;
    if ($filename eq '-') { return(*STDIN); } # <-- RETURN!!!
    my $reader;
    if    ($filename =~ /[.]gz$/)  { $reader = "zcat $filename |"; }     # Un-gzip a file and send it to STDOUT.
    elsif ($filename =~ /[.]bz2$/) { $reader = "bzcat $filename |"; }    # Un-bz2 a file and send it to STDOUT
    elsif ($filename =~ /[.]zip$/) { $reader = "unzip -p $filename |"; } # Un-regular-zip a file and send it to STDOUT with "-p": which is DIFFERENT from -c (-c is NOT what you want here). See 'man unzip'
    else                           { $reader = "$filename"; }  # Default: just read a file normally
    my $fh;
    open($fh, "$reader") or die("Couldn't read from <$filename>: $!");
    return $fh;
}

sub readIntoHash($$$$$) {
    my ($filename, $theDelim, $keyIndexCountingFromOne, $masterHashRef, $origCaseHashRef) = @_;
    my $numDupeKeys = 0;
    my $lineNum = 1;
    my $theFileHandle = openSmartAndGetFilehandle($filename);
    foreach my $line ( <$theFileHandle> ) {
	$lineNum++;
	if ($line =~ m/\r/) {
	    verboseWarnPrint("WARNING: The file <$filename> appears to have either WINDOWS-STYLE line endings or MAC-STYLE line endings ( with an '\\r' character) (as seen on line $lineNum)!\n         We are automatically REMOVING this malevolent character from the line.");
	    $line =~ s/(\r\n|\r)/\n/; # Turn PC-style \r\n, or Mac-style just-plain-\r into UNIX \n
	}
	chomp($line);
	#if(/\S/) { ## if there's some content that's non-spaces-only
	my @sp1 = split($theDelim, $line, $SPLIT_WITH_TRAILING_DELIMS);
	my $theKey = $sp1[ ($keyIndexCountingFromOne - 1) ]; # index from ZERO here!
	if ($theKey =~ /\s/) {
	    verboseWarnPrint("Warning: the key \"$theKey\" on line $lineNum of file 2 ($filename) had whitespace in it. This MAY be unintentional--beware!");
	}
	if (defined($masterHashRef->{$theKey})) {
	    ($numDupeKeys < $MAX_DUPE_KEYS_TO_REPORT) and verboseWarnPrint("Warning: the key <$theKey> appeared more than once in <$filename> (on line $lineNum). We are only keeping the FIRST instance of this key.");
	    ($numDupeKeys == $MAX_DUPE_KEYS_TO_REPORT) and verboseWarnPrint("Warning: suppressing any future duplicate key warnings.");
	    $numDupeKeys++;
	} else {
	    # Found a UNIQUE new key! ($isDebugging) && print STDERR "Added a line for the key <$theKey>.\n";
	    if ((0 == length($theKey)) and !$allowEmptyKey) {
		verboseWarnPrint("Warning: skipping an empty key on line <$lineNum>!");
	    } else {
		# Key was valid, OR we are allowing empty keys!
		if (defined($origCaseHashRef)) { $origCaseHashRef->{uc($theKey)} = $theKey; } # maps from the UPPER-CASE version of this key (KEY) back to the one we ACTUALLY put in the hash (VALUE)
		@{$masterHashRef->{$theKey}} = @sp1; # whole SPLIT UP line, even the key, goes into the hash!!!. # masterHashRef is a hash of ARRAYS: each line is ALREADY SPLIT UP by its delimiter
	    }
	}
    }
    ($numDupeKeys > 0) and verboseWarnPrint("Warning: $numDupeKeys duplicate keys were skipped in <$filename>.");
    if ($filename ne '-') { close($theFileHandle); } # close the file we opened in 'openSmartAndGetFilehandle' . This may not actually be necessary
}

#my %hash1 = readIntoHash($file1  , $delim1, $keyCol1);
#($isDebugging) && print STDERR ("Read in this many keys: " . scalar(keys(%hash1)) . " from primary file.\n");

sub arrayOfNonKeyElements(\@$) {
    # Returns everything EXCEPT the key! This is because by default, when joining, you move the key to the FRONT of the line, and then do not print it again later on the line.
    my ($inputArrayPtr, $inputKey) = @_;
    ($isDebugging) && (($inputKey >= 1) or die "Whoa, the input key was LESS THAN ONE, which is impossible, since numbering for keys starts from 1! Not zero-indexing!!!");
    my @nonKeyElements = (); # the final array with everything BUT the key. Apparently doesn't matter much whether we pre-allocate it to the right size or not, speed-wise.
    for (my $i = 0; $i < scalar(@{$inputArrayPtr}); $i++) {
	if ($i != ($inputKey-1)) {
	    # remember that the input key is indexed from ONE and not ZERO!!!
	    push(@nonKeyElements, $inputArrayPtr->[$i]); # It's not a key, so add it to the array
	    #debugPrint("Adding $i (index $i is not equal to the input key index $inputKey)...\n");
	} else {
	    # Huh, this IS a key item, so don't include it!!!
	    #debugPrint("OMITTING $i (index $i is EXACTLY EQUAL to the input key index $inputKey. Remember that one of them counts from zero!)...\n");
	}
    }
    return (@nonKeyElements); # The subset of the inputArrayPtr that does NOT include the non-key elements!
}

sub joinedUpOutputLine($$$$$$$) {
    my ($delim, $mainKey, $array1Ref, $k1col, $array2Ref, $k2col, $shouldNotMoveKey) = @_;
    if ($shouldNotMoveKey) {
	# do NOT move the key to the front of the line---this is a bit unusual! The key stays wherever it was on the line.
	return join($delim, @$array1Ref, @$array2Ref); # no newline!
    } else {
	# key gets moved to the front! -- like in unix join. This is the DEFAULT and UNIX-join-like way of doing it
	if (@{$array2Ref}) {
	    return join($delim, $mainKey, arrayOfNonKeyElements(@{$array1Ref}, $k1col), arrayOfNonKeyElements(@{$array2Ref}, $k2col)); # no newline!
	} else {
	    # array 2 was EMPTY
	    return join($delim, $mainKey, arrayOfNonKeyElements(@{$array1Ref}, $k1col)); # no newline!
	}
    }
}

my %hash2 = ();
my %originalCaseHash = (); # Hash: key = UPPER CASE version of key, value = ORIGINAL version of key
my $originalCaseHashRef = ($shouldIgnoreCase) ? \%originalCaseHash : undef; # UNDEFINED if we aren't ignoring case
readIntoHash($file2, $delim2, $keyCol2, \%hash2, $originalCaseHashRef);
debugPrint("Read in this many keys: " . scalar(keys(%hash2)) . " from secondary file.\n");

my $lineNumPrimary  = 0;
my $primaryFH       = openSmartAndGetFilehandle($file1);
my $prevLineCount1  = undef;
my $prevLineCount2  = undef;
my $numWeirdLengths = 0;
foreach my $line (<$primaryFH>) {
    $lineNumPrimary++; # Start it at ONE during the first iteration of the loop! (Was initialized to zero before!)
    if ($line =~ m/\r/) {
	    verboseWarnPrint("WARNING: The file <$file1> appears to have either WINDOWS-STYLE line endings or MAC-STYLE line endings ( with an '\\r' character) (as seen on line $lineNumPrimary).\n         We are automatically REMOVING this malevolent character from the line.");
	    $line =~ s/(\r\n|\r)/\n/; # Turn PC-style \r\n, or Mac-style just-plain-\r into UNIX \n
    }
    chomp($line); # Chomp each line of line endings no matter what. Even the header line!
    if ($lineNumPrimary <= $numHeaderLines) { # This is still a HEADER line, also: lineNumPrimary starts at 1, so this should be '<=' and not '<' to work properly!
	verboseWarnPrint("Note: directly printing $lineNumPrimary of $numHeaderLines header line(s) from file 1 (\"$file1\")...");
	print STDOUT $line . "\n"; # Print the input line, making sure to use a '\n' as the ending.
	next; # <-- skip to next iteration of loop!
    }
    #if(/\S/) { ## if there's some content that's non-spaces-only
    my @sp1 = split($delim1, $line, $SPLIT_WITH_TRAILING_DELIMS); # split-up line
    if (defined($prevLineCount1) and $prevLineCount1 != scalar(@sp1)) {
	($numWeirdLengths < $MAX_WEIRD_LINE_LENGTHS_TO_REPORT) and verboseWarnPrint("Warning: the number of elements in file 1 ($file1) is not constant. Line $lineNumPrimary had this many elements: " . scalar(@sp1) . " (previous line had $prevLineCount1)");
	($numWeirdLengths == $MAX_WEIRD_LINE_LENGTHS_TO_REPORT) and verboseWarnPrint("Warning: suppressing any further non-constant elements-per-line warnings.");
	$numWeirdLengths++;
    }
    $prevLineCount1 = scalar(@sp1);

    if ((($keyCol1-1) >= scalar(@sp1))) {
	# we're going to SKIP the line entirely if there was no key AT ALL (not even something blank!) at this location
	verboseWarnPrint("Warning: skipping line $lineNumPrimary---there was no key column on that line. (Key column: $keyCol1, Columns on line: " . scalar(@sp1) . ")");
	next; # <-- skip to next iteration of loop!ll
    }
    
    my $thisKey = $sp1[ ($keyCol1-1) ]; # index from ZERO here, that's why we subtract 1 from the key column
    if ($thisKey =~ /\s/) {
	verboseWarnPrint("Warning: key \"$thisKey\" on line $lineNumPrimary of file 1 ($file1) had whitespace in it. This is possibly unintentional--beware!");
    }
    if (!defined($thisKey)) {
	verboseWarnPrint("Warning: key at column number $keyCol1 (counting from 1) on line $lineNumPrimary of file 1 ($file1) was UNDEFINED. This is possibly a programming error in join.pl!");
    }
    
    my @sp2; # matching split-up line
    if ($shouldIgnoreCase) {
	my $keyInOrigCase = $originalCaseHash{uc($thisKey)}; # mutate the key so that it's in the SAME CASE as it was in the key we added
	#print "Found key \"$keyInOrigCase\" from upper-case " . uc($thisKey) . "...\n";
	if (defined($hash2{$thisKey})) {
	    @sp2 = @{$hash2{$thisKey}}; # () <-- empty list/array is the result of "didn't find anything"
	} else {
	    @sp2 = (defined($keyInOrigCase) and defined($hash2{$keyInOrigCase})) ? @{$hash2{$keyInOrigCase}} : (); # () <-- empty list/array is the result of "didn't find anything"
	}
    } else {
	@sp2 = (defined($hash2{$thisKey})) ? @{$hash2{$thisKey}} : (); # () <-- empty list/array is the result of "didn't find anything"
    }
    
    if (@sp2) {
	# Got a match for the key in question!
	if (defined($prevLineCount2) and $prevLineCount2 != scalar(@sp2)) {
	    ($numWeirdLengths < $MAX_WEIRD_LINE_LENGTHS_TO_REPORT) and verboseWarnPrint("Warning: the number of elements in file 2 ($file2) is not constant. Got a line with this many elements: " . scalar(@sp2) . " (previous line had $prevLineCount2)");
	    ($numWeirdLengths == $MAX_WEIRD_LINE_LENGTHS_TO_REPORT) and verboseWarnPrint("Warning: suppressing any further non-constant elements-per-line warnings.");
	    $numWeirdLengths++;
	}
	$prevLineCount2 = scalar(@sp2);
	
	if ($shouldNegate) { 
	    # Since we are NEGATING this, don't print the match when it's found (only when it isn't...)
	} else {
	    # Great, the OTHER file had a valid entry for this key as well! So print it... UNLESS we are negating.
	    print STDOUT joinedUpOutputLine($outputDelim, $thisKey, \@sp1, $keyCol1, \@sp2, $keyCol2, $shouldKeepKeyInOriginalPosition) . "\n";
	}
    } else {
	# Ok, there was NO MATCH for this key!
	debugPrint("WARNING: Hash2 didn't have the key $thisKey\n");
	if ($shouldNegate) {
	    # We didn't find a match for this key, but because we are NEGATING the output, we'll print this line anyway
	    print STDOUT joinedUpOutputLine($outputDelim, $thisKey, \@sp1, $keyCol1, \@sp2, $keyCol2, $shouldKeepKeyInOriginalPosition) . "\n";
	} else {
	    if (defined($stringWhenNoMatch)) {
		# We print the line ANYWAY, because the user specified an outer join, with the "-o SOMETHING" option.
		my $suffixWhenNoMatch = (length($stringWhenNoMatch)>0) ? "${outputDelim}${stringWhenNoMatch}" : "$stringWhenNoMatch"; # handle zero-length -ob SPECIALLY
		print STDOUT joinedUpOutputLine($outputDelim, $thisKey, \@sp1, $keyCol1, \@sp2, $keyCol2, $shouldKeepKeyInOriginalPosition) . $suffixWhenNoMatch . "\n";
		#print STDOUT join($outputDelim, $thisKey, arrayOfNonKeyElements(@sp1, $keyCol1)) . $suffixWhenNoMatch . "\n";
	    } else {
		# Omit the line entirely, since there was no match in the secondary file.
	    }
	}
    }
    #print "$line\n";
}
if ($file1 ne '-') { close($primaryFH); } # close the file we opened in 'openSmartAndGetFilehandle'

exit(0); # looks like we were successful


################# END MAIN #############################

__DATA__
syntax: join.pl [OPTIONS] LOOKUP_FILE  HUGE_DICTIONARY_FILE

join.pl goes through each line/key in the LOOKUP_FILE, and finds the *first* matching
key in the DICTIONARY_FILE. The data from those corresponding rows is then
printed out. It does not handle cross-products.

Unlike the UNIX "join", join.pl does NOT require sorted keys.

 * Note that join.pl reads the ENTIRE contents of the second file into memory! It may be unsuitable for joining very large (> 1000 MB) files.

DESCRIPTION:

This script takes two tables, contained in delimited files, as input and
produces a new table that is a join of FILE1 and FILE2.

By default, files are assumed to be tab-delimited with the keys in the first column.
(See the options to change these defaults.)

If FILE1 contains the tuple (VOWELS, A, O, A) and FILE2 contains the
tuple (VOWELS, I, U, U) then the joined output will be (VOWELS, A, O, A, I, U, U)

CAVEATS:

Every line of the LOOKUP_FILE is processed in the
order that it appears in the file.

The DICTIONARY_FILE is only for handling join operations.

If the DICTIONARY_FILE contains several lines with the same key,
only the *LAST* key read will actually ever be used.

Because of this, join.pl does NOT exactly duplicate the function of
GNU "join"--in particular it does not output the cross-product
in multiple-key situations.

*** NOTE: If you are trying to merge two files, or you are joining ***
*** multiple times with multiple files, try "join_multi.pl" .      ***

OPTIONS are:

-1 COLUMN: Include column COL from FILE1 as part of the key (default: 1).
        Only supports ONE key field, unlike unix join. Indexed starting at 1.

-2 COLUMN: Include column COL from FILE2 as part of the key (default: 1).
        Only supports ONE key field, unlike unix join. Indexed starting at 1.

-o TEXT:  Do a left outer join.  If a key in FILE1 is not in FILE2, then the
          tuple from FILE1 is printed along with the FILLER TEXT in place of
          a tuple from FILE2 (by default these tuples are not reported in the
          result). (See 'examples' section for usage.)
          Example:  join -o "NO_match_in_file2" FILE1.txt FILE2.txt > OUTPUT.txt

--ob: Same as -o ''---report blank entries for non-matching output.

-h INT: Number of header lines to print verbatim (without joining) from FILE1. (default: 0)

-neg: Negate output -- print keys that are in FILE1 but not in FILE2.
        These keys are the same ones that would be left out of the join,
        or those that would have a FILL tuple in a left outer join
        Cannot specify both this AND ALSO -ob or -o.

-t DELIM or -d DELIM or --delim=DELIM:
        Set the input delimiters for both FILE1 and FILE2 to DELIM (default: tab)
        Equivalent to setting both --d1 and --d2.

--d1=DELIM: Set the input delimiter for FILE1 to DELIM (default: tab).

--d2=DELIM: Set the input delimiter for FILE2 to DELIM (default: tab).

--do=DELIM: Set the OUTPUT delimiter to DELIM. (default: same as input delim)

--no-reorder-key or --nrk: (Default: DO move the key to column 1 -- the UNIX join default)
  If this is specified, then instead of having the key at the BEGINNING of the output line (the default),
  the key stays wherever it was (like if it was in column 10, it STAYS in column 10. Regular join would move
  it to column 1).

--allow-empty-key or --eok: (Default: do not allow it): Whether to allow an EMPTY key as valid.
  Note: even when the empty key is NOT allowed for matching, if we do an outer join or
  negation, we will still print items from FILE1 where the key was blank.

-i or --ignore-case: (Case-insensitive join)

-q or --quiet : No verbose output. May hide some useful warning messages!

Example:

join.pl -1 1 -2 2 file1--key_in_first_col.txt  file2--key_in_second_col.compressed.gz > join.output.txt
  Print lines from file1 that are ALSO in file2, and append the data from file2 in the output.
  This is the most standard-plain-vanilla join, and should be similar to the unix "join" results.

join.pl -o "NO_MATCH_HERE_I_SEE" -1 4 -2 1 file_with_key_in_fourth_col.compressed.bz2  file2.gz > join.with.unmatched.rows.txt
  Also print the un-matching lines from the first file (lines with no match will say "NO_MATCH_HERE_I_SEE")

cat myfile.txt | join.pl - b.txt | less -S
  Read from STDIN (use a '-' instead of a filename!), and pipe into the program "less".

join.pl --no-reorder-key a.txt b.txt > a_and_b.txt
join.pl --no-reorder-key a_and_b.txt c.txt > a_and_b_and_c.txt
  Since we are using --no-reorder-key, it makes it easier to join multiple files sometimes, since the order of columns
  does not keep moving around. This is only of interest if your key column is NOT the very first one already!

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

join.pl -o "NOPE!" File2 File1
Results in:
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

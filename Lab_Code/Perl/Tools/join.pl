#!/usr/bin/perl -w

# New version of Join.pl by Alex Williams. (This isn't related to the previous UCSC code at all... and probably produces different results! Note that both versions produce different results from UNIX join, even on properly sorted input!)

use strict; use warnings; use Getopt::Long;
$| = 1;  # Flush output to STDOUT immediately.

my $isDebugging = 0;

my $verbose = 1; # use 'q' (quiet) to suppress this

sub quitWithUsageError($) { print($_[0] . "\n"); printUsageAndQuit(); print($_[0] . "\n"); }
sub printUsageAndQuit() { printUsageAndContinue(); exit(1); }
sub printUsageAndContinue() {    print STDOUT <DATA>; }
sub debugPrint($) {   ($isDebugging) and print STDERR $_[0]; }
sub verboseWarnPrint($) {
    if ($verbose) { print STDERR $_[0] . "\n"; }
}
my $keyCol1 = 1; # indexed from ONE rather than 0!
my $keyCol2 = 1; # indexed from ONE rather than 0!

my $DEFAULT_DELIM_IF_NOTHING_ELSE_IS_SPECIFIED = "\t";

my $delim1 = $DEFAULT_DELIM_IF_NOTHING_ELSE_IS_SPECIFIED;
my $delim2 = $DEFAULT_DELIM_IF_NOTHING_ELSE_IS_SPECIFIED;
my $delimBoth = undef;
my $outputDelim = undef;

my $file1   = undef;
my $file2 = undef;

my $shouldNegate = 0; # whether we should NEGATE the output
my $shouldIgnoreCase = 0; # by default, case-sensitive

my $stringWhenNoMatch = undef;

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
	   , "i|ignore-case!" => \$shouldIgnoreCase
	   , "debug!" => \$isDebugging
    ) or printUsageAndQuit();

my $numUnprocessedArgs = scalar(@ARGV);
if ($numUnprocessedArgs != 2) {
    quitWithUsageError("Error in arguments! You must send exactly TWO filenames (or one filename and '-' for STDIN) to this program.");
}
$file1 = $ARGV[0]; # first un-processed join argument
$file2 = $ARGV[1]; # second un-processed join argument

foreach my $ff ($file1, $file2) {
    (($ff eq '-') or (-f $ff)) or die "File $ff did not exist!"; # Specified files must be either - (for stdin) or a real, valid, existing filename)
}

if ($shouldNegate && defined($stringWhenNoMatch)) {
    quitWithUsageError("Error in arguments! It doesn't make sense to both --neg (negate) the join AND ALSO specify -o or --ob -- the outer join specifies that we should print lines REGARDLESS of match, whereas the --neg specifies that we should ONLY print lines with no match. You cannot specifiy both of these options at the same time.");
}


if (defined($delimBoth)) { # If "delimBoth" was specified, then set both of the input delimiters accordingly.
    $delim1 = $delimBoth; $delim2 = $delimBoth;
}

if (!defined($outputDelim)) { # Figure out what the output delimiter should be, if it wasn't explicitly specified.
    if (defined($delimBoth)) { $outputDelim = $delimBoth; } # default: set the output delim to whatever the input delim was
    elsif ($delim1 eq $delim2) { $outputDelim = $delim1; } # or we can set it to the manually-specified delimiters, if they are the SAME only
    else { $outputDelim = $DEFAULT_DELIM_IF_NOTHING_ELSE_IS_SPECIFIED; } # otherwise, set it to the default delimiter
}

($keyCol1 != 0) or die "Key1 CANNOT BE ZERO! These indices are numbered from ONE and not zero!";
($keyCol2 != 0) or die "Key2 CANNOT BE ZERO! These indices are numbered from ONE and not zero!";

sub openSmartAndGetFilehandle($) {
    # returns a FILEHANDLE. Can be standard in, if the 'filename' was specified as '-'
    my ($filename) = @_;
    if ($filename eq '-') {
	return(*STDIN);
    } else {
	my $fh;
	open($fh, '<', $filename) or die("Couldn't read from file: $!");
	return $fh;
    }
}

sub readIntoHash($$$$$) {
    my ($filename, $theDelim, $keyIndexCountingFromOne, $masterHashRef, $uppercaseHashMapRef) = @_;
    my $numDupeKeys = 0;
    my $lineNum = 1;
    my $theFileHandle = openSmartAndGetFilehandle($filename);
    foreach my $line ( <$theFileHandle> ) {
	chomp($line);
	#if(/\S/) { ## if there's some content that's non-spaces-only
	my @sp1 = split($theDelim, $line);
	my $thisKey = $sp1[ ($keyIndexCountingFromOne - 1) ]; # index from ZERO here!
	if (exists($masterHashRef->{$thisKey})) {
	    print STDERR "Warning: the key <$thisKey> appeared more than once in <$filename> (on line $lineNum). We are only keeping the FIRST instance of this key.\n";
	    $numDupeKeys++;
	} else {
	    # Found a UNIQUE new key!
	    # ($isDebugging) && print STDERR "Added a line for the key <$thisKey>.\n";
	    if (defined($uppercaseHashMapRef)) { # apparently we want to deal with things case-insensitively
		$uppercaseHashMapRef->{uc($thisKey)} = $thisKey; # maps from the UPPER-CASE version of this key back to the one we ACTUALLY put in the hash
	    }
	    # masterHashRef is a hash of ARRAYS: each line is ALREADY SPLIT UP by its delimiter
	    @{$masterHashRef->{$thisKey}} = @sp1; # whole SPLIT UP line, even the key, goes into the hash!!!
	    #print "$line\n";
	}
	$lineNum++;
    }

    if ($numDupeKeys > 0) { print STDERR "Warning: $numDupeKeys duplicate keys were skipped in <$filename>.\n"; }
    if ($filename ne '-') { close($theFileHandle); } # close the file we opened in 'openSmartAndGetFilehandle' . This may not actually be necessary
}


#my %hash1 = readIntoHash($file1  , $delim1, $keyCol1);
#($isDebugging) && print STDERR ("Read in this many keys: " . scalar(keys(%hash1)) . " from primary file.\n");

sub arrayOfNonKeyElements(\@$) {
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

my %hash2 = ();
my %uppercaseHash = ();
my $uppercaseHashRef = ($shouldIgnoreCase) ? \%uppercaseHash : undef; # UNDEFINED if we aren't ignoring case
readIntoHash($file2, $delim2, $keyCol2, \%hash2, $uppercaseHashRef);
debugPrint("Read in this many keys: " . scalar(keys(%hash2)) . " from secondary file.\n");

my $lineNumPrimary = 1;
my $primaryFH = openSmartAndGetFilehandle($file1);
my $numElementsOnPreviousLineInPrimary = undef;
my $numElementsOnPreviousLineInSecondary = undef;
foreach my $line (<$primaryFH>) {
    chomp($line);
    #if(/\S/) { ## if there's some content that's non-spaces-only
    my @sp1 = split($delim1, $line); # split-up line
    my $thisKey = $sp1[ ($keyCol1-1) ]; # index from ZERO here, that's why we subtract 1 from the key column
    my @sp2; # matching split-up line

    if ($shouldIgnoreCase) {
	my $keyInSameCaseItWasInTheOriginalHash = $uppercaseHash{uc($thisKey)}; # mutate the key so that it's in the SAME CASE as it was in the key we added
	@sp2 = (exists($hash2{$keyInSameCaseItWasInTheOriginalHash})) ? @{$hash2{$keyInSameCaseItWasInTheOriginalHash}} : ();
    } else {
	@sp2 = (exists($hash2{$thisKey})) ? @{$hash2{$thisKey}} : ();
    }

    if (@sp2) {
	# Got a match for the key in question!
	if (defined($numElementsOnPreviousLineInSecondary) && $numElementsOnPreviousLineInSecondary != scalar(@sp2)) {
	    verboseWarnPrint("Warning: the number of elements in <$file2> is not constant. Got a line with this many elements: " . scalar(@sp2) . "\n");	    
	}
	$numElementsOnPreviousLineInSecondary = scalar(@sp2);
	
	if ($shouldNegate) { 
	    # Since we are NEGATING this, don't print the match when it's found (only when it isn't...)
	} else {
	    # Great, the OTHER file had a valid entry for this key as well! So print it... UNLESS we are negating.
	    print STDOUT join($outputDelim, $thisKey, arrayOfNonKeyElements(@sp1, $keyCol1), arrayOfNonKeyElements(@sp2, $keyCol2)) . "\n";
	}
    } else {
	# Ok, there was NO MATCH for this key!
	debugPrint("Hash2 didn't have the key $thisKey\n");
	if ($shouldNegate) {
	    # We didn't find a match for this key, but because we are NEGATING the output, we'll print this line anyway
	    print STDOUT join($outputDelim, $thisKey, @sp1[getAllNonKeyIndices(@sp1, $keyCol1)]) . "\n";
	} else {
	    if (defined($stringWhenNoMatch)) {
		# We print the line ANYWAY, because the user specified an outer join, with the "-o SOMETHING" option.
		my $suffixWhenNoMatch = (length($stringWhenNoMatch)>0) ? "${outputDelim}${stringWhenNoMatch}" : "$stringWhenNoMatch"; # handle zero-length -ob SPECIALLY
		print STDOUT join($outputDelim, $thisKey, getAllNonKeyIndices(@sp1, $keyCol1)) . $suffixWhenNoMatch . "\n";
	    } else {
		# Omit the line entirely, since there was no match in the secondary file.
	    }
	}
    }
    #print "$line\n";

    if (defined($numElementsOnPreviousLineInSecondary) && $numElementsOnPreviousLineInSecondary != scalar(@sp1)) {
	verboseWarnPrint("Warning: the number of elements in <$file1> is not constant. Got a line with this many elements: " . scalar(@sp1) . "\n");
    }
    $numElementsOnPreviousLineInPrimary = scalar(@sp1);
    $lineNumPrimary++;
}
if ($file1 ne '-') { close($primaryFH); } # close the file we opened in 'openSmartAndGetFilehandle'

exit(0); # looks like we were successful


################# END MAIN #############################

__DATA__
syntax: join.pl [OPTIONS] LOOKUP_FILE  DICTIONARY_FILE

join.pl goes through each line/key in the LOOKUP_FILE, and finds the *first* matching
key in the DICTIONARY_FILE. The data from those corresponding rows is then
printed out. It does not handle cross-products.

Unlike the UNIX "join", join.pl does NOT require sorted keys.

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

-1 COL: Include column COL from FILE1 as part of the key (default is 1).
        Only supports ONE key field.

-2 COL: Include column COL from FILE2 as part of the key (default is 1).

-o FILLER: Do a left outer join.  If a key in FILE1 is not in FILE2, then the
          tuple from FILE1 is printed along with the text in FILLER in place of
          a tuple from FILE2 (by default these tuples are not reported in the
          result).  See -of option also to supply FILLER from a file.
          (See below for an example of usage.)

-ob: Same as -o but do not actually fill with anything (identical to -o '').

-neg: Negate output -- print keys that are in FILE1 but not in FILE2.
        These keys are the same ones that would be left out of the join,
        or those that would have a FILL tuple in a left outer join
        Cannot specify both this AND ALSO -ob or -o.


-t DELIM or -d DELIM or --delim=DELIM: Set the input delimiters for both FILE1 and FILE2 to DELIM (default: tab)
          Equivalent to setting both --d1 and --d2.

--d1 DELIM: Set the input delimiter for FILE1 to DELIM (default: tab).

--d2 DELIM: Set the input delimiter for FILE2 to DELIM (default: tab).

-do DELIM: Set the OUTPUT delimiter to DELIM. (default: same as input delim)

-i or --ignore-case: (Case-insensitive join)

Example:

If you have the following two files:

File 1: (tab-delimited)
 AAA   avar    aard
 ZZZ   zebra
 BBB   beta    bead    been
 MMM   man     most

File 2: (tab-delimited)
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
...is:
 AAA    avar     aard      111
 ZZZ    zebra    NONE!
 BBB    beta     bead      been      222
 MMM    man      mouse     NONE!

Note that if you switch the order of File1 and File2 for an outer join...
   join.pl -o "NOPE!" File2 File1
...you get a different result (the first file specifies the keys):
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

#!/usr/bin/perl -w

use strict; use warnings; use Getopt::Long;
$| = 1;  # Flush output to STDOUT immediately.

my $isDebugging = 0;

sub quitWithUsageError($) { print($_[0] . "\n"); printUsageAndQuit(); print($_[0] . "\n"); }
sub printUsageAndQuit() { printUsageAndContinue(); exit(1); }
sub printUsageAndContinue() {    print STDOUT <DATA>; }
sub debugPrint($) {   ($isDebugging) and print STDERR $_[0]; }

my $key1 = 1; # indexed from ONE rather than 0!
my $key2 = 1; # indexed from ONE rather than 0!

my $delim1 = "\t";
my $delim2 = "\t";
my $delimBoth = undef;
my $outputDelim = "\t";

my $filePrimary   = undef;
my $fileSecondary = undef;

my $stringWhenNoMatch = undef;

$Getopt::Long::passthrough = 1; # ignore arguments we don't recognize in GetOptions, and put them in @ARGV
GetOptions("help|?|man" => sub { printUsageAndQuit(); }
	   , "1=s" => \$key1
	   , "2=s" => \$key2
	   , "d1=s" => \$delim1
	   , "d2=s" => \$delim2
	   , "o=s"  => \$stringWhenNoMatch
	   , "ob!"  => sub { $stringWhenNoMatch = ''; } # shortcut for just a blank when there's no match. Default is to OMIT lines with no match.
	   , "d=s" => \$delimBoth
	   , "do=s" => \$outputDelim
	   , "debug!" => \$isDebugging
    ) or printUsageAndQuit();

if (defined($delimBoth)) {
    $delim1 = $delimBoth; $delim2 = $delimBoth;
}

my $numUnprocessedArgs = scalar(@ARGV);
if ($numUnprocessedArgs != 2) {
    quitWithUsageError("Error in arguments! You must send exactly TWO filenames to this program.\n");
}

$filePrimary = $ARGV[0]; # first un-processed join argument
$fileSecondary = $ARGV[1]; # second un-processed join argument

($key1 != 0) or die "Key1 CANNOT BE ZERO! These indices are numbered from ONE and not zero!";
($key2 != 0) or die "Key2 CANNOT BE ZERO! These indices are numbered from ONE and not zero!";

foreach my $ff ($filePrimary, $fileSecondary) {
    # Specified files must be either - (for stdin) or a real, valid, existing filename)
    (($ff eq '-') or (-f $ff)) or die "File $ff did not exist!";
}

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

sub readIntoHash($$$) {
    my ($filename, $theDelim, $keyIndexCountingFromOne) = @_;
    my %tempHash = ();
    my $numDupeKeys = 0;
    my $lineNum = 1;
    my $theFileHandle = openSmartAndGetFilehandle($filename);
    foreach my $line ( <$theFileHandle> ) {
	chomp($line);
	#if(/\S/) { ## if there's some content that's non-spaces-only
	my @sp = split($theDelim, $line);
	my $thisKey = $sp[ ($keyIndexCountingFromOne - 1) ]; # index from ZERO here!
	if (exists($tempHash{$thisKey})) {
	    print STDERR "Warning: the key <$thisKey> appeared more than once in <$filename> (on line $lineNum). We are only keeping the FIRST instance of this key.\n";
	    $numDupeKeys++;
	} else {
	    # Found a UNIQUE new key!
	    # ($isDebugging) && print STDERR "Added a line for the key <$thisKey>.\n";
	    @{$tempHash{$thisKey}} = @sp; # whole SPLIT UP line, even the key!!!
	    #print "$line\n";
	}
	$lineNum++;
    }

    if ($numDupeKeys > 0) { print STDERR "Warning: $numDupeKeys duplicate keys were skipped in <$filename>.\n"; }
    if ($filename ne '-') { close($theFileHandle); } # close the file we opened in 'openSmartAndGetFilehandle'
    return %tempHash; # This is a hash of ARRAYS: each line is ALREADY SPLIT UP by its delimiter
}


#my %hash1 = readIntoHash($filePrimary  , $delim1, $key1);
#($isDebugging) && print STDERR ("Read in this many keys: " . scalar(keys(%hash1)) . " from primary file.\n");

sub getAllNonKeyIndices(\@$) {
    my ($inputArrayPtr, $inputKey) = @_;

    ($isDebugging) && (($inputKey >= 1) or die "Whoa, the input key was LESS THAN ONE, which is impossible, since numbering for keys starts from 1! Not zero-indexing!!!");
    my @nonKeyIndices = ();
    for (my $i = 0; $i < scalar(@{$inputArrayPtr}); $i++) {
	if ($i != ($inputKey-1)) {
	    # remember that the input key is indexed from ONE and not ZERO!!!
	    push(@nonKeyIndices, $i); #It's not a key, so add it to the array
	    debugPrint("Adding $i (index $i is not equal to the input key index $inputKey)...\n");
	} else {
	    debugPrint("OMITTING $i (index $i is EXACTLY EQUAL to the input key index $inputKey. Remember that one of them counts from zero!)...\n");
	}
    }
    return (@nonKeyIndices);
}

my %hash2 = readIntoHash($fileSecondary, $delim2, $key2);
debugPrint("Read in this many keys: " . scalar(keys(%hash2)) . " from secondary file.\n");

my $lineNumPrimary = 1;
my $primaryFH = openSmartAndGetFilehandle($filePrimary);
foreach my $line (<$primaryFH>) {
    chomp($line);
    #if(/\S/) { ## if there's some content that's non-spaces-only
    my @sp = split($delim1, $line);
    my $thisKey = $sp[ ($key1-1) ]; # index from ZERO here

    my @nonKeyIndicesPrimary = getAllNonKeyIndices(@sp, $key1);
    #($isDebugging) && print STDERR "Looking for the key <$thisKey>...\n";
    my @matchingSp = (exists($hash2{$thisKey})) ? @{$hash2{$thisKey}} : undef;
    if (defined(@matchingSp)) {
	# Great, the OTHER file had a valid entry for this key as well!
	my @nonKeyIndicesSecondary = getAllNonKeyIndices(@matchingSp, $key2);
	print STDOUT join($outputDelim, @sp[@nonKeyIndicesPrimary], @matchingSp[@nonKeyIndicesSecondary]) . "\n";
    } else {
	debugPrint("Hash2 didn't have the key $thisKey\n");
	if (defined($stringWhenNoMatch)) {
	    my $suffixWhenNoMatch = (length($stringWhenNoMatch)>0) ? "${outputDelim}${stringWhenNoMatch}" : "$stringWhenNoMatch"; # handle zero-length -ob SPECIALLY
	    print STDOUT join($outputDelim, @sp[@nonKeyIndicesPrimary]) . $suffixWhenNoMatch . "\n";
	} else {
	    # omit the line entirely, since there was no match in the secondary file
	}
    }
    #print "$line\n";
    $lineNumPrimary++;
}
if ($filePrimary ne '-') { close($primaryFH); } # close the file we opened in 'openSmartAndGetFilehandle'


################# BEGIN MAIN ###########################

# my (@key1)                                = (1);
# my (@key2)                                = (1);
# my ($beg,$end)                            = (0,0);
# my ($file1,$file2)                        = ('','');
# my ($key)                                 = '';
# my (%values)                              = ();
# my (%exists)                              = ();
# my $delim                                 = "\t";
# my ($delim_in1, $delim_in2)               = ("\t","\t");
# my ($delim_out)                           = ("\t");
# my ($value1, $value2)                     = ('','');
# my ($printable, $printable1, $printable2) = ('','','');
# my ($suppress1, $suppress2, $suppressk) = (0,0,0);
# my $negate=0;
# my ($numeric) = 0;
# my $empty = '';
# my $max_tuple_size = undef;
# my $skip_empty_lines = 0;
# my $fill = '';
# my $empty_placeholder = '___@@@_THIS_IS_A_BLANK_VALUE_@@@___';
# my ($outer) = 0;
# my ($reverse) = 0;
# my ($uppercase) = 0;
# my $hit = 0;
# my @fill_lines = ();
# my $merge=0;
# my $verbose=1;
# my $header=0;


# while(@ARGV) {
#     my $arg = shift @ARGV;
#     if($arg eq '--help') {
# 	print STDOUT <DATA>;
# 	exit(0);
#     } elsif($arg eq '-q') {
# 	$verbose = 0;
#     } elsif($arg eq '-f') {
# 	$arg = shift @ARGV;
# 	@key1 = &parseRanges($arg);
# 	@key2 = @key1;
#     }
#     elsif($arg eq '-1') {
#         $arg = shift @ARGV;
#         @key1 = &parseRanges($arg);
#     }
#     elsif($arg eq '-2') {
#         $arg = shift @ARGV;
#         @key2 = &parseRanges($arg);
#     } elsif($arg eq '-o' or $arg eq '-of') {
# 	$arg = shift @ARGV;
# 	if(-f $arg) { ## if it's a file...
# 	    my $FILE = openFile($arg);
# 	    (defined(openFile($arg))) or die "join.pl: Could not open file '$arg' to find outer text, skipping.";
# 	    while(<$FILE>) {
# 		if(/\S/) { ## if there's some content that's non-spaces-only
# 		    chomp;
# 		    push(@fill_lines, $_);
# 		}
# 	    }
# 	    close($FILE);
# 	    $outer = 1;
# 	} else {
# 	    $fill = $arg; ## fill with this text
# 	    $fill =~ s/[\\]t/\t/g; ## any time you see "slash t" replace it with an actual tab!
# 	    $outer = 1;
# 	}
#     } elsif ($arg eq '-ob' or $arg eq '--outer_blank') {
# 	$fill  = undef;
# 	$outer = 1;
#     } elsif ($arg eq '-h1') {
# 	$header = 1;
#     } elsif ($arg eq '-h2') {
# 	$header = 2;
#     } elsif ($arg eq '-e' or $arg eq '--empty') {
# 	$empty = shift @ARGV;
#     } elsif ($arg eq '-skip') {
# 	$skip_empty_lines = 1;
#     } elsif ($arg eq '-m') {
# 	$merge = 1;
#     } elsif ($arg eq '-num') {
#         $numeric = 1;
#     } elsif($arg eq '-neg') {
#         $negate = 1;
#     } elsif($arg eq '-t' or $arg eq '-d') {
#         $arg = shift @ARGV;
#         $delim_in1 = $delim_in2 = $arg;
#     } elsif($arg eq '-di1') {
#         $delim_in1 = shift @ARGV;
#     } elsif($arg eq '-di2') {
#         $delim_in2 = shift @ARGV;
#     } elsif(($arg eq '-di') or ($arg eq '-d')) {
#         $delim_in1 = shift @ARGV;
#         $delim_in2 = $delim_in1;
#     } elsif($arg eq '-do') {
#         $delim_out = shift @ARGV;
#     } elsif($arg eq '-s1') {
# 	# Suppress printing of values from table 1 (key will be printed however).
#         $suppress1 = 1;
#     } elsif($arg eq '-s2') {
# 	# Suppress printing of values from table 2 (key will be printed however).
#         $suppress2 = 1;
#     } elsif($arg eq '-sk') {
#         $suppressk = 1;
#     } elsif($arg eq '-r' or $arg eq '--reverse' or $arg eq '-rev') {
#         $reverse = 1;
#     } elsif($arg eq '-u') {
# 	$uppercase = 1;
#     } elsif(length($file1) < 1) {
#         $file1 = $arg;
#     } elsif(length($file2)<1) {
#         $file2 = $arg;
#     } else {
# 	print STDERR "join.pl: UNRECOGNIZED COMMAND LINE ARGUMENT: $arg\n\n";
#     }
# }

# if (scalar(@fill_lines) > 0) {
#     if (!defined($fill) or (length($fill) == 0)) {
# 	$fill = join($delim_out, @fill_lines);
#     } else {
# 	$fill .= ($delim_out . join($delim_out, @fill_lines));
#     }
# }

# if($reverse) {
#     ($suppress1, $suppress2) = ($suppress2, $suppress1); # swap them!
# }


# if ( (length($file1) < 1) or (length($file2) < 1) ) {
#   print STDERR "join.pl: ERROR: Two input files must be specified on the command line!\n\n";
#   print STDERR <DATA>;
#   exit(1);
# }

# for(my $i=0; $i<=$#key1; $i++) { $key1[$i]--; } # what is going on. We subtract one from all keys I guess.
# for(my $i=0; $i<=$#key2; $i++){ $key2[$i]--; } # what the heck. We subtract one from all keys I guess.

# @key1 = sort {$a <=> $b} @key1; # sorts @key1 alphabetically
# @key2 = sort {$a <=> $b} @key2; # sorts @key2 alphabetically

# # print STDERR "Key 1: [", join(',', @key1), "]\n",
# #         "Key 2: [", join(',', @key2), "]\n",
# #         "Input delimiter 1: [$delim_in1]\n",
# #         "Input delimiter 2: [$delim_in2]\n",
# #         "Output delimiter 1: [$delim_out1]\n",
# #         "Output delimiter 2: [$delim_out2]\n",
# #         "\n",
# #         ;

# # Read in the key-printable pairs from the second file:
# my ($loops) = 0;
# my ($passify) = 10000; # <-- this means "print out a dot to STDERR every *this* many iterations"

# ## Note: We TRANSPARENTLY handle gzipped / bzip2 files in "openFile" , in "libfile.pl"
# my $fileRef2 = openFile($file2) or die("join.pl: Could not open file <$file2>.");

# if($verbose) {
#     print STDERR "join.pl: Reading relations from ", ($file2 eq '-') ? "standard input" : "file <$file2>";
# }

# #my $commentChar = '#';

# my $numDuplicateKeysRead = 0; # how many duplicate keys were read from $file2?
# my $maxDuplicateKeyWarnings = 10; # how many times to nag the user about multiple keys?

# my $header_data='';
# my $line=0;
# while(<$fileRef2>) {
#     if((not($skip_empty_lines) or /\S/)) { # and not(/^\s*$commentChar/)) {
# 	$line++;
# 	if($line==1 and $header==2) { $header_data = $_; }
# 	my @tmp = split($delim_in2);
# 	chomp($tmp[$#tmp]);
	
# # print STDERR "\n2: tmp: [", join('|',@tmp), "]\n";
# # print STDERR "2: key cols: [", join('|',@key2), "]\n";
# 	$key='';
# 	for(my $i=$#key2; $i >= 0; $i--) {
# 	    my $key_part = splice(@tmp,$key2[$i], 1);
# 	    # $key .= length($key)>0 ? ($delim_out . $key_part) : $key_part;
	    
# 	    if (length($key) == 0) {
# 		$key = "$key_part";
# 	    } else {
# 		$key = "${key_part}${delim_out}${key}";
# 	    }
	    
# # print STDERR "1: key: [$key]\n";
# # print STDERR "Before splice: [", join('|',@tmp), "] [$key]\n";
# # print STDERR "After splice: [", join('|',@tmp), "] [$key]\n";
# 	    # $key .= splice(@tmp, $i-1, 1) . $delim_out;
# 	}
# 	# Get rid of the last delimiter we added:
# 	if ($numeric) { $key = int($key); }
# 	if ($uppercase) { $key =~ tr/a-z/A-Z/; }
# # print STDERR "2: tmp: [", join('|',@tmp), "]\n";
# # print STDERR "2: key: [$key]\n";
# 	# $tmp = join($delim_out, @tmp);
# 	my $tmpMini = &myJoin($delim_out, \@tmp, $empty_placeholder);
	
# 	if (defined($key) && (length($key) > 0) && defined($exists{$key})) {
# 	    $numDuplicateKeysRead++;
# 	    if ($verbose && ($numDuplicateKeysRead <= $maxDuplicateKeyWarnings)) {
# 		# only print this error if we are in VERBOSE mode.
# 		print STDERR "\njoin.pl: WARNING: Multiple lines with the key \"$key\" were found in <$file2>. We will only use the LAST row with this key\n";
# 	    }
# 	}
# 	$values{$key} = $tmpMini;
# 	$exists{$key} = 1;

# 	my $tuple_size = scalar(@tmp);
	
# 	if(not(defined($max_tuple_size)) or $tuple_size > $max_tuple_size) {
# 	    $max_tuple_size = $tuple_size;
# 	}
	
# 	$loops++;
# 	if($verbose and (($loops % $passify) == ($passify-1))) {
# 	    print STDERR '.';
# 	}
#     }
# }
# if($verbose) { print STDERR " done (" . scalar(keys(%values)) . ").\n"; }
# close($fileRef2);

# if($outer and not(defined($fill))) {
#     $fill = &replicate($max_tuple_size, $empty_placeholder, $delim);
# }

# # Read in the key-printable pairs from the first file and print out the joined key:
# my $fileRef1 = openFile($file1) or die("join.pl: Could not open file <$file1>.");

# if($verbose) { print STDERR "join.pl: Joining on file <$file1>\n"; }
# $loops = 0;
# my $found;
# while(<$fileRef1>) {
#     if((not($skip_empty_lines) or /\S/)) { # and not(/^\s*${commentChar}/)) {
# 	my @tmp = split($delim_in1);
# 	chomp($tmp[$#tmp]);
# 	$key='';
# 	# print STDERR "1: tmp: [", join('|',@tmp), "]\n";
# 	for(my $i = $#key1; $i >= 0; $i--) {
# 	    my $key_part = splice(@tmp,$key1[$i], 1);
# 	    $key = (length($key)>0) ? ($key_part . $delim_out . $key) : $key_part;
# 	    # print STDERR "1: key: [$key]\n";
# 	}
# 	# Get rid of the last delimiter we added:
# 	if($numeric)   { $key = int($key); }
# 	if($uppercase) { $key = uc($key); }  # uppercase-ify it
	
# 	$value1 = &myJoin($delim_out, \@tmp, $empty_placeholder);
# 	$value2 = $values{$key};
# 	$found = $exists{$key};
	
# # print STDERR "1: key: [$key]\n";
# # print STDERR "1: value1: [$value1]\n";
# # print STDERR "1: value2: [$value2]\n";
# # print STDERR "1: found: [$found]\n";
# 	if((not($negate) and $found) or ($negate and not($found))) {
# 	    $hit = 1;
# 	} elsif($outer) {
# 	    $value2 = $fill;
# 	    $hit = 1;
# 	} else {
# 	    $hit = 0;
# 	}

# 	if($merge and $hit) {
# 	    $value1 = $value2;
# 	    $value2 = $empty;
# 	} elsif($merge and not($hit)) {
# 	    $hit = 1;
# 	    $value2 = $empty;
# 	}

# 	if ($hit) {
# 	    if ($reverse) {
# 		($value1, $value2) = ($value2, $value1) # Swap the two values if we're supposed to print the second value before the first one
# 	    }
# 	    $printable  = $suppressk ? '' : $key;
# 	    $printable1 = ($suppress1 or length($value1)<1) ? '' : $value1;
# 	    $printable2 = ($suppress2 or length($value2)<1) ? '' : $value2;

# 	    if (length($printable1) > 0) {
# 		$printable = (length($printable) > 0)
# 		    ? ($printable . $delim . $printable1)
# 		    : $printable1;
# 	    }
# 	    if (length($printable2) > 0) {
# 		$printable = (length($printable) > 0)
# 		    ? ($printable . $delim . $printable2)
# 		    : $printable2;
# 	    }
# 	    $printable =~ s/$empty_placeholder/$empty/g;
# 	    print $printable, "\n";
# 	}
# 	$loops++;
# 	if ($verbose and (($loops % $passify)==($passify-1)) ) {
# 	    print STDERR '.';
# 	}
#     }
# }
# close($fileRef1);
# if ($verbose) { print STDERR "join.pl finished. (Joined <$file1> with <$file2>)\n"; }
# if ($verbose && ($numDuplicateKeysRead > 0)) {
#     print STDERR "join.pl: WARNING: Read $numDuplicateKeysRead lines with identical keys in <$file2>.\n";
#     print STDERR "         Only the *last* line in a file with a redundant key is actually used.\n"
# }


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

-m: Merge - if key exists in file 2, use value from file 2 else
    use the value from file 1.

-q: Quiet mode: turn verbosity off (default is verbose)
    Quiet mode also removes the "more than one line with the key..." warnings.

-1 COL: Include column COL from FILE1 as part of the key (default is 1).
         Multiple columns may be specified in which case keys are constructed
         by concaternating each column in NUMERICAL order (NOT the order specified,
         where the delimiter used is equal to the output delimiter, see -do flag).
 CAVEAT: Multiple join fields do NOT respect the order on the command line
         -1 1,2,3,4,5  is exactly the same as -1 4,3,5,1,2. Therefore you cannot
         join files with out-of-order keys in this manner!

-2 COL: Include column COL from FILE2 as part of the key (default is 1).
         See -1 option for discussion of multiple columns.

-o FILLER: Do a left outer join.  If a key in FILE1 is not in FILE2, then the
          tuple from FILE1 is printed along with the text in FILLER in place of
          a tuple from FILE2 (by default these tuples are not reported in the
          result).  See -of option also to supply FILLER from a file.
          (See below for an example of usage.)

-ob: Same as -o but fill with blanks.

-of FILE: Same as -o, but use the text in FILE as the text for FILLER.

-e EMPTY: Set the empty string to EMPTY (default is blank).  If both keys
           exist in FILE1 and FILE2 but one tuple is blank, then the empty
           character EMPTY will be printed.

-num: Treat the keys as numeric quantities (default is off: treat keys as
        strings).  If this is turned on, each key will be forced into an
        integer quantity.

-neg: Negative output -- print keys that are in FILE1 but not in FILE2.
        These keys are the same ones that would be left out of the join,
        or those that would have a FILL tuple in a left outer join (see -o
        option).

-di1 DELIM: Set the input delimiter for FILE1 to DELIM (default is tab).

-di2 DELIM: Set the input delimiter for FILE2 to DELIM (default is tab).

-di DELIM: Set the input delimiters for both FILE1 and FILE2 to DELIM (default
        is tab).  Equivalent to using both the -di1 and -di2 options.

-do DELIM: Set the output file delimiter to DELIM. (default is tab)

-s1: Suppress printing of tuples from FILE1.  The key is printed, followed by
        the tuple found in FILE2.

-s2: Suppress printing of tuples from FILE2.

-sk: Suppress printing of keys.

-r: Reverse
     Instead of printing: <key> <tuple_from_FILE1> <tuple_from_FILE2>,
                    print <key> <tuple_from_FILE2> <tuple_from_FILE1>

-u: Uppercase.  (Case-insensitive join)
     Converts non-numeric keys to uppercase, resulting in a case-insensitive join.
     Keys from both FILE1 and FILE2 will be uppercased before attempting the join.

-skip: Skip empty lines (default processes them).


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

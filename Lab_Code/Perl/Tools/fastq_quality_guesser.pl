#!/usr/bin/perl -w

# This script detects what quality format a FASTQ file is
# It can operate on .fastq, .fq, .gzip, .bz2 files

# Example usage: fastq_quality_guesser.pl myFile.fq myfile2.fq.bz2 another_file.fastq.gz

use POSIX      qw(ceil floor);
use List::Util qw(max min);

#use File::Basename;
use Getopt::Long;
use strict;  use warnings;  use diagnostics;

# See the long description of quality scores at: http://en.wikipedia.org/wiki/FASTQ_format

my $S_SANGER_MIN       = 33; #  33 == '!'
my $S_SANGER_ALT_MIN   = 35; #  35 == '#' # Somtimes the lowest value is a '#' instead of a '!'... not sure why.
my $S_SANGER_MAX       = 73; #  73 == 'I'

my $X_SOLEXA_MIN       = 59;  #  59 == ';'
my $X_SOLEXA_MAX       = 105; #104; # 104 == 'h' (possibly sometimes 'i' ?)

my $I_ILLUMINA_1_3_MIN = 64;  #  64 == '@' ('at' sign)
my $I_ILLUMINA_1_3_MAX = 105; #104; # 104 == 'h' (possibly sometimes 'i' ?)

my $J_ILLUMINA_1_5_MIN = 66;  #  66 == 'B' (letter B)
my $J_ILLUMINA_1_5_MAX = 105; #104; # 104 == 'h' (possibly sometimes 'i' ?)

my $L_ILLUMINA_1_8_MIN     = 33; #  33 == '!'
my $L_ILLUMINA_1_8_ALT_MIN = 35;  # 35 == '#'  # Somtimes the lowest value is a '#' instead of a '!'... not sure why.
my $L_ILLUMINA_1_8_MAX     = 74; # 74 == 'J'

my $M_IONTORRENT_MIN = 36; # '$'
my $M_IONTORRENT_MAX = 71; # 'G'

my $NUM_QUAL_LINES_BEFORE_MAKING_A_GUESS = 100; # Default: check this many lines before making a guess
my $OUT_DELIM                            = "\t";
my $SHOULD_REPORT_ALTERNATIVES           = 0;

my $numInvalidFiles = 0; # To start with, there are no invalid files... yet.

sub main();

sub quitWithUsageError($) { print($_[0] . "\n"); printUsageAndQuit(); print($_[0] . "\n"); }
sub printUsageAndQuit() { printUsage(); exit(1); }

sub printUsage() {
    print STDOUT <DATA>;
    exit(0);
}


sub numPossibleQualScoresRemaining($) { # takes a hash as input
    # Returns the sum of the items in the hash
    # In this case, this sum here tells us how many possible quality score types have NOT YET
    # been ruled out.
    # The hash must contain ONLY 1s and 0s
    my ($theHashRef) = @_;
    my $sum = 0;
    while (my ($key, $value) = each(%{$theHashRef})) {
	#print "$key => $value\n";
	($value == 0 or $value == 1) or die "What, this hash can only contain 1 and 0, anything else is a programming error!";
	$sum += $value;
    }
    #print "Returning value of $sum...\n\n";
    return ($sum);
}

sub isProbably($$$$) {
    my ($lo, $hi, $expectedLo, $expectedHi) = @_;
    if (defined($lo) && defined($hi) && $lo==$expectedLo &&    ( $hi==$expectedHi || $hi==($expectedHi-1) )    ) {
	return 1; # note that we allow the HIGH value to not exactly match the expected Wikipedia high value. This is because in practice we see the quality score 105 ('i'), but Wikipedia says the max score is actually 104. Not sure what is up with that. The low score matches what Wikipedia says, however.
    } else {
	return 0; # the ACTUALLY SEEN values were not in this exact range
    }
}

sub getStringWithFilename($$$$$) {
    my ($theFilename, $theHashRef, $theKey, $theLowChar, $theHighChar) = @_;
    return ("FILE=" . $theFilename . ${OUT_DELIM} . (exists($$theHashRef{$theKey}) ? $$theHashRef{$theKey} : $$theHashRef{"INVALID"}) . ${OUT_DELIM} . "LOWEST=${theLowChar}" . ${OUT_DELIM} . "HIGHEST=${theHighChar}");
}




# ==1==
sub main() { # Main program
    $Getopt::Long::passthrough = 1; # ignore arguments we don't recognize in GetOptions, and put them in @ARGV. Important for reading in files!

    my $should_warn_the_user_about_alternatives = 0; # don't tell them about the "alternative" flag unless there's a reason to

    GetOptions("help|?|man" => sub { printUsageAndQuit(); }
	       , "delim|d=s" => \$OUT_DELIM
	       , "a|report-alternatives!" => \$SHOULD_REPORT_ALTERNATIVES
	       , "n|lines=i" => \$NUM_QUAL_LINES_BEFORE_MAKING_A_GUESS
	) or printUsageAndQuit();

    my $numUnprocessedArgs = scalar(@ARGV);
    ($numUnprocessedArgs >= 1) or quitWithUsageError("Error in arguments! You must send AT LEAST ONE fastq (or bzip / gzip compressed!) file to this program!\n");
    ($NUM_QUAL_LINES_BEFORE_MAKING_A_GUESS >= 1) or quitWithUsageError("You can't set --lines (-n) to a value less than 1.\n");

    my %outputAnnotation = (       "S sanger"         => join($OUT_DELIM, ("QUALITY=SANGER",       "SCALE=Phred+33",  "OFFSET=33",      "TOPHAT_PARAM=(nothing, phred33 is default)",   "RANGE=[0,40]",            "LETTERS=[! or #,I]") )
				   , "X solexa"       => join($OUT_DELIM, ("QUALITY=SOLEXA",       "SCALE=Solexa+64", "OFFSET=64",      "TOPHAT_PARAM=--solexa-quals",                  "RANGE=[-5,40] or [0,45]", "LETTERS=[;,h]") )
				   , "I illumina 1.3" => join($OUT_DELIM, ("QUALITY=ILLUMINA_1.3", "SCALE=Phred+64",  "OFFSET=64",      "TOPHAT_PARAM=--solexa1.3-quals",               "RANGE=[0,40]",            "LETTERS=[\@,h]") )
				   , "J illumina 1.5" => join($OUT_DELIM, ("QUALITY=ILLUMINA_1.5", "SCALE=Phred+64",  "OFFSET=64",      "TOPHAT_PARAM=--solexa1.3-quals",               "RANGE=[3,40]",            "LETTERS=[B,h]") )
				   , "L illumina 1.8" => join($OUT_DELIM, ("QUALITY=ILLUMINA_1.8", "SCALE=Phred+33",  "OFFSET=33",      "TOPHAT_PARAM=(nothing, phred33 is default)",   "RANGE=[0,41]",            "LETTERS=[! or #,J]") )
				   , "M iontorrent"   => join($OUT_DELIM, ("QUALITY=IONTORRENT"  , "SCALE=Phred???",  "OFFSET=36?",     "TOPHAT_PARAM=(Probably phred)",                "RANGE=[36,71]",           "LETTERS=[\$,G]") )
				   , "INVALID"        => join($OUT_DELIM, ("QUALITY=INVALID",       "SCALE=INVALID",  "OFFSET=INVALID", "TOPHAT_PARAM=INVALID",                         "RANGE=INVALID",           "LETTERS=INVALID") )
	); # <-- to initialize the hash properly, these MUST be parens and not braces!!!

    foreach my $fname (@ARGV) { # these were arguments that were not understood by GetOptions
	# =============== Below: RESET certain variables that are file-specific ===================
	my $lowestValue = undef; # stores the INTEGER VALUE of the lowest score
	my $highestValue = undef;
	my $lowestChar = undef; # stores the CHARACTER value of the lowest score
	my $highestChar = undef;

	my %canBe = (); # hash
	foreach my $key (keys(%outputAnnotation)) {
	    $canBe{$key} = 1; # <-- Reset the "can be hash" to 1 every time (meaning that a new file could be this option)
	}
	# ================ [FINISHED] resetting variables that are specific to each file ===========
	chomp($fname);
	
	#print STDERR "Handling filename: $fname\n";

	my $READER = "cat"; # default way to read a file is just to "cat" it
	if ($fname =~ /.gz$/  or $fname =~ /.gzip$/) {	    $READER = "gzip --decompress --stdout "; }
	if ($fname =~ /.bz2$/ or $fname =~ /.bzip2$/) {     $READER = "bzip2 --decompress --stdout "; }
	
	my $lineNum = 0; # <-- counts ALL lines in this one file, not just quality score lines
	my $numQualScoreLinesRead = 0; # <-- counts ONLY quality score lines in this one file
	my $numPossibleScoreTypes = scalar(keys %canBe);
	open(FILE, "$READER $fname | ") or die("cat failed\n");
	while ( defined(my $line = <FILE>) and ($numPossibleScoreTypes != 1) and ($numQualScoreLinesRead < $NUM_QUAL_LINES_BEFORE_MAKING_A_GUESS)) {
	    $lineNum++; # <-- Increment this at the BEGINNING; the very first line of the file will see this as a value of 1.
	    if (0 != $lineNum % 4) { # Every 4th line (i.e. remainder 3 when modulo-divided by 4) is a quality score line! 
		next; # <-- This is NOT a qual score line! There are four lines in each FASTQ record, but only the last item in each group is the quality score.
	    }

	    if ($line =~ m/^\s+$/) {
		next; # Skip any blank lines! Sometimes seen at the end of the file.
	    }
	    
	    # Ok, this is a quality score line!
	    chomp($line);
	    #print STDERR "Debugging: " . $line . "\n";
	    my @qualitiyChars = split(//, $line);
	    my $qLength = scalar(@qualitiyChars);
	    for (my $i = 0; $i < $qLength; $i++) {
		my $theChar  = $qualitiyChars[$i];
		my $theValue = ord($theChar); # numeric value of $theChar
		#print $theChar . "\t" . $theValue; print "\n"; # debugging
		if (!defined($lowestValue) or $theValue < $lowestValue) {
		    $lowestValue = $theValue;
		    $lowestChar = $theChar;
		}
		if (!defined($highestValue) or $theValue > $highestValue) {
		    $highestValue = $theValue;
		    $highestChar = $theChar;
		}
		#print "Hi: $highestValue   Low: $lowestValue\n"; print "Hi: $highestChar  Low: $lowestChar\n"; # debugging

		if ($lowestValue < $S_SANGER_MIN       || $highestValue > $S_SANGER_MAX) {       $canBe{"S sanger"} = 0; } # can't be Sanger anymore...
		if ($lowestValue < $X_SOLEXA_MIN       || $highestValue > $X_SOLEXA_MAX) {       $canBe{"X solexa"} = 0; } # can't be Solexa anymore...
		if ($lowestValue < $I_ILLUMINA_1_3_MAX || $highestValue > $I_ILLUMINA_1_3_MAX) { $canBe{"I illumina 1.3"} = 0; } # can't be Illumina 1.3+ anymore...
		if ($lowestValue < $J_ILLUMINA_1_5_MIN || $highestValue > $J_ILLUMINA_1_5_MAX) { $canBe{"J illumina 1.5"} = 0; } # can't be Illumina 1.5+ anymore...
		if ($lowestValue < $L_ILLUMINA_1_8_MIN || $highestValue > $L_ILLUMINA_1_8_MAX) { $canBe{"L illumina 1.8"} = 0; } # can't be Illumina 1.8+ anymore...
		if ($lowestValue < $M_IONTORRENT_MIN   || $highestValue > $M_IONTORRENT_MAX) {   $canBe{"M iontorrent"} = 0; } # can't be Illumina 1.8+ anymore...
	    }
	    #print $line . "...\n";
	    $numPossibleScoreTypes = numPossibleQualScoresRemaining(\%canBe);
	    $numQualScoreLinesRead++; # Every FOUR lines is a qual score: this should be one! Note that this is separate from "$lineNum", which counts ALL lines (even non-quality score lines)
	}
	close(FILE); # Done reading the file, now we need to report back...
	if ($numPossibleScoreTypes < 1) {
	    print STDERR "\n------------------------------------------\n";
	    print STDERR "WARNING: fastq_quality_guesser.pl is skipping the file <$fname>...\n";
	    print STDERR "<$fname>: ALL score types were disqualified by the time we got to line $lineNum in the file.\n";
	    print STDERR "<$fname>: Perhaps the file is malformed, or there is a problem in the code for 'fastq_which_quality_score.pl' OR possibly this is a new score type we've never seen before!\n";
	    print STDERR "<$fname>: The lowest  value was <$lowestValue>, which is character '$lowestChar'.\n";
	    print STDERR "<$fname>: The highest value was <$highestValue>, which is character '$highestChar'.\n";
	    print STDERR "------------------------------------------\n";
	    print STDOUT getStringWithFilename($fname, \%outputAnnotation, undef, $lowestChar, $highestChar) . "\n"; # Print an "INVALID" entry
	    $numInvalidFiles++;
	} elsif ($numPossibleScoreTypes == 1) {
	    # Great, there is exactly one possibility!
	    while (my ($key, $value) = each(%canBe)) {
		($value == 1) && print STDOUT getStringWithFilename($fname, \%outputAnnotation, $key, $lowestChar, $highestChar) . "\n";
	    }
	} else {
	    # Ok, but although there are actually MULTIPLE possibilities for this score, it is most likely to be a specific one...
	    # ===================== HANDLE THE GUESSING, when there are multiple options ===================
	    my %guessHash = (); #%canBe;
	    foreach my $key (keys(%canBe)) {
		$guessHash{$key} = 0; # Set all the hash VALUES to zero (the keys stay the same)
	    }
	    if (isProbably($lowestValue, $highestValue, $S_SANGER_MIN, $S_SANGER_MAX) # <-- Note that we will allow both '!' and '#' as options for the min score.
		or isProbably($lowestValue, $highestValue, $S_SANGER_ALT_MIN, $S_SANGER_MAX)) {      $guessHash{"S sanger"} = 1; }
	    
	    if (isProbably($lowestValue, $highestValue, $X_SOLEXA_MIN, $X_SOLEXA_MAX)) {             $guessHash{"X solexa"} = 1; }
	    if (isProbably($lowestValue, $highestValue, $I_ILLUMINA_1_3_MIN, $I_ILLUMINA_1_3_MAX)) { $guessHash{"I illumina 1.3"} = 1; }
	    if (isProbably($lowestValue, $highestValue, $J_ILLUMINA_1_5_MIN, $J_ILLUMINA_1_5_MAX)) { $guessHash{"J illumina 1.5"} = 1; }

	    if (isProbably($lowestValue, $highestValue, $L_ILLUMINA_1_8_MIN, $L_ILLUMINA_1_8_MAX)  # <-- Note that we will allow both '!' and '#' as options for the min score.
		or isProbably($lowestValue, $highestValue, $L_ILLUMINA_1_8_ALT_MIN, $L_ILLUMINA_1_8_MAX)) { $guessHash{"L illumina 1.8"} = 1; }

	    if (isProbably($lowestValue, $highestValue, $M_IONTORRENT_MIN, $M_IONTORRENT_MAX))     { $guessHash{"M iontorrent"} = 1; }

	    #print ("Can be...\n");
	    if ($guessHash{"S sanger"} and $guessHash{"L illumina 1.8"}) { # The sanger and illumina 1.8 only differ by one quality base--sanger tops out at 'I', while illumina tops out at J.
		#print ("Did this trigger\n");
		#print ("Highest was $highestChar .. $highestValue \n");
		if ($highestChar eq 'I') { $guessHash{"L illumina 1.8"} = 0; } # If it can be both, but we only saw (at most) an I, then we will call it a SANGER sequence.
		if ($highestChar eq 'J') { $guessHash{"S sanger"}       = 0; } # If it can be both, but we only saw (at most) an I, then we will call it an ILLUMINA 1.8 sequence.
		# If it can be both, but we only saw (at most) a J, then we will call it an ILLUMINA sequence.
	    }
	    
	    my $numBestGuesses = numPossibleQualScoresRemaining(\%guessHash);
	    if ($numBestGuesses != 1) {
		my $noGuessMessage = "NO_GUESS: Bad news, we could NOT make a best quality score guess for <$fname>! We had a total of $numPossibleScoreTypes possible scores to guess from, and $numBestGuesses guesses. The lowest character was '$lowestChar' ($lowestValue), the highest was '$highestChar' ($highestValue)";
		print STDOUT $noGuessMessage . "\n";
		print STDERR "STDERR: $noGuessMessage\n";
		if (!$SHOULD_REPORT_ALTERNATIVES) { $should_warn_the_user_about_alternatives = 1; }
	    }
	    while (my ($key, $value) = each(%guessHash)) {
		# Best guess!
		($value == 1) && print STDOUT getStringWithFilename($fname, \%outputAnnotation, $key, $lowestChar, $highestChar) . "\n";
	    }
	    # ===================== DONE HANDLING THE GUESSING, when there are multiple options ===================
	    if ($SHOULD_REPORT_ALTERNATIVES) {
		while (my ($key, $value) = each(%canBe)) {
		    next if ($numBestGuesses > 0 and $key eq "INVALID"); # If there are ANY valid guesses at all, then DO NOT also print "INVALID" as a possible guess.
		    ($value == 1) && print STDOUT ("ALTERNATIVE_OPTION_FOR_" . getStringWithFilename($fname, \%outputAnnotation, $key, $lowestChar, $highestChar)) . "\n";
		}
	    }
	}
    } # <-- end of handling each individual file

    if ($should_warn_the_user_about_alternatives) { print STDERR "STDERR: *** Note that you can specify -a (--report-alternatives) to report all the alternative possible scores that this file could theoretically be.\n"; }
} # end main()

main();

if ($numInvalidFiles > 0) {
    print STDERR "Exit code 1: WARNING: There were this many files with apparently-invalid quality scores: $numInvalidFiles.\n";
    exit(1);
} else {
    exit(0);
}
# ====

__DATA__

fastq_guess_quality_score.pl [OPTIONS]  [FASTQ_FILENAMES]
    * by Alex Williams, Dec. 2013
    * Reads .bz2 and .gz files automatically. Guesses compression from file extensions.
    * Does not support reading from STDIN; fastq data must be read from the filesystem.
    * See examples below.

OPTIONS:

  -d DELIMITER or --delim=DELIMITER   (Default: tab)
     Sets the OUTPUT delimiter to DELIMITER.

  -a or --report-alternatives (Default: FALSE)
     Intead of just reporting the BEST guess (for fastq files that could theoretically have multiple options),
     reports ALL possibilities. Normally this is unnecessary. Keep in mind that since some of the score types
     are strict subsets of others, it is impossible to 100% know for certain which score type was chosen in
     some cases.

  -n or --lines (Default: 100)
     Check THIS many quality score lines (note: not TOTAL lines in the file, just counting the quality score lines!)
     before making a guess, when the scores are ambiguous up to that point.
 
EXAMPLES:

fastq_quality_guesser.pl  -d "|"  file1.fq  file2.fq.gz  file3.fq.bz2
  Reports the quality of the THREE files specified, one per line. Pipe-delimited ("|").

fastq_quality_guesser.pl  --report-alternatives  file1.fq  file2.fq.bz2
  Reports quality scores PLUS alternate score possibilities in case there are any.
  Note that --report-alternatives is usually NOT A USEFUL OPTION!

fastq_quality_guesser.pl --help
  Displays this help

Score descriptions from Wikipedia:
   S - Sanger        Phred+33,  raw reads typically (0, 40)
       ASCII values are between 33 ("!") and 73 ("I")
   
   X - Solexa        Solexa+64, raw reads typically (-5, 40)
       ASCII values are between 59 (";") and 104 ("h")
   
   I - Illumina 1.3+ Phred+64,  raw reads typically (0, 40)
       ASCII values are between 64 ("@") and 104 ("h") or possibly 105 ("i")
       (Note from Alex: I have seen 'i' as a valid value in files detected
        as Illumina 1.3+ Phred+64 data--even though Wikipedia says that 'h' is
        the highest value, apparently sometimes it is 'i'.)
   
   J - Illumina 1.5+ Phred+64,  raw reads typically (3, 40)
       with 0=unused, 1=unused, 2=Read Segment Quality Control Indicator
       ASCII values are between 66 ("B") and 104 ("h")
   
   L - Illumina 1.8+ Phred+33,  raw reads typically (0, 41)
       ASCII values are between 35 ("#") and 74 ("J")

   M - IonTorrent **unknown** quality score (probably Phred+36 ?)
       ASCII values are between 36 ("$") and 71 ("G")

CAVEATS:

None so far.

KNOWN ISSUES:

  Does not support piping in via STDIN (i.e., you CANNOT do 'cat file.fq | fastq_quality_guesser.pl -')

--------------

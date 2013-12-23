#!/usr/bin/perl -w

# This script detects what quality format a FASTQ file is
# It can operate on .fastq, .fq, .gzip, .bz2 files

# Example usage: fastq_what_quality.pl myFile.fq myfile2.fq.bz2 another_file.fastq.gz

use POSIX      qw(ceil floor);
use List::Util qw(max min);

#use File::Basename;
use Getopt::Long;
use strict;  use warnings;  use diagnostics;


my $S_SANGER_MIN       = 33; #  33 == '!'
my $S_SANGER_MAX       = 73; #  73 == 'I'

my $X_SOLEXA_MIN       = 59;  #  59 == ';'
my $X_SOLEXA_MAX       = 105; #104; # 104 == 'h' (possibly sometimes 'i' ?)

my $I_ILLUMINA_1_3_MIN = 64;  #  64 == '@'
my $I_ILLUMINA_1_3_MAX = 105; #104; # 104 == 'h' (possibly sometimes 'i' ?)

my $J_ILLUMINA_1_5_MIN = 66;  #  66 == '@'
my $J_ILLUMINA_1_5_MAX = 105; #104; # 104 == 'h' (possibly sometimes 'i' ?)

my $L_ILLUMINA_1_8_MIN = 2;  #  2 == '#'
my $L_ILLUMINA_1_8_MAX = 74; # 74 == 'J'

my $NUM_QUALS_TO_CHECK_BEFORE_MAKING_A_GUESS = 10; # Check this many lines before making a guess

sub main();

sub quitWithUsageError($) { print($_[0] . "\n"); printUsageAndQuit(); print($_[0] . "\n"); }
sub printUsageAndQuit() { printUsage(); exit(1); }

sub printUsage() {
    print STDOUT <DATA>;
    exit(0);
}


# ==1==
sub main() { # Main program
    my ($delim) = "\t";
    my ($decimalPlaces) = 4; # How many decimal places to print, by default
    $Getopt::Long::passthrough = 1; # ignore arguments we don't recognize in GetOptions, and put them in @ARGV

    GetOptions("help|?|man" => sub { printUsageAndQuit(); }
	) or printUsageAndQuit();

    if (1 == 0) {
	quitWithUsageError("1 == 0? Something is wrong!");
    }

    my $numUnprocessedArgs = scalar(@ARGV);
    if ($numUnprocessedArgs < 1) {
	quitWithUsageError("Error in arguments! You must send AT LEAST ONE fastq (or bzip / gzip compressed!) file to this program!\n");
    }

    my $lineNum = 0;
    my $numQualScoreLinesRead = 0;

    my %canBe = (       "S sanger"         => 1   # It's a hash that keeps track of how many possible quality scores haven't been disqualified yet.
			, "X solexa"       => 1   # <-- note that the values are all expected to be set to ONE initially!
			, "I illumina 1.3" => 1  # <-- if something CAN BE this quality type, it's set to 1
			, "J illumina 1.5" => 1  # <-- if this quality score type has been ruled out, it becomes ZERO
			, "L illumina 1.8" => 1
	); # <-- to initialize the hash properly, these MUST be parens and not braces!!!

    my %outputAnnotation = (       "S sanger"         => "SANGER\tPhred+33\tRange: (0,40)\t('!','I')"
				   , "X solexa"       => "SOLEXA\tSolexa+64\tRange: (-5,40) or (0,45)\t(';','h')"
				   , "I illumina 1.3" => "ILLUMINA_1.3\tPhred+364\tRange: (0,40)\t('@','h')"
				   , "J illumina 1.5" => "ILLUMINA_1.5\tPhred+64\tRange: (3,40)\t('B','h')"
				   , "L illumina 1.8" => "ILLUMINA_1.8\tPhred+33\tRange: (0,41)\t('!','I')"
	); # <-- to initialize the hash properly, these MUST be parens and not braces!!!


    #sort(keys(%canBe))
    #for (my $i = 0; $i < scalar(keys %canBe); $i++);
    

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
	return ($sum);
    }


    foreach my $fname (@ARGV) { # these were arguments that were not understood by GetOptions
	chomp($fname);
	
	print STDERR "Handling filename: $fname\n";

	my $READER = "cat"; # default way to read a file is just to "cat" it
	if ($fname =~ /.gz$/  or $fname =~ /.gzip$/) {	    $READER = "gzip --decompress --stdout "; }
	if ($fname =~ /.bz2$/ or $fname =~ /.bzip2$/) {     $READER = "bzip2 --decompress --stdout "; }
	
	my $lowestValue = undef; # stores the INTEGER VALUE of the lowest score
	my $highestValue = undef;

	my $lowestChar = undef; # stores the CHARACTER value of the lowest score
	my $highestChar = undef;

	my $qualStillNeedsChecking = 1;

	open(FILE, "$READER $fname | ") or die("cat failed\n");

	my $numPossibleScoreTypes = scalar(keys %canBe);
	while ( (my $line = <FILE>) || $qualStillNeedsChecking) {
	    if (3 == $lineNum % 4) {
		chomp($line);
		print $line . "\n";
		my @qualitiyChars = split(//, $line);
		my $qLength = scalar(@qualitiyChars);
		for (my $i = 0; $i < $qLength; $i++) {
		    my $theChar  = $qualitiyChars[$i];
		    my $theValue = ord($theChar); # numeric value of $theChar
		    
		    print $theChar . "\t" . $theValue;
		    #print "\n";

		    if (!defined($lowestValue) or $theValue < $lowestValue) {
			$lowestValue = $theValue;
			$lowestChar = $theChar;
		    }

		    if (!defined($highestValue) or $theValue > $highestValue) {
			$highestValue = $theValue;
			$highestChar = $theChar;
		    }

		    #print "Hi: $highestValue   Low: $lowestValue\n";
		    #print "Hi: $highestChar  Low: $lowestChar\n";

		    if ($lowestValue < $S_SANGER_MIN || $highestValue > $S_SANGER_MAX) { $canBe{"S sanger"} = 0; } # can't be Sanger anymore...
		    if ($lowestValue < $X_SOLEXA_MIN || $highestValue > $X_SOLEXA_MAX) { $canBe{"X solexa"} = 0; } # can't be Solexa anymore...
		    if ($lowestValue < $I_ILLUMINA_1_3_MAX || $highestValue > $I_ILLUMINA_1_3_MAX) { $canBe{"I illumina 1.3"} = 0; } # can't be Illumina 1.3+ anymore...
		    if ($lowestValue < $J_ILLUMINA_1_5_MIN || $highestValue > $J_ILLUMINA_1_5_MAX) { $canBe{"J illumina 1.5"} = 0; } # can't be Illumina 1.5+ anymore...
		    if ($lowestValue < $L_ILLUMINA_1_8_MIN || $highestValue > $L_ILLUMINA_1_8_MAX) { $canBe{"L illumina 1.8"} = 0; } # can't be Illumina 1.8+ anymore...
		}

		print "Num types remaining is: " . (numPossibleQualScoresRemaining(\%canBe)) . "\n";

		$numPossibleScoreTypes = numPossibleQualScoresRemaining(\%canBe);
		if (1 == $numPossibleScoreTypes) {
		    print STDOUT "Great, we got the answer! Printing it to STDOUT.";
		    #exit(0);
		    #print
		}
		if (0 == $numPossibleScoreTypes) {
		    print STDOUT "Uh oh, ALL score types were disqualified!";
		    exit(1);
		}
		
		$numQualScoreLinesRead++;
		# Every FOUR lines is a qual score: this should be one!
	    } else {
		# This is NOT a qual score line! There are four lines in each FASTQ record, but only the last item in each group is the quality score.
	    }


	    sub isProbably($$$$) {
		my ($lo, $hi, $expectedLo, $expectedHi) = @_;
		if ($lo == $expectedLo and ($hi == $expectedHi or ($hi == ($expectedHi - 1)))) {
		    return 1; # note that we allow the HIGH value to not exactly match the expected Wikipedia high value. This is because in practice we see the quality score 105 ('i'), but Wikipedia says the max score is actually 104. Not sure what is up with that. The low score matches what Wikipedia says, however.
		} else {
		    return 0; # the ACTUALLY SEEN values were not in this exact range
		}
	    }

	    sub getStringWithFilename($$$) {
		my ($theFilename, $theHashRef, $theKey) = @_;
		return ($theFilename . "\t" . $$theHashRef{$theKey});
	    }

	    if ($numPossibleScoreTypes == 1 or $numQualScoreLinesRead >= $NUM_QUALS_TO_CHECK_BEFORE_MAKING_A_GUESS) {
		# Ok, so we haven't TOTALLY narrowed it down to exactly one option only.
		# But we probably have a good guess.

		print STDOUT "It is still theoretically possible that the following are the score types:\n";

		while (my ($key, $value) = each(%canBe)) {
		    ($value == 1) && print STDOUT getStringWithFilename($fname, \%outputAnnotation, $key) . "\n";
		}
		
		if ($numPossibleScoreTypes > 1) {
		    # Ok, but although there are actually MULTIPLE possibilities for this score, it is most likely to be a specific one:
		    
		    my %guessHash = %canBe;
		    foreach my $key (keys(%guessHash)) {
			$guessHash{$key} = 0; # Set all the hash VALUES to zero (the keys stay the same)
		    }
		    # if (isProbably($lowestValue, $highestValue, $S_SANGER_MIN, $S_SANGER_MAX))             $guessHash{"S sanger"} = 1;
		    # if (isProbably($lowestValue, $highestValue, $X_SOLEXA_MIN, $X_SOLEXA_MAX))             $guessHash{"X solexa"} = 1;
		    # if (isProbably($lowestValue, $highestValue, $I_ILLUMINA_1_3_MIN, $I_ILLUMINA_1_3_MAX)) $guessHash{"I illumina 1.3"} = 1; 
		    # if (isProbably($lowestValue, $highestValue, $J_ILLUMINA_1_5_MIN, $J_ILLUMINA_1_5_MAX)) $guessHash{"J illumina 1.5"} = 1;
		    # if (isProbably($lowestValue, $highestValue, $L_ILLUMINA_1_8_MIN, $L_ILLUMINA_1_8_MAX)) $guessHash{"L illumina 1.8"} = 1;
		    
		    # print STDOUT "But our best guess is:\n";
		    # while (my ($key, $value) = each(%guessHash)) {
		    # 	($value == 1) && print STDOUT getStringWithFilename($fname, \%outputAnnotation, $key) . "\n";
		    # }
		}

	    }

	    $lineNum++; # <-- done checking this line
	}
	close(FILE);
    } # <-- end of handling each individual file
} # end main()

main();

exit(0);
# ====

__DATA__

fastq_guess_quality_score.pl  [FASTQ_FILENAMES]
    * .bz2 and .gz files can also be read automatically
    * Does not support reading from pipes ('|'); fastq data must be read from the filesystem.

by Alex Williams, Dec. 2013

# 
# Description from Wikipedia:
# S - Sanger        Phred+33,  raw reads typically (0, 40)
# ASCII values are between 33 ("!") and 73 ("I")

# X - Solexa        Solexa+64, raw reads typically (-5, 40)
# ASCII values are between 59 (";") and 104 ("h")

# I - Illumina 1.3+ Phred+64,  raw reads typically (0, 40)
# ASCII values are between 64 ("@") and 104 ("h") or possibly 105 ("i") (I have seen 'i' as a valid value in files detected as Illumina 1.3+ Phred+64 data, despite what Wikipedia says)

# J - Illumina 1.5+ Phred+64,  raw reads typically (3, 40)
#     with 0=unused, 1=unused, 2=Read Segment Quality Control Indicator
# ASCII values are between 66 ("B") and 104 ("h")
#

# L - Illumina 1.8+ Phred+33,  raw reads typically (0, 41)
# ASCII values are between 35 ("#") and 74 ("J")


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

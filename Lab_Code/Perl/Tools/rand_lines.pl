#!/usr/bin/perl

## By Alex Williams, 2010.

## Named the same as a previous script from UCSC, but shares *none* of the same code.

# Randomly choose a certain number of lines from a file.
# Supports multi-line records, as long as they are always the SAME line.

use strict;
use Getopt::Long;
Getopt::Long::Configure(qw(pass_through)); # pass_through: "Don't flag options we don't process as errors"

my $numToPrint = undef;
my $proportion = undef;
my $numLinesPerRecord = 1;
my $allowDupes = 0;
my $numHeaders    = 0;
my $blanks        = 1;
my $verbose       = 1;
my $randSeed      = undef;
my $suppressPrintingOfHeaders = 0;

GetOptions("help|man|?" => sub { print STDERR <DATA>; exit(0); }
	   , "n=s"            => \$numToPrint
	   , "p|proportion=s" => \$proportion
	   , "r|num_lines_per_record=i" => \$numLinesPerRecord
	   , "q|quiet"     => sub { $verbose = 0; }
	   , "wr"          => sub { $allowDupes = 1; } ## with replacement
	   , "wor"         => sub { $allowDupes = 0; } ## without replacement
	   , "h|headers=i" => \$numHeaders
	   , "srand|seed|randseed=i" => \$randSeed
	   , "nb|noblanks|no_blank|no_blanks"          => sub { $blanks = 0; }
	   , "sh|suppress_header|suppress_headers"     => sub { $suppressPrintingOfHeaders = 1; }
	   );

(($numLinesPerRecord >= 1) or die "The number of lines per record cannot be less than 1.");

if (defined($randSeed)) {
    $verbose && print STDERR ("Using the value " . $randSeed . " to seed the random number generator. This way we will get the same \'random\' lines each time...\n");
    srand($randSeed); ## <-- if you set srand, you get THE SAME "random" values each time
}

my @lines = <>;

my @headerLines = ();

if ($numHeaders > 0) {
    @headerLines = @lines[0..($numHeaders-1)]; ## save the headers
    @lines = @lines[$numHeaders..(scalar(@lines)-1)]; ## remove the headers from "@lines"
}

if (!$blanks) {
    ## Remove blank lines
    $verbose && print STDERR "Removing all blank lines from input, due to the -nb option...\n";
    @lines = grep/\S/, @lines; ## Find all non-blank lines only
}

my $numRecords = scalar(@lines) / $numLinesPerRecord;

if (scalar(@lines) % $numLinesPerRecord != 0) {
    die "Uh oh! The number of lines in the file *NOT INCLUDING THE HEADER* was " . scalar(@lines) . ", which is not evenly divisible by the number of rows per record, which you specified with the -r option as being " . $numLinesPerRecord . ".\nMaybe you should specify no blank lines (with --noblanks)--there could be a trailing blank line.\nIs each \"record\" in the file really a group of " . $numLinesPerRecord . "?\nOr maybe you specified the wrong number of header lines, or forgot that this file had a header. (Default number of header lines is 0. Specify a single header line with --headers=1.\n\n";
}




((!defined($numToPrint) || ($numToPrint >= 0)) or die "The number of lines to obtain cannot be less than 0.");

((!defined($proportion) || ((0 <= $proportion) && ($proportion <= 1.0))) or die "Hey, the proportion of lines to print must be between 0.0 and 1.0, inclusive.");

if (!defined($numToPrint) && defined($proportion)) {
    $numToPrint = ($proportion * $numRecords);    
}

if (!defined($numToPrint)) { $numToPrint = $numRecords; }

if ($numToPrint > $numRecords && !$allowDupes) {
    print STDERR "Warning: you specified more lines to select randomly than there are lines in the input to rand_lines.pl, but you also said you didn't want lines to be repeated (the -wor option, which is the default).\nSo to solve this problem, we set the actual number of lines to print to the number of lines in the file. You could also specify -wr (with replacement) to allow repeated lines.\n\n";
    $numToPrint = $numRecords;
}


$verbose && print STDERR ("Number of lines to print: " . $numToPrint . "\n");

$verbose && print STDERR ("Number of lines in file: " . $numRecords . "\n");

$verbose && $allowDupes && print STDERR ("Allowing duplicate lines to be picked. Specify -wor to pick lines randomly *without* replacement.\n");

my @x = (); ## array of indices of which line to pick

if (!$allowDupes) {

    #use List::Util qw(first max maxstr min minstr reduce shuffle sum);
    ## We could have, instead, just used "shuffle" to perl-style shuffle the array! Note that we would have to
    ## figure out how srand interacts with shuffle. (To make sure we get the same results each time.)"

    $#x = $numRecords;
    for (my $i = 0; $i < $numRecords; $i++) { $x[$i] = $i; }

    for (my $i = 0; $i < $numRecords; $i++) {
	my $switchWith = int($numRecords*rand());
	my $temp = $x[$i];
	$x[$i] = $x[$switchWith];
	$x[$switchWith] = $temp;
    }
} else {
    ## Allow duplicates!
    $#x = $numToPrint;
    for (my $i = 0; $i < $numToPrint; $i++) {
	my $thisRandNum = int($numRecords*rand());
	$x[$i] = $thisRandNum;
    }
    #for (my $i = 0
}



## DEBUGGING: print the indices that we chose randomly
#print STDERR "Here are the indices that we chose (numbering starts at 0 rather than 1):\n";
#print STDERR (join("\n",@x),"\n");

if (!$suppressPrintingOfHeaders) {
    print STDOUT @headerLines;
}

## x is the INDICES in order
for (my $i = 0; $i < $numToPrint; $i++) {
    my $theBaseIndexForThisRecord = ($x[$i] * $numLinesPerRecord); ## Note that this is counting in LINES, not records.
    for (my $rec = 0; $rec < $numLinesPerRecord; $rec++) {
	## When a record is multi-lined, we want to go through each of the lines and print it.
	print STDOUT ($lines[$theBaseIndexForThisRecord + $rec]);
    }
}




__DATA__

syntax: rand_lines.pl FILENAME

Randomly shuffles a file around and chooses the first N random lines to print out.

Examples:
        rand_lines.pl -n 100 MY_BIG_FILE.txt > a_hundred_random_lines_from_the_file.txt

        rand_lines.pl --randseed=123 -n 10 SOME_FILE.txt > always_same_10_random_lines.txt

        rand_lines.pl -
Supports more-than-one-line-per-record files. So it can deal with, for example, FASTA files.

Options:
  -n INTEGER  : Print *NUMBER* randomly-chosen records. Default: print all lines in the file.

  -p NUMBER or --proportion=NUMBER   : Print a *NUMBER* fraction of randomly-picked records.
             1.0 prints all records, 0.5 prints half of them, etc.

  -r INTEGER or --num_lines_per_record=INTEGER : Specifies that records are more than one line.
             Good for FASTA files. This means that every INTEGER lines will be taken as being a single
             mega-line.

  -wor       : Without replacement---choose lines without replacement (default)

  -wr        : With replacement---allow lines to be chosen more than once.

  -h INTEGER or --headers=INTEGER : How many header lines are in the file? These will be
             printed to the output, unless --suppress_headers is specified.

  -sh or --suppress_headers: Do not print any header lines to the output.

  -nb or --noblanks:  Eliminate all blank lines.

  --randseed=INTEGER : Use INTEGER as a "seed" for the random number generator. This means
             the same "random" lines will be selected each time. Uses perl's srand(...) function
             to do the random seed. Note: might NOT be portable between machines!

  -q         : Quiet. (Not verbose.) Suppress debugging info. Default is verbose.

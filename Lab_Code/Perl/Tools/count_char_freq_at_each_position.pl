#!/usr/bin/perl

use strict;
use warnings;
use diagnostics;

use List::Util qw(sum);

## Counts character frequency at each position in a file
## No delimiters are accepted.
## Example:
## File is:
### ABCDA
### ZZZZZZ

## Output would be:
#         	1	2	3	4	5	6	TOTAL_LETTERS
#A        	1	0	0	0	1	0	2
#B        	0	1	0	0	0	0	1
#C        	0	0	1	0	0	0	1
#D        	0	0	0	1	0	0	1
#Z        	1	1	1	1	1	1	6
#POSITION_TOTAL	2	2	2	2	2	1	11

print STDERR "### count_char_freq_at_each_position.pl is now running...\n";
print STDERR "###  Example usage:  cat myfile | count_char_freq_at_each_position.pl > OUTPUT.txt\n";
print STDERR "###  Note that this program accepts a single file argument OR a STDIN pipe.\n";
print STDERR "###  About to count up the frequencies of letters at each position in the input.\n";
print STDERR "###  Output will be a tab-delimited table of positions (columns) by the letters found at each position (rows).\n";
print STDERR "###  If a literal tab is found, it will be reported as the letter \"<tab>\".\n";
print STDERR "###  Now attempting to read from STDIN...\n";
print STDERR "###  (If no more text appears, please double-check that you specified an input file!)\n";

my $lineLength = undef;

my $lineNum = 1;

my @lineStats = ();
my @lineTotals = ();

my $numLinesWereDifferentLengths = 0;
while (my $line = <>) {
    if ($lineNum == 1) { print STDERR ("### [OK]: Successfully reading from STDIN.\n"); }
    chomp($line);
    if (!defined($lineLength)) {
	$lineLength = length($line);
    } else {
	if ((length($line) != $lineLength) && (length($line) > 0)) {
	    $numLinesWereDifferentLengths++;

	    if ($numLinesWereDifferentLengths < 100) {
		print STDERR "### NOTE: Line $lineNum was length ", length($line), " but we expected it to be length ", $lineLength, ". Line $lineNum is: $line\n";
	    } elsif ($numLinesWereDifferentLengths == 100) {
		print STDERR "Not outputting any more lines-were-different-length errors. Check the input file to see if this is a problem.\n...\n...Continuing with the input...\n";
	    }
	}
    }
    
    my @la = split('', $line);
    for (my $i = 0; $i < scalar(@la); $i++) {
	
	my $c = $la[$i];
	if ($c eq "\t") {
	    $c = "<tab>";
	}

	if (!exists($lineStats[$i])) {
	    %{$lineStats[$i]} = ();
	    $lineTotals[$i] = 0;
	}

	if (!exists($lineStats[$i]{$c})) {
	    ${$lineStats[$i]}{$c} = 1;
	} else {
	    ${$lineStats[$i]}{$c}++;
	}

	$lineTotals[$i]++;
    }
    $lineNum++;
    
    if ($lineNum % 100000 == 0) { print STDERR ("### count_char_freq_at_each_position.pl has processed $lineNum input lines so far...\n"); }
}

## ============================ PRINT HEADER ===================
## ============================ and initialize %allKeys ========
my %allKeys = ();
for (my $i = 0; $i < scalar(@lineStats); $i++) {
    print STDOUT ("\t" . ($i+1)); ## <-- header line!
    foreach my $key (sort (keys(%{$lineStats[$i]}))) {
	$allKeys{$key} = 1;
    }
}
print STDOUT ("\t" . "COLUMN_SUM" . "\n"); ## <-- end of header line
## ============================ DONE PRINTING HEADER ===========


if (scalar(@lineTotals) != scalar(@lineStats)) { die "Problem in programming: count_char_freq_at_each_position.pl expected lineTotals and lineStats to be the same length!\n"; }

foreach my $key (sort (keys(%allKeys))) {
    print STDOUT ($key);
    my $sumForThisLetter = 0;
    for (my $i = 0; $i < scalar(@lineStats); $i++) {
	my $numberOfThisLetter = (exists($lineStats[$i]{$key})) ? $lineStats[$i]{$key} : 0;
	$sumForThisLetter += $numberOfThisLetter;
	print STDOUT ("\t" . $numberOfThisLetter);
    }
    print STDOUT "\t" . $sumForThisLetter . "\n"; ## print the total, then the end of this line
}

## ============================ PRINT POSITION TOTALS ===========
print STDOUT "RSUM";
for (my $i = 0; $i < scalar(@lineTotals); $i++) {
    print STDOUT ("\t" . $lineTotals[$i]);
}
print STDOUT "\t" . sum(@lineTotals) . "\n"; ## end of totals line

if ($numLinesWereDifferentLengths) {
    print STDERR "### Warning: $numLinesWereDifferentLengths lines in the input that we just counted were different length from the first line in the input! Check the STDERR output to see the list of lines.\n";
}

print STDERR "### \n";
print STDERR "### Processed input and wrote the results to STDOUT. Columns indicate character positions in the file, rows give the characters found at each position.\n";
print STDERR "### Finished: Read a total of $lineNum lines.\n";

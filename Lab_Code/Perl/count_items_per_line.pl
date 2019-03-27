#!/usr/bin/env perl

# By Alex Williams, July 2006
# alexgw@soe.ucsc.edu

# Input: pipe in a file
# Output: output to STDOUT the number of delimiters per line

use strict;
use Getopt::Long;

my $delimiter = "\t";
my $printTitle = 0;
my $printUniqOnly = 0;
my $stop_if_different_n_fields = 0;

sub printUsage {
    print STDOUT <DATA>;
    exit(0);
}

GetOptions("help|?|man"  => sub { printUsage(); }
	   , "delim|d:s" => \$delimiter
	   , "title|t"   => sub { $printTitle = 1; }
	   , "uniq|u"    => sub { $printUniqOnly = 1; }
	   , "strict|x"  => sub { $stop_if_different_n_fields = 1; }
	   ) or printUsage();

my $prev_n                = -1; # must be an invalid number to start
my $i                     = 1;  # lines in a file are numbered from 1
while (<>) {
    my $line = $_;
    my @lineArr = split(/$delimiter/, $line);
    my $numDelimiters = (scalar(@lineArr) - 1); #() = $line =~ /${delimiter}/g; # count how many delimiters are on this line

    if (!$printUniqOnly or ($prev_n != $numDelimiters)) {
	    # If we are printing everything, OR the number of fields on this line is different, then print it.
	    if ($printTitle) {
		    my $title = $lineArr[0];
		    chomp($title); # just in case it was the only item in its row, so it has a newline
		    print STDOUT $title . "\t"; # Note that $lineArr[1] is the FIRST element ([0] contains the entire string)
	    }
	    print STDOUT ($numDelimiters) . "\n";

    	    if ($stop_if_different_n_fields and defined($prev_n) and $prev_n >= 0 and $prev_n != $numDelimiters) {
		    print STDERR "STRICT MODE: Failure on line $i: expected ${prev_n} delimiters, found ${numDelimiters}\n";
		    print STDERR "Offending line number $i is shown below:\n${line}\n\n";
		    exit(1);
	    }
    }
    $i++;
    $prev_n = $numDelimiters;
}

exit(0);


__DATA__

count_items_per_line:
 Counts the number of items per line, AFTER the initial item, which is assumed to be a header for
 the row (although you can change that with -addamt=1) . Counts based on the delimiter on that line.

Does NOT transparently work on gzipped files; you should use 'gzcat' or 'zcat' for those.

IMPORTANT NOTE:
  You may want to use "row_stats.pl -count" instead of this program!

Usage:
  count_items_per_line.pl [OPTIONS]  <FILE|STDIN>


Options:
 -d=DELIMITER or  --delim=DELIMITER   Default: tab
    Sets the column delimiter. Can be a regular expression, but make sure to format it properly.
    Uses the perl "split" command, so it can take any regular expression that split can understand.
    Normal usage example: -d=',' would count the number of items on a column-delimited file

-x or --strict: Exits with error code 1 if a differing number of items on a row was found, after the first row.
  Useful for checking a file that you may want to load into a database.

-u or --uniq: For each run of numbers, only print a single representative. E.g., could print 40, 30, 40, 30, 40
    if 100 lines had 40 elements, 1 had 30, 100 more had 40, etc.

-t  or --title  (Default: do not print titles, just counts)
    Print the item title at the beginning of each line, then a tab, then the count of remaining items.

Q: If you are using this program, have you considered using row_stats.pl -h 0 -count instead? Note that the output is slightly different due to headers!

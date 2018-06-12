#!/usr/bin/env perl
#@COMMENT@ transpose.pl transposes a table of data (flips the rows/columns). Example:  cat myfile.txt | transpose.pl > flipped.txt. Frequency-of-use rating: 10/10.
# Original version by Josh Stuart, 2000. Modified by Alex Williams, 2017.
use strict; use warnings; use Getopt::Long;
my ($verbose, $delim) = (0, "\t");    #echo -e "\n\n\na\nb\n1,2,3,4,5,6\na,b,c,d,e,f\ng,h,i\n,,z"
$Getopt::Long::passthrough = 1; # ignore arguments we don't recognize in GetOptions, and put them in @ARGV
GetOptions("help|?|man"   => sub { print STDOUT <DATA>;exit(0); }
	   , "delim|d=s"  => \$delim
	   , "v|verbose!" => \$verbose
	  ) or do { print STDOUT <DATA>;exit(1); };
#sub max($$) { $_[0] > $_[1] ? $_[0] : $_[1]; } # return max
my @final; # final output, this will be a two-dimensional table
my ($ncols, $nrows) = (-1, 0); # save this outside the loop below!!!
foreach my $line (<>) { # <-- read STDIN or a passed-in file
	chomp($line);
	my @row = split($delim, $line, -1); # -1 is important! otherwise we don't split on things where a column is entirely empty
	@{$final[$nrows]} = @row; # there's probably some way to do this with 'push'
	if (scalar(@row) > $ncols) { $ncols = scalar(@row); }
	$nrows++;
}
#my $nrows = scalar(@final);
$verbose and print STDERR "Transposing $nrows ROWS by $ncols COLUMNS...\n";
for (my $cc = 0; $cc < $ncols; $cc++) {
	for (my $rr = 0; $rr < $nrows; $rr++) {
		($rr > 0) and print STDOUT "$delim";
		print STDOUT defined($final[$rr][$cc]) ? $final[$rr][$cc] : '';
	}
	print STDOUT "\n";
}
exit(0);
__DATA__
transpose.pl [-d DELIM] [FILE | < FILE]

Transposes a table -- flips the rows and columns, writes to STDOUT.
The input file/STDIN is assumed to have rows delimited by newlines
and columns delimited by tabs.
Possibly requires UNIX-format newlines (\n).

FILE: a file containing a tabular data plain text file.
   e.g.:         MyHeader      Col2     Col3
                 Row1_yes      a        b
                 Row2_maybe    c        d

Example usage:     cat myfile.tab.txt | transpose.pl > out.txt
		 transpose.pl -d ',' file.csv > transposed.csv

OPTIONS:
   -d DELIM        : Change column delimiter (default: tab).
   or --delim=DELIM: Example: --delim=',' for a CSV file.
   -v      : More verbose output to STDERR
================================================================

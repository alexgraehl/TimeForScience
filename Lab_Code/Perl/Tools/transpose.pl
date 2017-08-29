#!/usr/bin/perl

#@COMMENT@ transpose.pl transposes a table of data (flips the rows/columns). Example:  cat myfile.txt | transpose.pl > flipped.txt. Frequency-of-use rating: 10/10.
use strict; use warnings; use Getopt::Long;
my $verbose = 0;
my $delim   = "\t";
$Getopt::Long::passthrough = 1; # ignore arguments we don't recognize in GetOptions, and put them in @ARGV
GetOptions("help|?|man"   => sub { print STDOUT <DATA>;exit(0); }
	   , "delim|d=s"  => \$delim
	   , "v|verbose!" => \$verbose
	  ) or do { print STDOUT <DATA>;exit(1); };
my $ncols = -1;
my @x;
#sub max($$) { $_[0] > $_[1] ? $_[0] : $_[1]; } # return max
my $r     = 0; # save this outside the loop below!!!
foreach my $line (<>) { # <-- read from STDOUT or from a passed-in file
    my @row = split($delim, $line);
    chomp($row[$#row]);
    my $c = 0;
    for ($c = 0; scalar(@row) > 0; $c++) {
	    $x[$r][$c] = shift(@row);
    }
    if ($c > $ncols) { $ncols = $c; }
    $r++;
}
my $nrows = $r;
$verbose and print STDERR " done ($nrows by $ncols)!\n";
$verbose and print STDERR "Transposing to $ncols by $nrows...\n";
for(my $cc = 0; $cc < $ncols; $cc++) {
    for(my $rr = 0; $rr < $nrows; $rr++) {
	my $x = defined($x[$rr][$cc]) ? $x[$rr][$cc] : '';
	print "$x";
	if ($r < $nrows - 1) { print $delim; }
    }
    print "\n";
}
exit(0);
__DATA__
----------------------------------------------------------------
transpose.pl [-d DELIM] [FILE | < FILE]
----------------------------------------------------------------
Transposes a table -- flips the rows and columns.
The original file is assumed to have rows delimited by newlines
and columns delimited by tabs.
Possibly requires UNIX-format newlines (\n).

FILE: a file containing a tabular data plain text file.
   Example:      MyHeader      Col2     Col3
                 Row1_yes      a        b
                 Row2_maybe    c        d

Example usage:     cat myfile.tab.txt | transpose.pl > out.txt
		 transpose.pl -d ',' file.csv > transposed.csv

OPTIONS:
   -d DELIM        : Change column delimiter (default: tab).
   or --delim=DELIM: Example: --delim=',' for a CSV file.

   -v      : More verbose output to STDERR
================================================================

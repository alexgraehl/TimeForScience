#!/usr/bin/perl

##############################################################################
##############################################################################
##
## edges2matrix.pl
##
##############################################################################
##############################################################################
##
## Written by Josh Stuart in the lab of Stuart Kim, Stanford University.
##
##  Email address: jstuart@stanford.edu
##          Phone: (650) 725-7612
##
## Postal address: Department of Developmental Biology
##                 Beckman Center Room B314
##                 279 Campus Dr.
##                 Stanford, CA 94305
##
##       Web site: http://www.smi.stanford.edu/people/stuart
##
##############################################################################
##############################################################################
##
## Written: 00/00/02
## Updated: 00/00/02
##
##############################################################################
##############################################################################

use lib "$ENV{MYPERLDIR}/lib"; use lib "$ENV{TIME_FOR_SCIENCE_DIR}/Lab_Code/Perl/LabLibraries"; require "libfile.pl";
#use lib "$ENV{MYPERLDIR}/lib"; use lib "$ENV{TIME_FOR_SCIENCE_DIR}/Lab_Code/Perl/LabLibraries"; require "libset.pl";

use strict;
use warnings;

print STDERR "Note: sets.pl does something similar to what this script does, but with more options!\n";

# Two options: give this program TWO COLUMNS of data or THREE COLUMNS.

# If you feed this program two columns, and it gives you a matrix.
# For example, if you have a list of gene interactions of some kind,
# it prints out "1" when the interaction is in the list, and "0" otherwise.

# If you feed this program THREE columns, it gives you a matrix with the cell value being whatever was in the third column.
# For example, if you have a list of gene interactions of some kind,
# it prints out "1" when the interaction is in the list, and "0" otherwise.

# Give it TWO columns,:
# A  B
# A  C
# A  D
# B  C
# D  E

# And it will give you back:
#    A  B  D
# B  1  0  0
# C  1  1  0
# D  0  0  0
# E  0  0  1

# If you give it THREE columns instead, then it gives you the value in the third column!

# Flush output to STDOUT immediately.
#$| = 1;

my @flags   = (
                  [    '-q', 'scalar',     0,     1]
                , [   '-k1', 'scalar',     1, undef]
                , [   '-k2', 'scalar',     2, undef]
                , [    '-v',  'scalar',    3, undef]
                , [    '-b',  'scalar',    0,     1]
                , [    '-d', 'scalar',  "\t", undef]
                , [    '-s', 'scalar',     0,     1]
                , ['--missing', 'scalar', "", undef]
                , ['--file', 'scalar',   '-', undef]
              );

my %args = %{&parseArgs(\@ARGV, \@flags)};

if(exists($args{'--help'}))
{
   print STDOUT <DATA>;
   exit(0);
}

my $verbose   = not($args{'-q'});

my $missingVal = $args{'--missing'};

my $key_col1  = $args{'-k1'}; # Note: these are numbered from ONE and not zero!
my $key_col2  = $args{'-k2'}; # Note: these are numbered from ONE and not zero!
my $val_col   = $args{'-v'};  # Note: these are numbered from ONE and not zero!

my $should_output_binary_options = $args{'-b'};
my $delim     = $args{'-d'};
my $symmetric = $args{'-s'};
my $file      = $args{'--file'};

my %sets;

my $filep;
open($filep, $file) or die("Could not open file '$file' for reading");

my $numItemsSeenBefore = undef;
my $lineNum = 0;

my %seen;

my %keySet1 = ();
my %keySet2 = ();

while(my $line = <$filep>) {
    $lineNum++;
    next if $line =~ /^\s+$/; # Skip ENTIRELY blank lines
    chomp($line);
    my @lar = split($delim, $line);
    my $numItems = scalar(@lar);
    my $item1 = "";
    my $item2 = "";
    my $value = "";

    if (defined($numItemsSeenBefore) && ($numItemsSeenBefore != $numItems)) { print STDERR "[WARNING]: matrixFromEdges.pl: on line $lineNum of <$file>, we only detected a total of $numItems delimited items, whereas we had earlier seen $numItemsSeenBefore items! This may be a serious error.\n"; }

    if ($key_col1 > $numItems) { print STDERR "[WARNING: matrixFromEdges.pl: on line $lineNum of <$file>. Num items on line: $numItems, but we were looking for a first key at position $key_col1 (counting from 1, not zero, but which is still out of bounds)!\n";
    } else { $item1 = $lar[$key_col1 - 1]; }

    if ($key_col2 > $numItems) { print STDERR "[WARNING: matrixFromEdges.pl: on line $lineNum of <$file>. Num items on line: $numItems, but we were looking for a second key at index $key_col2! (counting from 1, not zero, but which is still out of bounds)!\n";
    } else { $item2 = $lar[$key_col2 - 1 ]; }

    if ($should_output_binary_options) {
	$value = 1;
    } else {
	if (defined($val_col) && $val_col <= $numItems) { $value = $lar[$val_col - 1]; }
    }

    if (!exists($seen{$item1})) { %{$seen{$item1}} = (); }
    print "Value was $value...\n";
    print "Value col was $val_col and num items was $numItems...\n";
    $seen{$item1}{$item2} = $value;
    $keySet1{$item1} = 1; $keySet2{$item2} = 1;
    
    if ($symmetric) {
	if (!exists($seen{$item2})) { %{$seen{$item2}} = (); }
	$seen{$item2}{$item1} = $value;
	$keySet2{$item1} = 1; $keySet1{$item2} = 1; # note: numbers should be backwards!
    }

    $numItemsSeenBefore = $numItems;
}
close($filep);

# &setsPrint(\%sets, \*STDOUT, 0, 1);

my $outdelim = "\t";

my @sortedKeys1 = sort(keys(%keySet1));
my @sortedKeys2 = sort(keys(%keySet2));

# ================ HEADER LINE =======================
print "\t"; # <-- Empty cell in the top left corner of the output
foreach my $k2 (@sortedKeys2) {
    print $k2 . $outdelim;
}
print "\n"; # Done with the header line
# ================ [DONE] WITH HEADER LINE ============

foreach my $k1 (@sortedKeys1) {
    print $k1 . $outdelim;
    foreach my $k2 (@sortedKeys2) {
	my $cell = (exists($seen{$k1}{$k2})) ? $seen{$k1}{$k2} : $missingVal;
	print $cell . $outdelim;
    }
    print "\n";
}
#setsPrintMatrix(\%seen, \*STDOUT, "\t"); # <-- in libset.pl

exit(0);


__DATA__
syntax: edges2matrix.pl [OPTIONS]

OPTIONS are:

-q: Quiet mode (default is verbose)

-k1 COL: Set the column of item 1 to COL (default is 1).

-k2 COL: Set the column of item 2 to COL (default is 2).

VALUE options (default is -v 3):
   1)  -v COL: Set the CELL VALUE column to COL (default is 3, IF there is a third column!).
   2) -b (binary presence/absense output--no values are printed!)
            Instead of -v COL, you can also specify that you want BINARY output only (1 or 0) instead.
            This shows whether a pair occurred, but doesn't print a value in the cell besides a 0 or 1.

-d DELIM: Set the field delimiter to DELIM (default is tab).

-s: Symmetric.  Assume edges are undirected so add symmetric relationships.

--missing STRING (Default: (blank))
        * Sets the value for "not found." Examples: --missing='NA' or --missing=NONE.

USAGE EXAMPLES:

matrixFromEdges.pl --missing NA  INPUTFILE.txt
  * Output column 3 values (if present) based on keys in columns 1 and 2. Output "NA" for non-present cases.

matrixFromEdges.pl --missing NONE  -k1 2  -k2 3  -v 4 INPUTFILE.txt
  * Output column 4 values based on keys in columns 2 and 3.

matrixFromEdges.pl --missing NONE  -k1 1  -k2 2  -b   INPUTFILE.txt
  * Output binary "1 / 0" values based on keys in columns 1 and 2.


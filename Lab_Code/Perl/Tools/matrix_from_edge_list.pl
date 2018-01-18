#!/usr/bin/perl

#@COMMENT@ matrix_from_edge_list.pl can turn a 2- or 3-column file into a matrix. The matrix will either be an adjacency matrix (2 column input) or will have the values of each edge (3 column input). Frequency-of-use rating: 5/10.

##############################################################################
##############################################################################
## matrix_from_edge_list.pl or edges2matrix.pl
##############################################################################
##############################################################################
## Written by Josh Stuart in the lab of Stuart Kim, Stanford University.
##############################################################################
##############################################################################
## Written: 00/00/02
## Updated: 00/00/02
## Updated by AGW: Feb 10 2015
##############################################################################
##############################################################################

use lib "$ENV{MYPERLDIR}/lib"; use lib "$ENV{TIME_FOR_SCIENCE_DIR}/Lab_Code/Perl/LabLibraries"; require "libfile.pl";
use strict;
use warnings;

# Two options: give this program TWO COLUMNS of data or THREE COLUMNS.

# If you feed this program two columns, and it gives you a matrix.
# For example, if you have a list of gene interactions of some kind,
# it prints out "1" when the interaction is in the list, and "0" otherwise.

# If you feed this program THREE columns, it gives you a matrix with the cell value being whatever was in the third column.
# For example, if you have a list of gene interactions of some kind,
# it prints out "1" when the interaction is in the list, and "0" otherwise.

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
if (exists($args{'--help'})) {
   print STDOUT <DATA>;
   exit(0);
}
my $verbose   = not($args{'-q'});
my $key_col1  = $args{'-k1'}; # Note: these are numbered from ONE and not zero!
my $key_col2  = $args{'-k2'}; # Note: these are numbered from ONE and not zero!
my $val_col   = $args{'-v'};  # Note: these are numbered from ONE and not zero!
my $should_output_binary_options = $args{'-b'};
my $missingVal = ($should_output_binary_options) ? 0 : $args{'--missing'};
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

    if (defined($numItemsSeenBefore) && ($numItemsSeenBefore != $numItems)) { print STDERR "[WARNING]: matrix_from_edge_list.pl: on line $lineNum of <$file>, we only detected a total of $numItems delimited items, whereas we had earlier seen $numItemsSeenBefore items! This may be a serious error.\n"; }

    if ($key_col1 > $numItems) { print STDERR "[WARNING]: matrix_from_edge_list.pl: on line $lineNum of <$file>. Num items on line: $numItems, but we were looking for a first key at position $key_col1 (counting from 1, not zero, but which is still out of bounds)!\n";
    } else { $item1 = $lar[$key_col1 - 1]; }

    if ($key_col2 > $numItems) { print STDERR "[WARNING]: matrix_from_edge_list.pl: on line $lineNum of <$file>. Num items on line: $numItems, but we were looking for a second key at index $key_col2! (counting from 1, not zero, but which is still out of bounds)!\n";
    } else { $item2 = $lar[$key_col2 - 1 ]; }

    if ($should_output_binary_options) {
	$value = 1;
    } else {
	if (defined($val_col) && $val_col <= $numItems) { $value = $lar[$val_col - 1]; }
    }

    if (!exists($seen{$item1})) { %{$seen{$item1}} = (); }
    #print STDERR "Value was $value...\n";
    #print STDERR "Value col was $val_col and num items was $numItems...\n";

    #if (exists($seen{$item1}{$item2})) { print STDERR "[WARNING]: matrix_from_edge_list.pl: on line $lineNum of <$file>, we saw a duplicate key/value pair ($item1 and $item2) that we had also seen earlier!\n"; }

    $seen{$item1}{$item2} = $value;
    $keySet1{$item1} = 1; $keySet2{$item2} = 1;
    
    if ($symmetric) {
	if (!exists($seen{$item2})) { %{$seen{$item2}} = (); }
	$seen{$item2}{$item1} = $value;
	$keySet2{$item1} = 1; $keySet1{$item2} = 1; # note: numbers should be backwards!
    }

    $numItemsSeenBefore = $numItems;
}

# &setsPrint(\%sets, \*STDOUT, 0, 1);

my $outdelim = "\t";

my @sortedKeys1 = sort(keys(%keySet1));
my @sortedKeys2 = sort(keys(%keySet2));

# ================ HEADER LINE =======================
print STDOUT "\t"; # <-- Empty cell in the top left corner of the output
foreach my $k2 (@sortedKeys2) {
    print STDOUT $k2 . $outdelim;
}
print STDOUT "\n"; # Done with the header line
# ================ [DONE] WITH HEADER LINE ============

foreach my $k1 (@sortedKeys1) {
    print STDOUT $k1 . $outdelim;
    foreach my $k2 (@sortedKeys2) {
	my $cell = (exists($seen{$k1}{$k2})) ? $seen{$k1}{$k2} : $missingVal;
	print STDOUT $cell . $outdelim;
    }
    print STDOUT "\n";
}

exit(0); # Done

__DATA__
syntax: matrix_from_edge_list.pl [OPTIONS]

Note: related in some ways to 'flatten.pl' and 'expand.pl'--those may be useful
  for getting a list of edges into the right format.

Runs in two modes: TWO COLUMN mode and THREE COLUMN mode:
# Give it TWO columns,:
# A  B
# A  C
# A  D
# B  C
# D  E

# And it will give you back:
#    A  B  D
# B  1      
# C  1  1   
# D         
# E        1     <-- '1's indicate the spaces that had a pairing in the input file

# If you give it THREE columns instead, then instead of getting a 1, you get the value from that third column!
# In the case of duplicates, you get the LAST value seen. This is not a well-defined behavior.

* Note: "sets.pl" is very similar to this script, and has some additional features too!

OPTIONS are:

-b: Print 1 or 0 (present/absent) rather than trying to print a third column.
 
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

--missing STRING (Default: (NA))
        * Sets the value for "not found." Examples: --missing 'NA' or --missing NONE.
        * Note: Do not put an equal sign between <missing> and the string!

USAGE EXAMPLES:

matrix_from_edge_list.pl --missing NA  INPUTFILE.txt
  * Output column 3 values (if present) based on keys in columns 1 and 2. Output "NA" for non-present cases.

matrix_from_edge_list.pl --missing NONE  -k1 2  -k2 3  -v 4 INPUTFILE.txt
  * Output column 4 values based on keys in columns 2 and 3.

matrix_from_edge_list.pl -b   --missing NONE  -k1 1  -k2 2   INPUTFILE.txt
  * Output binary "1 / 0" values based on keys in columns 1 and 2.


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
use lib "$ENV{MYPERLDIR}/lib"; use lib "$ENV{TIME_FOR_SCIENCE_DIR}/Lab_Code/Perl/LabLibraries"; require "libset.pl";

use strict;
use warnings;

print STDERR "Note: sets.pl does something similar to what this script does, but with more options!";



# If you feed this program two columns, and it gives you a matrix.
# For example, if you have a list of gene interactions of some kind,
# it prints out "1" when the interaction is in the list, and "0" otherwise.


# If you feed this program THREE columns, it gives you a matrix with the cell value being whatever was in the third column.
# For example, if you have a list of gene interactions of some kind,
# it prints out "1" when the interaction is in the list, and "0" otherwise.


# There are TWO MODES for this:

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
$| = 1;

my @flags   = (
                  [    '-q', 'scalar',     0,     1]
                , [   '-k1', 'scalar',     1, undef]
                , [   '-k2', 'scalar',     2, undef]
                , [   '-v',  'scalar',     3, undef]
                , [    '-d', 'scalar',  "\t", undef]
                , [    '-s', 'scalar',     0,     1]
                , ['--file', 'scalar',   '-', undef]
              );

my %args = %{&parseArgs(\@ARGV, \@flags)};

if(exists($args{'--help'}))
{
   print STDOUT <DATA>;
   exit(0);
}

my $verbose   = not($args{'-q'});
my $key_col1  = $args{'-k1'};
my $key_col2  = $args{'-k2'};
my $val_col   = $args{'-v'};
my $delim     = $args{'-d'};
my $symmetric = $args{'-s'};
my $file      = $args{'--file'};

$key_col1--;
$key_col2--;

my %sets;

my $filep;
open($filep, $file) or die("Could not open file '$file' for reading");
while(<$filep>) {
   my @tuple = split($delim, $_);
   chomp($tuple[$#tuple]);
   my $item1 = $tuple[$key_col1];
   my $item2 = $tuple[$key_col2];

   if($item1 =~ /\S/ and $item2 =~ /\S/) {
      &addEntry(\%sets, $item1, $item2);
      if($symmetric) {    &addEntry(\%sets, $item2, $item1); }
   }
}
close($filep);

# &setsPrint(\%sets, \*STDOUT, 0, 1);
&setsPrintMatrix(\%sets, \*STDOUT);

exit(0);

sub addEntry {
   my ($sets, $from, $to) = @_;
   if(not(exists($$sets{$from}))) {
      my %set;
      $$sets{$from} = \%set;
   }
   my $set = $$sets{$from};
   $$set{$to} = 1;
}

__DATA__
syntax: edges2matrix.pl [OPTIONS]

OPTIONS are:

-q: Quiet mode (default is verbose)

-k1 COL: Set the column of item 1 to COL (default is 1).

-k1 COL: Set the column of item 2 to COL (default is 2).

-d DELIM: Set the field delimiter to DELIM (default is tab).

-s: Symmetric.  Assume edges are undirected so add symmetric relationships.




#!/usr/bin/perl

# Alex Williams

# Splits up a file from a matrix into its constituent columns:

# Inputs: first argument is the matrix file, second is the directory to put the output files into.

# Example:   ./split_up_matrix.pl   MY_MATRIX.txt   OUT_DIRECTORY

# TITLE  KEY1 KEY2 KEY3
# DATA1  0.4  0.2  0.14
# DATA2  4.1  4.5  3.2

# Into this:
# TITLE KEY1
# DATA1 0.4
# DATA2 4.1

# And this:
# TITLE KEY2
# DATA1 0.2
# DATA2 4.5

# Etc. for all of the columns

use strict;

my $matrixFile = $ARGV[0];
my $outDir = $ARGV[1];

my $firstLine = `head -n 1 $matrixFile`;
my @firstLineTokens = split(/\t/, $firstLine);
my $numColumns = scalar(@firstLineTokens);

# note, hard coded here to number of columns...
for (my $x=2; $x <= $numColumns; $x++) {
    system("mkdir -p $outDir");
    system("cut -f 1,$x $matrixFile > $outDir/tmp");
    my $title = `cut -f 2 $outDir/tmp | head -n 1`;
    chomp($title);
    system("mv $outDir/tmp $outDir/$title");
    print "$x ($title)\n";
}


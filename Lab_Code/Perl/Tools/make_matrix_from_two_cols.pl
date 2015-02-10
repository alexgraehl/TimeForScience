#!/usr/bin/perl -w

# note: sets.pl will probably do this, only better


die "Note: sets.pl probably does what this script does, only better.";

# You feed this program two columns, and it gives you a matrix.
# For example, if you have a list of gene interactions of some kind,
# it prints out "1" when the interaction is in the list, and "0" otherwise.

# Give it this:
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

# The FIRST column becomes the top of the matrix.
# The SECOND column becomes the index on the left of the matrix.

use strict;

sub printMatrix;

printMatrix();

sub printMatrix {
    my(%matrix); # The hash that will have the IDs as keys

    my(%first);
    my(%second);

    my($line);
    while ($line = <STDIN>) { # Read through the file...
	chomp($line); # Take off the newline, if there is one.
	if ($line =~ /^(\S+)\s(\S+)/) {
	    $matrix{$1 . "\t" . $2} = 1;
 	    $first{$1} = 1;
	    $second{$2} = 1;
	} else {
	    
	}
    } # end while
     
    my($isTopLine) = 1;
    foreach my $x (sort (keys(%second))) {
	if (!$isTopLine) {
	    print $x;
	}
	print "\t";
	foreach my $y (sort (keys(%first))) {
	    if ($isTopLine == 1) {
		print $y . "\t";
	    } else {
		if (defined($matrix{$x . "\t" . $y}) ||
		    defined($matrix{$y . "\t" . $x})) {
		    print '1';
		} else {
		    print '0';
 		}
		print "\t";
	    }
	}
	$isTopLine = 0;
	print "\n";
    }
    
}

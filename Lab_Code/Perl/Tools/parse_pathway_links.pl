#!/usr/bin/perl

# Alex Williams

# Takes a file like this:

# ORF_1 <TAB> pathway1 <TAB> pathway2 <TAB> pathway3
# ORF_2 <TAB> pathway2
# ORF_3 <TAB> pathway2 <TAB> pathway3

# And prints out this:

# ORF_1 <TAB> pathway1
# ORF_1 <TAB> pathway2
# ORF_1 <TAB> pathway3
# ORF_2 <TAB> pathway2
# ORF_3 <TAB> pathway2
# ORF_3 <TAB> pathway3


use strict;

while (my $line = <>) {
    chomp($line);
    my @a = split('\t', $line);

    my $ORF = shift(@a); #shift:  Shifts the first value of the array off and returns it, shortening the array by 1 and moving everything down. 
    # Save & remove the first item from the array (the first item is the ORF)

    if (defined($ORF) && length($ORF) > 0) {
	foreach my $pathway (@a) {
	    print $ORF . "\t" . $pathway . "\n";
	}
    }
}

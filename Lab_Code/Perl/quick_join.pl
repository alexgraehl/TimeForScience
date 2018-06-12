#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;

my $firstFileKey = 1;
my $secondFileKey = 1;


sub printUsage() {
    print STDOUT <DATA>;
    exit(0);
}


Getopt::Long::Configure("pass_through"); # don't give an error if there are extra unused arguments (those should be our two files)

GetOptions("help|?|man" => sub { printUsage(); }
	   , "1=i" => \$firstFileKey
	   , "2=i" => \$secondFileKey
	   ) or printUsage();



if (scalar(@ARGV) < 2) { die "Must include two files to join!"; }

my $primary   = $ARGV[0];
my $secondary = $ARGV[1];

my $suffix = ".sorted.quick.tmp";

my $sorted1 = "/tmp/${primary}${suffix}";
my $sorted2 = "/tmp/${secondary}${suffix}";

system("sort -k ${firstFileKey} ${primary}    > ${sorted1}");
system("sort -k ${secondFileKey} ${secondary} > ${sorted2}");

system("join -1 $firstFileKey -2 $secondFileKey   $sorted1 $sorted2");

#system("rm -f $sorted1 $sorted2");

exit(0);

__DATA__

trim_lines.pl -n LENGTH [< INPUT_FILE]

Trims lines to length LENGTH. -n=3 would trim every line to three characters.

Trims off the END of a line by default. Thus, the earlier characters on the line are preserved.

I wrote this script because the sed command (with parentheses and \1) for trimming lines is actually very slow.
This perl version is at least twice as fast.

#!/usr/bin/perl

use strict;  use warnings;  use diagnostics;
use Getopt::Long;
use Carp; # backtrace on errors. Has the "confess" function. Use this instead of "die" if you want useful information!   

$| = 1; # Always flush text output IMMEDIATELY to the console, don't wait to buffer terminal output! Setting this to zero can cause STDERR and STDOUT to be interleaved in weird ways.

sub main();
sub quitWithUsageError($) { print($_[0] . "\n"); printUsageAndQuit(); print($_[0] . "\n"); }
sub printUsageAndQuit() { printUsage(); exit(1); }

sub printUsage() {
    print STDOUT <DATA>;
    exit(0);
}

# ==1==
sub main() { # Main program
    my ($delim) = "\t";
    my ($decimalPlaces) = 4; # How many decimal places to print, by default
    $Getopt::Long::passthrough = 1; # ignore arguments we don't recognize in GetOptions, and put them in @ARGV

    GetOptions("help|?|man" => sub { printUsageAndQuit(); }
	       , "delim|d=s" => \$delim
	       , "dp=i" => \$decimalPlaces
	) or printUsageAndQuit();

    my $numUnprocessedArgs = scalar(@ARGV);
    if ($numUnprocessedArgs == 0) {
	quitWithUsageError("Error in arguments! You must send ONE OR MORE bam or sam filenames to this program.\n");
    }

    foreach my $in (@ARGV) { # these were arguments that were not understood by GetOptions

	(-f $in) or die "Could not find the file '$in'. Exiting now.";

	my $isSAM = 0;
	my $isBAM = 0;
	if ($in =~ /.sam$/i) {
	    $isSAM = 1;
	} elsif ($in =~ /.bam$/i) {
	    $isBAM = 1;
	} else {
	    die "The input file did not end in .sam or .bam, but it has to be a sam or bam file!";
	}
	
	my $out1 = $in;
	$out1 =~ s/.[sb]am$/.pair1.bam/i;
	my $out2 = $in;
	$out2 =~ s/.[sb]am$/.pair2.bam/i;

	print STDERR "Splitting the hopefully-paired-end file <$in> into the files <$out1> and <$out2>...\n";
	
	my $SAM_SPECIFIC_OPTIONS = ($isSAM) ? " -S " : " ";
	
	system("samtools view -b $SAM_SPECIFIC_OPTIONS -h -f 0x40 $in > $out1"); # always output a bam file!
	system("samtools view -b $SAM_SPECIFIC_OPTIONS -h -f 0x80 $in > $out2"); # always output a bam file!
    }
} # end main()


main();

exit(0);
# ====

__DATA__

sam_split_paired_end.pl  FILENAME.bam FILENAME2.sam

by Alex Williams, 2014

Takes a number of input PAIRED END sam or bam files, and splits them up by paired end-ness.

I.e., input is:

 sam_split_paired_end.pl  myfile.sam

Now the output will be two new files:

 myfile.pair1.bam  and  myfile.pair2.bam

Note that the output is ALWAYS bam format, even if the input is sam format.

OPTIONS:

 --help: Displays this help.
 
 No other options.

EXAMPLES:

sam_split_paired_end.pl  myfile.sam

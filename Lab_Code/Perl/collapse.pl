#!/usr/bin/env perl

##############################################################################
##############################################################################
## collapse.pl: related to 'expand.pl'. Makes a single convenient table out of a mega-list of pairs
# NOTE : not the same as 'expand.pl', which is the opposite of 'flatten.pl'

##############################################################################
##############################################################################
## Based on original version by Josh Stuart.
## Modified to reduce library dependencies 2017 by Alex Williams.
##############################################################################
##############################################################################
use strict; use warnings; use diagnostics;
use Getopt::Long;
use Carp; # backtrace on errors. Has the "confess" function. Use this instead of "die" if you want useful information!
sub printUsageAndQuit() { printUsage(); exit(1); }
sub printUsage() { print STDOUT <DATA>; }

sub main() {
	my $keycol    = 1; # key column
	my $delim_out = ",";	# default is ',', like "MyKey    10,2,42,91"
	my $delim_in  = "\t";   # Assumes the input file is tab-delimited
	my $n_header_lines   = 0; # Print and skip header lines, if there are any
	my $function  = undef;
	my $file      = '-';
	$Getopt::Long::passthrough = 1; # ignore arguments we don't recognize in Ge
	GetOptions("help|?|man" => sub { printUsageAndQuit(); } #, "q|quiet!" => \$quiet
		   , "k=i" => $keycol
		   , "i|di=s" => \$delim_in
		   , "o|do=s" => \$delim_out
		   , "f|function=s"  => \$function
		   , "file" => \$file
		   , "h|headers=i"  => \$n_header_lines
		   , "dot!" => sub { $delim_out = "\t"; }
		  ) or printUsageAndQuit();
	$keycol = $keycol-1;	 # subtract one since "1" really means "column 0" in perl-land
	my $filep;
	open($filep, $file) or die("[ERROR in collapse.pl]: Could not open file '$file' for reading: $!");
	while ($n_header_lines-- > 0) { $_ = <$filep>; print($_); } # print header lines (if any) without collapsing them
	my @order;
	my %data;
	while (my $line = <$filep>) {
		my @tuple = split($delim_in, $line);
		chomp($tuple[$#tuple]);
		my $key = splice(@tuple, $keycol, 1);
		if (not(exists($data{$key}))) {
			$data{$key} = \@tuple;
			push(@order, $key);
		} else {
			my $list = $data{$key};
			my $n    = scalar(@{$list}); # list length
			my $m    = scalar(@tuple);   # tuple length
			my $min  = ($n < $m) ? $n : $m;
			for (my $i = 0; $i < $min; $i++) {
				$$list[$i] .= $delim_out . $tuple[$i];
			}
			for (my $i = $m; $i < $n; $i++) {
				$$list[$i] .= $delim_out;
			}
		}
	}
	close($filep);

	# Print results to console
	foreach my $key (@order) {
		my $list = $data{$key};
		if (defined($function)) {
			# Only if 'function' is defined to we need to handle vec_eval / libstats.pl
			use lib "$ENV{MYPERLDIR}/lib"; use lib "$ENV{TIME_FOR_SCIENCE_DIR}/Lab_Code/Perl/LabLibraries"; require "libstats.pl";
			my @tuple;
			foreach my $field (@{$list}) {
				my @vector = split($delim_out, $field);
				my $val    = vec_eval(\@vector, $function); # <-- 'vec_eval' is a function from libstats.pl
				push(@tuple, defined($val) ? $val : '');
			}
			$list = \@tuple;
		}
		print $key, $delim_in, join($delim_in, @{$list}), "\n";
	}
}

main();

__DATA__
syntax: collapse.pl [OPTIONS]  < INPUT_FILE.pairs.txt

Turns a list of pairs into a more convenient table.
# NOTE : not the same as 'expand.pl', which is the opposite of 'flatten.pl'
  
For example, if we have this input in a file:
A  1
A  2
B  3
B  4
Then it will become this after running 'collapse.pl' with the defaults:
A  1,2
B  3,4

Note: the --function option is VERY USEFUL and lets you easily summarize your data (mean, etc.)
So make sure to check that out.

EXAMPLE:

cat myfile.pairs.txt | collapse.pl --dot > out.txt

OPTIONS:

-k COL: Compare the values in column COL to the threshold in the file (default is 1).

-i or --di DELIM: "Delimiter, input". The fields in the input file are seperated by DELIM (default is tab).

-o or --do DELIM: "Delimiter, output": Use this to seperate within a field when duplicates are found (default is comma).
	Note that there will ALWAYS be a tab to separate the key (first column) from the subsequent columns.
--dot: Make the output delimiter a TAB. Equivalent to --do=$'\t'

-f or --function FUNCTION: Instead of concatenating entries, apply the function
             FUNCTION to the vector (and print one result per key).
   Possible functions:
    -f mean    - Compute the mean   (synonym: "-f ave")
    -f min     - Compute the minimum
    -f max     - Compute the maximum
    -f median  - Compute the median
    -f sum     - Compute the sum
    -f count   - Compute the number of non-empty values
    -f std     - Compute the standard deviation
    -f var     - Compute the variance
    -f entropy - Compute the Shannon entropy


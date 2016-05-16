#!/usr/bin/perl

# By Alex Williams, 2016
# Usage: cat "fastaFile.fa" | fasta2gtf_one_line_per_chr.pl > out.gtf

use strict;  use warnings;  use diagnostics;
use POSIX      qw(ceil floor);
use List::Util qw(max min);
use Getopt::Long;
use Carp; # backtrace on errors. Has the "confess" function. Use this instead of "die" if you want useful information!   




sub main();
sub quitWithUsageError($) { print($_[0] . "\n"); printUsageAndQuit(); print($_[0] . "\n"); }
sub printUsageAndQuit() { printUsage(); exit(1); }

sub printUsage() {
    print STDOUT <DATA>;
    exit(0);
}




my $feature = undef;
my $len     = undef;

sub printFeatAsGtfLine($$) {
	my ($feat_name, $feat_len) = @_; # feature name, feature length
	if ($feat_len < 1) { confess "Got a FASTA feature named '$feat_name' with a length less than 1 (it was $feat_len), which is an error!"; }
	my $delim = "\t";
	my $source = "AGW_autogenerated_GTF";
	my $start = 1; # ("Start position of the feature, with sequence numbering starting at 1.": http://uswest.ensembl.org/info/website/upload/gff.html)
	my $end   = $start+$feat_len-1; # We do NOT add one to the length---in fact, we SUBTRACT one, because "1 to 1" is length 1, "1 to 2" is length 2, etc.  "<start> must be less than or equal to <end>"
	if ($end < $start) { confess "End cannot be less than start. Feature length must have been 0, which should be impossible."; }
	my $score = "0.00";
	my $strand = "+";
	# From GTF docs: http://mblab.wustl.edu/GTF22.html
	# The following feature types are required: "CDS", "start_codon", "stop_codon". The features "5UTR", "3UTR", "inter", "inter_CNS", "intron_CNS" and "exon" are optional. All other features will be ignored. The types must have the correct capitalization shown here.
	# The "start_codon" feature is up to 3bp long in total and is included in the coordinates for the "CDS" features. The "stop_codon" feature similarly is up to 3bp long and is excluded from the coordinates for the "3UTR" features, if used.
	my $ann = qq{gene_id "$feat_name"; transcript_id "$feat_name";};
	print STDOUT join($delim, "$feat_name", "$source", "CDS"        , "$start", "$end"  , "$score", "$strand", ".", $ann) . "\n";
	print STDOUT join($delim, "$feat_name", "$source", "exon"       , "$start", "$end"  , "$score", "$strand", ".", $ann) . "\n";
	print STDOUT join($delim, "$feat_name", "$source", "start_codon", "$start", "$start", "$score", "$strand", ".", $ann) . "\n";
	print STDOUT join($delim, "$feat_name", "$source", "stop_codon" , "$end"  , "$end"  , "$score", "$strand", ".", $ann) . "\n";
}

sub main() {
	$Getopt::Long::passthrough = 1; # ignore arguments we don't recognize in GetOptions, and put them in @ARGV

	GetOptions("help|?|man" => sub { printUsageAndQuit(); }
		  ) or printUsageAndQuit();

	my %feats = ();		# key = feature, value = length
	my $curFeature = undef;
	my $lineNum = 0;
	while (my $line = <>) {
		$lineNum++;
		chomp($line);
		if ($line =~ m/^>(.*)/) {
			$curFeature = $1;
			#if (defined($curFeature)) { printFeatAsGtfLine($feature, $len); } # print the previously-found feature
			if (exists($feats{$curFeature})) {
				confess "ERROR: Duplicate FASTA entries! There were multiple lines with the 'name' field (the part after the '>') as <$curFeature>, which is forbidden! Each FASTA name must be UNIQUE.";
			}
			$feats{$curFeature} = 0	# New feature found! $feature = $1;
			  #$len = 0; # reset it to zero!
		} else {
			if (!defined($curFeature)) {
				confess "ERROR: Odd... found something that is NOT in a feature, on line $lineNum...\n";
			}
			$line =~ s/\s//g; # remove **all** whitespace from this line!
			$feats{$curFeature} += length($line);
		}
	}

	foreach my $k (sort(keys(%feats))) { 
		printFeatAsGtfLine($k, $feats{$k});
	}
}


main();
exit(0);


# ====

__DATA__

fasta2gtf_one_line_per_chr.pl  INPUT.fasta   [ >   output.gtf]
   * Writes to STDOUT by default

by Alex Williams, 2016.

This program just makes a GTF file from an input FASTA file.
It makes one feature per '>' line in the fasta file.

See the examples below for more information.

CAVEATS:

Does not really follow the GTF specs properly, as it does not make a CDS start/end, for example.

OPTIONS:

  No options.

EXAMPLES:

fasta2gtf_one_line_per_chr.pl    my_seqs.fasta    >   out.gtf

KNOWN BUGS:

  None known.

--------------








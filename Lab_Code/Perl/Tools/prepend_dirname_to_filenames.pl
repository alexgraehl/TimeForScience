#!/usr/bin/perl -w
#

# Prepends the directory name to any filenames matching this pattern. Useful for uniquely-named directories with tons of identically-named files.

use strict;

use Getopt::Long;
use Cwd;
use File::Spec;
#Getopt::Long::Configure('bundling');

my ($verbose, $dryrun, $force);

die "\"prepend_dirname_to_filenames\" requires input file arguments (NOT directory arguments!).\nUsage: rename [-v] [-n] [-f] [filenames]\nExample: prepend_dirname_to_filenames.pl */*.txt\n\n"
    unless GetOptions(
	'v|verbose' => \$verbose,
	'n|d|no-act|dry|dry_run|dryrun|dry-run'  => \$dryrun,
	'f|force'   => \$force,
    );

if ($dryrun) { $verbose = 1; } ## always be verbose when we are using a dry-run

if (!@ARGV) {
    print STDERR "Since no filenames were specified, we are reading filenames directly from STDIN.\n";
    @ARGV = <STDIN>;
    chop(@ARGV);
}


if ($verbose) { print "About to examine " . scalar(@ARGV) . " items that might be renamed if they match the pattern. . .\n"; }

for my $f (@ARGV) {
    my @splitp = File::Spec->splitpath($f);
    my $dir = $splitp[1];
    
    my $noslash = $f; $noslash =~ s/[\/]/_/g;

    my $newf = "${dir}${noslash}";

    if ($dryrun) { print "Dry run: "; }
    print "$f --> ${newf}" . "\n"; 

    if (-e $newf and !$force) {
	warn "<$f> was NOT renamed, because <${newf}> already exists\n";
    } else {
	if (!$dryrun) {
	    if (!rename($f, $newf)) { ## <-- now try to actually rename the file on disk!
		warn "ERROR: Can't rename $f to $newf: $!\n"; ## any error codes are stored in "$!" automatically
	    }
	}
    }

}


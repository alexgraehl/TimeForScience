#!/usr/bin/env perl

#use lib "$ENV{MYPERLDIR}/lib"; use lib "$ENV{TIME_FOR_SCIENCE_DIR}/Lab_Code/Perl/LabLibraries"; require "libfile.pl";
#use lib "$ENV{MYPERLDIR}/lib"; use lib "$ENV{TIME_FOR_SCIENCE_DIR}/Lab_Code/Perl/LabLibraries"; require "libstats.pl";
#use lib "$ENV{MYPERLDIR}/lib"; use lib "$ENV{TIME_FOR_SCIENCE_DIR}/Lab_Code/Perl/LabLibraries"; require "libstring.pl";
#use lib "$ENV{MYPERLDIR}/lib"; use lib "$ENV{TIME_FOR_SCIENCE_DIR}/Lab_Code/Perl/LabLibraries"; require "libsystem.pl";

#use List::Util 'shuffle';
#@shuffled = shuffle(@list);
# Check out: perldoc -q array

use strict;  use warnings;  use diagnostics;
use POSIX      qw(ceil floor);
use List::Util qw(max min);
use Getopt::Long;
use Carp; # backtrace on errors. Has the "confess" function. Use this instead of "die" if you want useful information!

#use File::Basename;

$| = 1; # Always flush text output IMMEDIATELY to the console, don't wait to buffer terminal output! Setting this to zero can cause STDERR and STDOUT to be interleaved in weird ways.

#no warnings 'numeric';
#use Scalar::Util;
#print Scalar::Util::looks_like_number($string), "\n";

sub tryToLoadModule($) {
    my $x = eval("require $_[0]");
    if ((defined($@) && $@)) {
	warn "We FAILED to load module $_[0]. Skipping this module, but continuing with the program.";
	return 0; # FAILURE
    } else {
	$_[0]->import();
	return 1; # SUCCESS
    }
}

my $SHOULD_USE_COLORS = tryToLoadModule("Term::ANSIColor");
if ($SHOULD_USE_COLORS) { use Term::ANSIColor; }

sub warnPrint($) { chomp($_[0]); warn(safeColor("[WARNING]: " . $_[0] . "", "yellow on_black")); } # regarding "warn": if the string ends with a newline it WON'T print the line number!

sub safeColor($;$) { # one required and one optional argument
    ## Returns colored text, but only if $SHOULD_USE_COLORS is set.
    ## Allows you to totally disable colored printing by just changing $SHOULD_USE_COLORS to 0 at the top of this file
    # Colorstring is OPTIONAL, and can be something like "red on_blue" or "red" or "magenta on_green"
    # Example usage:
    #    *    print STDERR safeColor("This warning message is red on yellow", "red on_yellow");
    my ($message, $color) = @_;
    return (($SHOULD_USE_COLORS && defined($color)) ? (Term::ANSIColor::colored($message, $color) . Term::ANSIColor::color("reset")) : $message);
}

sub printColorStderr($;$) {
    # prints color to STDERR *UNLESS* it is re-directed to a file, in which case NO COLOR IS PRINTED.
    my ($msg, $col) = @_; # Only prints in color if STDERR is to a terminal, NOT if it is redirected to an output file!
    if (! -t STDERR) { $col = undef; } # no coloration if this isn't to a terminal
    print STDERR safeColor($msg, $col);
}

sub printColorStdout($;$) {
    # prints color to STDOUT *UNLESS* it is re-directed to a file, in which case NO COLOR IS PRINTED.
    my ($msg, $col) = @_; # Only prints in color if STDOUT is to a terminal, NOT if it is redirected to an output file!
    if (! -t STDOUT) { $col = undef; } # no coloration if this isn't to a terminal
    print STDOUT safeColor($msg, $col);
}

sub dryNotify(;$) { # one optional argument
    my ($msg) = @_;
    $msg = (defined($msg)) ? $msg : "This was only a dry run, so we skipped executing a command.";
    print STDERR safeColor("[DRY RUN]: $msg\n", "black on_yellow");
}

sub notify($) { # one required argument
    my ($msg) = @_;
    warn safeColor("[DRY RUN]: $msg\n", "cyan on_blue");
}

sub main();

sub printUsage()        { print STDOUT <DATA>; }
sub printUsageAndQuit() { printUsage(); exit(1); }
sub quitWithUsageError($) { print($_[0] . "\n"); printUsageAndQuit(); print($_[0] . "\n"); }

# ==1==
sub main() { # Main program
    my ($delim) = "\t";
    my ($decimalPlaces) = 4; # How many decimal places to print, by default
    $Getopt::Long::passthrough = 1; # ignore arguments we don't recognize in GetOptions, and put them in @ARGV

    GetOptions("help|?|man" => sub { printUsageAndQuit(); }
	       , "delim|d=s" => \$delim
	       , "dp=i" => \$decimalPlaces
	) or printUsageAndQuit();

    print STDERR "===============================\n";
    print STDOUT safeColor("test color\n", "blue on_yellow");
    printColorStdout("hey", "red");
    printColorStdout("what", "red");
    printColorStderr("stderr what", "blue on_white");
    printColorStdout("is", "red on_blue");
    printColorStdout("this", "red on_green");
    printColorStderr("stderr what", "red on_white");
    print STDERR "hellow\n";
    print STDERR "Test color!\n";
    print STDERR "===============================\n";

    if (1 == 0) {
	quitWithUsageError("1 == 0? Something is wrong!");
    }

    my $numUnprocessedArgs = scalar(@ARGV);
    if ($numUnprocessedArgs != 2) {
	quitWithUsageError("Error in arguments! You must send TWO filenames to this program.\n");
    }

    my $filename1 = undef;
    my $filename2 = undef;

    foreach (@ARGV) { # these were arguments that were not understood by GetOptions
	if (!defined($filename1)) { $filename1 = $_; }
	elsif (!defined($filename2)) { $filename2 = $_; }
	else {
	    print STDERR "Unprocessed argument: $_\n";
	}
    }

    # You can use 'confess' instead of 'die' if you want a backtrace!

    #print "\t" . sprintf("%.${numDecimalPointsToPrint}f", $levFraction);
    #my @col1 = @{readFileColumn($filename1, 0, $delim)};

    print STDERR safeColor("===============================\n", "green");
    #print STDERR "Done!\n";
} # end main()


main();


END {
    # Runs after everything else.
    # Makes sure that the terminal text is back to its normal color.
    if ($SHOULD_USE_COLORS) { print STDERR Term::ANSIColor::color("reset"); print STDOUT Term::ANSIColor::color("reset"); }
}

exit(0);
# ====

__DATA__

MYPROGRAM.pl  [OPTIONS]

by YOUR NAME, THE_YEAR

THIS PROGRAM DOES SOMETHING. YOU GIVE IT ONE THING AND THIS OTHER THING,
AND THE OUTPUT IS THIS FINAL THING.

See the examples below for more information.

CAVEATS:

MAYBE IT TAKES 30 MINUTES TO RUN.

OPTIONS:

  --delim = DELIMITER   (Default: tab)
     Sets the input delimiter to DELIMITER.

EXAMPLES:

MYPROGRAM.pl --help
  Displays this help

MYPROGRAM.pl  --works=yes --bugs=4  -q
  Does nothing. -q indicates "quiet" operation.


KNOWN BUGS:

  None known.

TO DO:

  Add ???.

--------------



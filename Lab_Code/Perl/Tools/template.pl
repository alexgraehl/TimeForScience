#!/usr/bin/perl

#use lib "$ENV{MYPERLDIR}/lib"; use lib "$ENV{TIME_FOR_SCIENCE_DIR}/Lab_Code/Perl/LabLibraries"; require "libfile.pl";
#use lib "$ENV{MYPERLDIR}/lib"; use lib "$ENV{TIME_FOR_SCIENCE_DIR}/Lab_Code/Perl/LabLibraries"; require "libstats.pl";
#use lib "$ENV{MYPERLDIR}/lib"; use lib "$ENV{TIME_FOR_SCIENCE_DIR}/Lab_Code/Perl/LabLibraries"; require "libstring.pl";
#use lib "$ENV{MYPERLDIR}/lib"; use lib "$ENV{TIME_FOR_SCIENCE_DIR}/Lab_Code/Perl/LabLibraries"; require "libsystem.pl";

#use List::Util 'shuffle';
#@shuffled = shuffle(@list);
# Check out: perldoc -q array

use POSIX      qw(ceil floor);
use List::Util qw(max min);

sub tryToLoadModule($) {
    my $x = eval("require $_[0]");
    if ((defined($@) && $@)) {
	warn "Module loading of $_[0] FAILED. Skipping this module.";
	return 0; # FAILURE: return 0(false)
    } else {
	$_[0]->import();
	return 1; # SUCCESS: return 1(true)
    }
}

my $SHOULD_USE_COLORS = tryToLoadModule("Term::ANSIColor");

sub warnPrint($) { chomp($_[0]); warn(safeColor("[WARNING]: " . $_[0] . "", "yellow on_black")); } # regarding "warn": if it ends with a newline it WON'T print the line number

sub safeColor($;$) { # one required and one optional argument
    ## Prints colored text, but only if USER_COLORS_CONSTANT is set.
    ## Allows you to totally disable colored printing by just changing USE_COLORS_CONSTANT to 0 at the top of this file
    my ($str, $color) = @_;
    return ((USE_COLORS_CONSTANT) ? colored($str, $color) : $str);
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



if (USE_COLORS_CONSTANT) {
    use Term::ANSIColor;
}

#use File::Basename;
use Getopt::Long;

#no warnings 'numeric';
#use Scalar::Util;
#print Scalar::Util::looks_like_number($string), "\n";

use strict;  use warnings;  use diagnostics;

sub main();

#print colorString("blue");
#print "Arr";
#colorResetString();

sub safeColor($$) {
    my ($msg, $colorString) = @_;
    # Colorstring should be something like "red on_blue" or "red" or "magenta on_green"
    if (USE_COLORS_CONSTANT) {
	return(Term::ANSIColor::colored($msg, $colorString));
    } else {
	return($msg); # no color support apparently
    }
}

sub colorString($) {
    # Requires "use Term::AnsiColor". Note: you have to PRINT the result of this function!
    # It also only sets the output color if the output is a TERMINAL.
    # If the output is NOT a terminal, then colorizing output
    # results in lots of garbage characters (color control characters) written to the screen,
    # so we don't do it.
    # Example usage: print colorString("red"); print "something red"; print colorString("reset");
    return ((-t STDOUT && USE_COLORS_CONSTANT) ? (Term::ANSIColor::color($_[0])) : ""); # <-- checks to see if STDOUT goes directly to the terminal (instead of, say, outputting to a file with a redirect)
    ## probably should actually use "colored( ... , "red")" as an example
}

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

    print STDERR colorString("green");
    print STDERR "===============================\n";
    print STDERR "Test color!\n";
    print STDERR "===============================\n";
    print STDERR colorString("reset");

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

    #print "\t" . sprintf("%.${numDecimalPointsToPrint}f", $levFraction);
    #my @col1 = @{readFileColumn($filename1, 0, $delim)};

    print STDERR colorString("green");
    print STDERR "===============================\n";
    print STDERR "Done!\n";
    print STDERR colorString("reset");
} # end main()


main();


END {
    # Runs after everything else.
    # Makes sure that the terminal text is back to its normal color.
    print STDERR colorString("reset");
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

#!/usr/bin/perl

use strict;  use warnings;
use Getopt::Long;

sub main();
sub quitWithUsageError($) { print($_[0] . "\n"); printUsageAndQuit(); print($_[0] . "\n"); }
sub printUsageAndQuit() { printUsage(); exit(1); }
sub printUsage() { print STDOUT <DATA>; exit(0); }

# ==1==
sub main() { # Main program
    my ($delim) = ":"; # default delim for a path is ':'
    $Getopt::Long::passthrough = 1; # ignore arguments we don't recognize in GetOptions, and put them in @ARGV
    GetOptions("help|h|man" => sub { printUsageAndQuit(); }
	       , "delim|d=s" => \$delim
	      ) or printUsageAndQuit();

    my $numUnprocessedArgs = scalar(@ARGV);
    ($numUnprocessedArgs == 1) or quitWithUsageError("You need to leave ONE argument for this script--a path! Preferably in quotes, unless you want spaces to be mangled.");

    my $arg = $ARGV[0];
    chomp($arg);

    my @items = split($delim, $arg);
    my %seen = ();
    my @newItems = ();
    for my $x (@items) {
	    if (!defined($seen{$x})) { push(@newItems, $x); }
	    $seen{$x} = 1;
    }
    print STDOUT join($delim, @newItems);
    print STDOUT "\n";
} # end main()

main();
END { }
exit(0);
# ====

__DATA__

dedupe_path.pl [SOME PATH HERE]

by Alex Williams, 2016.

You give it ONE argument (a path), and it removes duplicate entries after the first one.

Make sure to quote your inputs. May behave unpredictably with spaces, just like any UNIX tool.

kOPTIONS:

  --delim or -d = DELIMITER   (Default: ':')
     Sets the input delimiter to DELIMITER. Normally a ':' for UNIX paths.

EXAMPLES:

dedupe_path.pl "/bin:/bin:/home/bin/:" <-- note the two "/bin" entries
Result:   "/bin:/home/bin/"  <-- note that the redundant "/bin" was taken out

You might want to use it like this, to prevent your UNIX path from getting super long if you re-source your ~/.bashrc.

PATH=$(dedupe_path.pl  /some/new/stuff/and/also:/bin/whatever:$PATH)

--------------

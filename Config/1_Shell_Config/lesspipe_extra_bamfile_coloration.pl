#!/usr/bin/perl
use strict;
use warnings;

# BEWARE WHEN CHANGING THIS PROGRAM:
# Note: if you make any errors in this file, "less" will stop coloring fastq files WITHOUT ANY ERROR REPORTING.
# You will have to actually run 'perl -c' with this program to check it, or better yet, actually debug it with real input files.

my $STO=qq{\033[37m} ; # "default" white color. Assumes dark terminal!
my $CCC=qq{\033[44m\033[1;37m}; # colon delimiter color; 34 = blue (NOT USED RIGHT NOW)

# The most common standard for colors (in sequence logos) is: A = green, C = blue, T = red, G = yellow

# F_COLOR = FOREGROUND IS THIS COLOR
my $F_RED    = qq{\033[31m};
my $F_GREEN  = qq{\033[32m};
my $F_YELLOW = qq{\033[33m};
my $F_CYAN   = qq{\033[36m};

my $STA=$F_GREEN; # start color A: 32 = green 31 = red  34 = blue
my $STT=$F_RED; # start color T: 36 = cyan. 35 = magenta
my $STC=$F_CYAN; # start color C: 36 = cyan; 31 = red, 36 = cyan. 32 = green. 35 = magenta
my $STG=$F_YELLOW; # start color G: 33 = yellow

# Coloring the BACKGROUND instead of the foreground (this sounds like it would work well, but it looks really ugly!)
#my $STA=qq{\033[42m}; # start color A: 32 = green 31 = red  34 = blue
#my $STT=qq{\033[41m}; # start color T: 36 = cyan. 35 = magenta
#my $STC=qq{\033[44m}; # start color C: 36 = cyan; 31 = red, 36 = cyan. 32 = green. 35 = magenta
#my $STG=qq{\033[43m}; #  \033[45m}; # start color N: magenta background, black foreground

# Note: currently we assume that ONLY THE 'N' has a background
# color. If that is not true, then we need to change the
# line with /[Nn]/ to push the 'reset' color for other bases
# with a background color. Otherwise the background "leaks"
# beyond where it is supposed to be.

my $STN=qq{\033[45m\033[1;37m}; #  \033[45m}; # start color N: magenta background, black foreground
my $RES=qq{\033[0m} ; # RESET color

my $MIN_DNA_SEQ_LENGTH_TO_COLOR = 5; # Color a run of "ACGTN-" if they are at least THIS long (or longer).

my %colors = ( 'A' => $STA
	       , 'a' => $STA
	       , 'T' => $STT
	       , 't' => $STT
	       , 'G' => $STG
	       , 'g' => $STG
	       , 'C' => $STC
	       , 'c' => $STC
	       , 'N' => $STN
	       , 'n' => $STN
	       , '-' => $RES # Hyphen does not get colored--reset it to the terminal default!
	       
);

while (my $line = <>) {
    my @potentiallyInteresting = ();
    my @boringVersion = ();
    my $inDNA = 0;

    my @arr = split(//, $line);

    my $prevC = '';
    for my $c (@arr) {
	if ($c =~ /[-ACGTNacgtn]/) {
	    # Hyphen (gap) is allowed to be "maybe" part of a DNA sequence. It gets whatever color the previous base had, if any.
	    # We could fix this by resetting the color after every base, if we really wanted to.
	    $inDNA = 1;

	    my $toPush;
	    if ($prevC eq $c) {
		$toPush = $c; # if the previous character was ALSO the same as this one, then no need to change the colors again... unless it's an N
	    } else {
		$toPush = $colors{$c} . $c;
	    }
	    if ($c =~ /[Nn]/) { $toPush .= $RES; } # note that we add the $RES reset character AFTER AN 'N' but not anything else. That's because N needs its background reset, as it is the only letter with a colored background.
	    
	    push(@potentiallyInteresting, $toPush);
	    push(@boringVersion, $c);
	} else {
	    if ($inDNA) {
		if (scalar(@potentiallyInteresting) >= $MIN_DNA_SEQ_LENGTH_TO_COLOR) {
		    print join('',@potentiallyInteresting) . $RES;
		} else {
		    print join('',@boringVersion);
		}
		@potentiallyInteresting = @boringVersion = (); # clear it!
		$inDNA = 0;
	    }
	    #if ($c eq ':') { print "$CCC:$RES"; }
	    #else { print $c; }
	    print $c;
	}
    }
    
}

#!/usr/bin/perl
use strict;
use warnings;

$| = 1; # <-- ALWAYS print directly out, NEVER buffer!

my $STO=qq{\033[1m} ; # start color

my $CCC=qq{\033[44m\033[1;37m}; # colon delimiter color; 34 = blue

my $STA=qq{\033[32m}; # start color A: 32 = green   31 = red  34 = blue
my $STT=qq{\033[31m}; # start color T: 36 = cyan. 35 = magenta
my $STC=qq{\033[36m}; # start color C: 36 = cyan; 31 = red, 36 = cyan. 32 = green. 35 = magenta
my $STG=qq{\033[33m}; # start color G: 33 = yellow

my $STN=qq{\033[45m\033[1;37m}; #  \033[45m}; # start color N: magenta background, black foreground
my $RES=qq{\033[0m} ; # RESET color

my $MIN_DNA_SEQ = 5; # Color a run of ACGTN if they are at least THIS long (or longer).

my $MAX_LINE_LENGTH_TO_COLOR = 7; # In order to keep 'less' from getting too slow, only specially process the first 500bp of a line. After that, no coloration. This keeps things from breaking on long .fasta files!

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
);

my $lineNum = 0;
while (my $line = <>) {
    $lineNum++;
    my @potentiallyInteresting = ();
    my @boringVersion = ();
    my $inDNA = 0;
    my $lenThisLine = 0;
    my @arr = split(//, $line, $MAX_LINE_LENGTH_TO_COLOR); # arr can only be at most 500 bp

    #print join("|", @arr);
    my $prevC = '';
    my $arlast = (scalar(@arr)-1);
    for (my $i = 0; $i <= $arlast; $i++) {
	my $c = $arr[$i]; #for my $c (@arr) {
	$lenThisLine++;
	my $stillHaveMoreToProcessAfterThisChar = ($i < $arlast);
	if ($stillHaveMoreToProcessAfterThisChar && $c =~ /[ACGTNacgtn]/) { #$MAX_LINE_LENGTH_TO_COLOR) {
	    $inDNA = 1;
	    my $toPush;
	    if ($prevC eq $c) {
		$toPush = $c; # if the previous character was ALSO the same as this one, then no need to change the colors again... unless it's an N
	    } else {
		$toPush = $colors{$c} . $c;
	    }
	    if ($c =~ /[Nn]/) { $toPush .= ${RES}; } # note that we add the $RES reset character AFTER AN 'N' but not anything else. That's because N needs its background reset, as it is the only letter with a colored background.
	    push(@potentiallyInteresting, $toPush);
	    push(@boringVersion, $c);
	} else {
	    if ($inDNA) { # was previously in a colored DNA set, so let's clear it!
		if (scalar(@potentiallyInteresting) >= $MIN_DNA_SEQ) {
		    print join('',@potentiallyInteresting) . $RES;
		} else {
		    print join('',@boringVersion);
		}
		@potentiallyInteresting = @boringVersion = (); # clear it!
		$inDNA = 0;
	    }
	    print $c;
	}
    }

    if ($lenThisLine >= $MAX_LINE_LENGTH_TO_COLOR) {
	print substr($line, $MAX_LINE_LENGTH_TO_COLOR); # print the rest of the line! ("substr(STRING, OFFSET, LENGTH)")
    }
    #print "\n";
}

#!/usr/bin/perl
use strict;
use warnings;

#$| = 1; # <-- ALWAYS print directly out, NEVER buffer!

my $STO=qq{\033[1m} ; # start color

my $CCC=qq{\033[44m\033[1;37m}; # colon delimiter color; 34 = blue

my $STA=qq{\033[32m}; # start color A: 32 = green   31 = red  34 = blue
my $STT=qq{\033[31m}; # start color T: 36 = cyan. 35 = magenta
my $STC=qq{\033[36m}; # start color C: 36 = cyan; 31 = red, 36 = cyan. 32 = green. 35 = magenta
my $STG=qq{\033[33m}; # start color G: 33 = yellow

my $STN=qq{\033[45m\033[1;37m}; #  \033[45m}; # start color N: magenta background, black foreground
my $RES=qq{\033[0m} ; # RESET color

my $MIN_DNA_SEQ_LENGTH_TO_COLOR = 10; # Color a run of ACGTN if they are at least THIS long (or longer).

my $MAX_LINE_LENGTH_TO_COLOR = 10; # In order to keep 'less' from getting too slow, only specially process the first 500bp of a line. After that, no coloration. This keeps things from breaking on long .fasta files!

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

my @potentiallyInteresting = ();
my @boringVersion = ();
my $inDNA = 0;
my $lenThisLine = 0;
my $numColoredBases = 0;
my $prevC = ''; # previous color

open(AFILE, "-");
my $c;
while (sysread(AFILE, $c, 1)) {
#while(my $c = getc) {
    #print "$c.";
    if ($c eq "\n") {
	print "RESETTING\n";
	@potentiallyInteresting = ();
	@boringVersion = ();
	$inDNA = 0;
	$lenThisLine = 0;
	$numColoredBases = 0;
	$prevC = ''; # previous color
	next;
    }
    $lenThisLine++;
    if (1 == 0 && $lenThisLine > $MAX_LINE_LENGTH_TO_COLOR) {
	if ($inDNA) {
	    if (scalar(@potentiallyInteresting) >= $MIN_DNA_SEQ_LENGTH_TO_COLOR) {
		print join('',@potentiallyInteresting) . $RES;
	    } else {
		print join('',@boringVersion);
	    }
	}
	#print substr($line, $lenThisLine); # print the rest of the line! ("substr(STRING, OFFSET, LENGTH)")
	#last; # <-- Go to the NEXT LINE!
    }

    if ($c =~ /[ACGTNacgtn]/) {
	next;
	$inDNA = 1;
	my $toPush;
	if ($prevC eq $c) {
	    $toPush = $c; # if the previous character was ALSO the same as this one, then no need to change the colors again... unless it's an N
	} else {
	    $toPush = $colors{$c} . $c;
	}
	if ($c =~ /[Nn]/) { $toPush .= ${RES}; } # note that we add the $RES reset character AFTER AN 'N' but not anything else. That's because N needs its background reset, as it is the only letter with a colored background.
	
	push(@potentiallyInteresting, $toPush);
	push(@boringVersion, $c); # this is really not necessary, we could use substring instead
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
    #}
    
}

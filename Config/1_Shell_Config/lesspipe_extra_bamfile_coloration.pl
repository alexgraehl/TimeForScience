#!/usr/bin/perl
use strict;
use warnings;

# BEWARE WHEN CHANGING THIS PROGRAM:
# Note: if you make any errors in this file, "less" will stop coloring fastq files WITHOUT ANY ERROR REPORTING.
# You will have to actually run 'perl -c' with this program to check it, or better yet, actually debug it with real input files.

my $filename = "NONE_SPECIFIED";
if (scalar(@ARGV) == 0) {
	warning("This script expects a filename (or blank placeholder argument) to be passed in.");
} else {
	$filename = $ARGV[0]; # first argument
}

my $CCC=qq{\033[44m\033[1;37m}; # colon delimiter color; 34 = blue (NOT USED RIGHT NOW)

# The most common standard for colors (in sequence logos) is: A = green, C = blue, T = red, G = yellow

# F_COLOR = FOREGROUND IS THIS COLOR
my $STO=qq{\033[37m} ; # "default" white color. Assumes dark terminal!???
my $F_DEFAULT = qq{\033[00m}; # <-- or is THIS the default?? 
my $F_RED     = qq{\033[31m};
my $F_GREEN   = qq{\033[32m};
my $F_YELLOW  = qq{\033[33m};
my $F_BLUE    = qq{\033[34m};
my $F_MAGENTA = qq{\033[35m};
my $F_CYAN    = qq{\033[36m};

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

#my $INPUT_FILENAME = undef;
#if (scalar(@ARGV) == 1) {
#	$INPUT_FILENAME = $ARGV[0];
#} elsif (scalar(@ARGV) >= 2) {
#	die "Error in lesspipe_extra_bamfile_coloration.pl: more than two input filenames!"
#}

my $STN=qq{\033[45m\033[1;37m}; #  \033[45m}; # start color N: magenta background, black foreground
my $RES=qq{\033[0m} ; # RESET color

my $MIN_DNA_SEQ_LENGTH_TO_COLOR = 5; # Color a run of "ACGTN-" if they are at least THIS long (or longer).

my %dna_colors = ( 'A' => $STA  , 'a' => $STA
	       , 'T' => $STT, 't' => $STT
	       , 'G' => $STG, 'g' => $STG
	       , 'C' => $STC, 'c' => $STC
	       , 'N' => $STN, 'n' => $STN
	       , '-' => $RES # Hyphen does not get colored--reset it to the terminal default!
);

# Amino acid colors!
my $ASP = $F_RED;   my $GLU = $F_RED;
my $LYS = $F_BLUE;  my $ARG = $F_BLUE;
my $PHE = $F_CYAN;  my $TYR = $F_CYAN;
my $GLY = $F_DEFAULT;
my $ALA = $F_BLUE;
my $HIS = $F_CYAN;
my $CYS = $F_YELLOW; my $MET = $F_YELLOW;
my $LEU = $F_GREEN;  my $VAL = $F_GREEN; my $ILE = $F_GREEN;
my $SER = $F_YELLOW; my $THR = $F_YELLOW;
my $ASN = $F_GREEN;  my $GLN = $F_GREEN;
my $TRP = $F_RED;
my $PRO = $F_RED;
my $XXX = $STN; # custom background color for N and X

# ASP,GLU   bright red [230,10,10]     CYS,MET     yellow [230,230,0]
#  LYS,ARG   blue       [20,90,255]     SER,THR     orange [250,150,0]
#  PHE,TYR   mid blue   [50,50,170]     ASN,GLN     cyan   [0,220,220]
#  GLY       light grey [235,235,235]   LEU,VAL,ILE green  [15,130,15]
#  ALA       dark grey  [200,200,200]   TRP         pink   [180,90,180]
#  HIS       pale blue  [130,130,210]   PRO         flesh  [220,150,130]
my %aa_colors = ( 'F' => $PHE
		  , 'L' => $LEU
		  , 'I' => $ILE
		  , 'M' => $MET
		  , 'V' => $VAL
		  , 'S' => $SER
		  , 'P' => $PRO
		  , 'T' => $THR
		  , 'A' => $ALA
		  , 'Y' => $TYR
		  , 'H' => $HIS
		  , 'Q' => $GLN
		  , 'N' => $ASN
		  , 'K' => $LYS
		  , 'D' => $ASP
		  , 'E' => $GLU
		  , 'C' => $CYS
		  , 'W' => $TRP
		  , 'R' => $ARG
		  , 'G' => $GLY
		  , 'X' => $XXX, '*' => $XXX, '-' => $XXX
		);
my $aa_candidate_string = join("", sort(keys(%aa_colors)));

# Try to guess if a section of text is DNA.
# (Maybe also try to guess if it's an amino acid? (AA)? Not currently implemented.)

my $is_fastq = ($filename =~ m/[.](fastq|fq)/i);

#print("FASTQ status for file ${filename} is: " . ($is_fastq ? "IS A FASTQ" : "NOT A FASTQ") . "\n");
while (my $line = <STDIN>) {
	my @potentiallyInteresting = ();
	my @boringVersion = ();
	my $inDNA = 0; # DNA bases
	my $inAA  = 0; # amino acids!
	my @arr = split(//, $line);
	#my $is_amino_acids_only_no_spaces_entire_line = ($line =~ m/^[${aa_candidate_string}]+$/);
	#print "Maybe --> $may_be_aa_only_line ($aa_candidate_string)" . "\n\n";
	my $prevC = '';
	for my $c (@arr) {
		if ($c =~ m/[-ACGTNacgtn]/) {# or $is_amino_acids_only_no_spaces_entire_line) {
			# Hyphen (gap) is allowed to be "maybe" part of a DNA sequence. It gets whatever color the previous base had, if any.
			# We could fix this by resetting the color after every base, if we really wanted to.
			$inDNA = 1;
			#my $toPush = "a"; #"$aa_colors[$c]$c";
			my $toPush = ($prevC eq $c) ? $c : "$dna_colors{$c}$c"; # <-- if the previous character was ALSO the same as this one, then no need to change the colors again... unless it's an N?
			if (($c eq 'N') or ($c eq 'n')) { $toPush .= $RES; } # note that we add the $RES reset character AFTER AN 'N' but not anything else. That's because N needs its background reset, as it is the only letter with a colored background.
			push(@potentiallyInteresting, $toPush);
			push(@boringVersion, $c);
		} else {
			if ($inDNA) {
				if (scalar(@potentiallyInteresting) >= $MIN_DNA_SEQ_LENGTH_TO_COLOR) {
					print join('',@potentiallyInteresting) . $RES; # $RES is RESET COLOR, which we do after each block of characters (NOT individual characters!)
				} else {
					print join('',@boringVersion);
				}
				@potentiallyInteresting = @boringVersion = (); # clear it!
				$inDNA = 0;
			}
			#if ($c eq ':') { print "$CCC:$RES"; }
			#else { print $c; }
			print $c;

			# does this properly handle the last element of each line, or is it RELYING on the newline to detect line-is-done?
		}
    }
}

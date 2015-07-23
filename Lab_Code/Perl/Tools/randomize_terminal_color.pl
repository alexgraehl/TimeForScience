#!/usr/bin/perl -w
use strict;

#@COMMENT@ randomize_terminal_color.pl only works with the Apple Mac OS X Terminal.app. Does NOT do anything on Linux or in iTerm on the Mac. It changes the background of your terminal when invoked with no arguments. Frequency-of-use rating: 1/10.

### Alex Williams' Terminal.app color randomizer, with minor additions by Charlie Vaske
## Usage:  ./randomize_terminal_color.pl -cycle   or randomize_terminal_color.pl -generate 1

if (not (exists($ENV{'TERM_PROGRAM'}) && ($ENV{'TERM_PROGRAM'} eq "Apple_Terminal"))) {
	print STDERR "randomize_terminal_color.pl only works on the Mac OS X Terminal.app. Exiting now.\n";
	exit(0);
}

# parameters which determine the range of random colors chosen
my  ($min_h,$max_h) =  ( 0,    360   );   # hue        (0-360) (degrees of angle on color wheel)
my  ($min_s,$max_s) =  ( 0.90, 0.91  );   # saturation (0-1) ("how washed out" things are)
my  ($min_v,$max_v) =  ( 0.10, 0.11  );   # brightness (0-1)
my  ($min_a,$max_a) =  ( 1.0,  1.0  );    # alpha      (0-1) (1.0 = opaque, to 0.0 = transparent)

# Order of colors is: cursor, text, bold, background. Specified in HUE (0-360), SATURATION(0-1), and VALUE(0-1)
my @profiles = ([ [0,0,1/3], [0,0,0], [0,0,0], [0,0,1] ], # black on white
		[ [0,0,1],   [0,0,1], [0,0,1], [0,0,0] ], # white on black
		[ [120,1,1], [120,1,1], [120,1,1], [0,0,0] ], # gren on blck
		[ [0,0,1/3], [0,0,0], [0,0,0], [60,0.29,1] ]  # blck on ylw
		);

my @colors_for_cycle = (0, 45, 140, 200, 240, 290);
#my @colors_for_cycle = (0, 60, 120, 180, 240, 300, 30, 90, 150, 210, 270);
#30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330);

sub oppositeHue {
    my ($inputHue) = @_;
    my $maxHue = 360;
    return ($maxHue/2.0 - $inputHue) + (($inputHue > ($maxHue/2.0)) ? $max_h : 0);
}

my @cycle_profiles = ();
foreach my $c (@colors_for_cycle) {
    # Order of colors is: cursor, text, bold, background. Specified in HUE (0-360), SATURATION(0-1), and VALUE(0-1)
    my $cursorSaturation = 0.70;
    my $cursorValue      = 1.0;
    my $textSat = 0.15;
    my $textVal = 1.00;
    my $boldSat = 0.50;
    my $boldVal = 1.00;
    my $bgSat   = 0.90;
    my $bgVal   = 0.10;
    push (@cycle_profiles, [ [oppositeHue($c), $cursorSaturation, $cursorValue], ## <-- Cursors
			     [$c, $textSat, $textVal],   ## <-- Regular text
			     [$c, $boldSat, $boldVal],   ## <-- Bold text
			     [$c, $bgSat,   $bgVal]  ## <-- Background
	  ]);
}

sub randbetween { my ($min,$max) = @_; return $min + rand($max-$min); } # generate random number between two values

# HSV to RGB conversion algorithm from  http://www.cs.rit.edu/~ncs/color/t_convert.html
sub hsv2rgb {
    my ($h,$s,$v) = @_;
    if ($s == 0) { return ($v,$v,$v); }  # handle greyscale case
    my ($i, $f, $p, $q, $t);
    $h /= 60;  # convert to sector between 0 and 5
    $i = int($h);
    $f = $h - $i;
    $p = $v * (1-$s);
    $q = $v * (1-$s*$f);
    $t = $v * (1-$s*(1-$f));
    return ($v, $t, $p) if $i == 0;
    return ($q, $v, $p) if $i == 1;
    return ($p, $v, $t) if $i == 2;
    return ($p, $q, $v) if $i == 3;
    return ($t, $p, $v) if $i == 4;
    return ($v, $p, $q);
}

sub rgb2string { return join ", ", map(int(65535*$_),@_); } # Turns a list of floats into a string suitable for applescript

#main program

my ($bg_color, $text_color, $bold_color, $cursor_color);
if (scalar(@ARGV) == 0) {
    my $randHue = randbetween($min_h, $max_h);
    my @rgb;
    my $alpha = int(65535); #int(65535*randbetween($min_a, $max_a));
    
    #print $randHue . " is the random hue.\n";

    my $oppositeHue = oppositeHue($randHue);

    #print "Cursor (opposite) hue is: " . $oppositeHue . "\n";

    @rgb = hsv2rgb(($max_h - $randHue), 0.70, 1.0); # Hue Saturation Value(lightness)
    $cursor_color = '{ ' . rgb2string(@rgb). ' }';
    
    @rgb = hsv2rgb($randHue, 0.05, 1.0); # Hue Saturation Value(lightness)
    $text_color = '{ ' . rgb2string(@rgb). ' }';
    
    @rgb = hsv2rgb($randHue, 0.50, 1.0); # Hue Saturation Value(lightness)
    $bold_color = '{ ' . rgb2string(@rgb). ' }';
    
    @rgb = hsv2rgb($randHue, randbetween($min_s,$max_s), randbetween($min_v,$max_v));
    $bg_color   = '{ ' . rgb2string(@rgb). ", $alpha }";
    
} elsif ($ARGV[0] eq '-generate') {
    my ($index);
    if (defined[$ARGV[1]]) {
	# User passed in a parameter to specify which color scheme to pick
	$index = $ARGV[1];
    } else {
	$index = int(rand(scalar(@profiles)));
    }
    my $profile = $profiles[$index];
    ($cursor_color, $text_color, $bold_color, $bg_color) =
	map('{ ' . rgb2string(hsv2rgb(@$_)) . ' }', @$profile);

} elsif ($ARGV[0] eq '-cycle') { # cycle through the colors
    my $prefFile = $ENV{'HOME'} . '/' . '.rand_term_color_setting'; # May generate a preference dot file. Ugh! Pollutes your home directory.
    my $termCount = 0;
    #print "Does $prefFile exist? We think " . (-e $prefFile) . "\n";
    if (-e $prefFile) {
	my $readFromFile = `head -n 1 $prefFile`;
	chomp($readFromFile);
	if (length($readFromFile) > 10) {
	    print "Removing the malformed preference file (it\'s too long)...\n";
	    system("rm -f $prefFile");
	} else {
	    $termCount = $readFromFile;
	}
	if ($termCount >= scalar(@cycle_profiles)) {
	    $termCount = 0;
	}
    }
    open PREFFILE, "> $prefFile" or die $!;
    print PREFFILE ($termCount+1) . "\n";
    close PREFFILE;
    #system("echo " . ($termCount+1) . " > " $prefFile)
    
    my $profile = $cycle_profiles[$termCount];
    ($cursor_color, $text_color, $bold_color, $bg_color) = map('{ ' . rgb2string(hsv2rgb(@$_)) . ' }', @$profile);
    
} else {
    die "Unknown options";
}

open(AS, "|osascript");
print AS "tell application \"Terminal\"
	set the  background  color of window 1 to $bg_color
	set the normal text  color of window 1 to $text_color
	set the   bold text  color of window 1 to $bold_color
	set the      cursor  color of window 1 to $cursor_color
end tell";
close(AS)

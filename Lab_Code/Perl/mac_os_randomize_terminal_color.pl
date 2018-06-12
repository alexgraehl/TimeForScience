#!/usr/bin/env perl
use strict;

exit(0) unless (exists($ENV{'TERM_PROGRAM'}) 
		&& ($ENV{'TERM_PROGRAM'} eq "Apple_Terminal"));

### 
### Alex Williams' Terminal.app color randomizer, with minor additions by Charlie Vaske
### 

# parameters which determine the range of random colors chosen
my  ($min_h,$max_h) =  ( 0,    360   );   # hue        (0-360) (degrees of angle on color wheel)
my  ($min_s,$max_s) =  ( 0.90, 0.91  );   # saturation (0-1) ("how washed out" things are)
my  ($min_v,$max_v) =  ( 0.14, 0.15  );   # brightness (0-1)
my  ($min_a,$max_a) =  ( 1.0,  1.0  );    # alpha      (0-1) (1.0 = opaque, to 0.0 = transparent)


my @profiles = ([ [0,0,1/3], [0,0,0], [0,0,0], [0,0,1] ], # black on white
		[ [0,0,1],   [0,0,1], [0,0,1], [0,0,0] ], # white on black
		[ [120,1,1], [120,1,1], [120,1,1], [0,0,0] ], # gren on blck
		[ [0,0,1/3], [0,0,0], [0,0,0], [60,0.29,1] ]  # blck on ylw
    );

# generate random number between two values
sub randbetween {
    my ($min,$max) = @_;
    return $min + rand($max-$min);
}



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

# turn list of floats into a string suitable for applescript
sub rgb2string {
    return join ", ", map(int(65535*$_),@_);
}

#main program

my ($bg_color, $text_color, $bold_color, $cursor_color);
if (scalar(@ARGV) == 0) {
    my $randHue = randbetween($min_h, $max_h);
    my @rgb = hsv2rgb($randHue, randbetween($min_s,$max_s), randbetween($min_v,$max_v));
    my $alpha = int(65535); #int(65535*randbetween($min_a, $max_a));
    
    $bg_color   = '{ ' . rgb2string(@rgb). ", $alpha }";
    
    @rgb = hsv2rgb($randHue, 0.20, 1.0); # Hue Saturation Value(lightness)
    $text_color = '{ ' . rgb2string(@rgb). ' }';
    
    @rgb = hsv2rgb($randHue, 0.50, 1.0); # Hue Saturation Value(lightness)
    $bold_color = '{ ' . rgb2string(@rgb). ' }';
    
    @rgb = hsv2rgb(180 - $randHue, 0.70, 1.0); # Hue Saturation Value(lightness)
    $cursor_color = '{ ' . rgb2string(@rgb). ' }';
    
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
} else {
    die "Unknown options";
}

open(AS, "|osascript");
print AS "tell application \"Terminal\"
	set the background color of window 1 to $bg_color
	set the normal text color of window 1 to $text_color
	set the bold text color of window 1 to $bold_color
	set the cursor color of window 1 to $cursor_color
end tell";
close(AS)

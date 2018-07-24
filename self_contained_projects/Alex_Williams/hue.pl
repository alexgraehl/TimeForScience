#!/usr/bin/perl -w

# The 'curl' part of the code (the command that controls the lights over the network) is from:
#   * http://www.everyhue.com/vanilla/discussion/136/perl-script-to-control-hue-light

use strict;
use warnings;

use Getopt::Long;
use Term::ANSIColor;
use Time::HiRes; # Allows usleep for sleeping at sub-second intervals
use List::Util qw(max min);

my $baseIP     = $ENV{"HUE_BASE_IP"};
if (not $baseIP) { die "You need to set the HUE_BASE_IP address in your .bashrc. For example, try adding code for the correct address (which you will need to 1) find in your router configuration and 2) make sure it is a STATIC IP ALLOCATION! Then add it like this to your ~/.bashrc:     export HUE_BASE_IP=10.1.1.99   or  export HUE_BASEIP=192.168.1.111   , or whatever you have picked for your IP address."; }

my $deviceType = "Talk to Hue";
my $user       = "61"x16; # apparently needs to be 16 digits represented in hex. "a" = 61, "b" = 62 ... "o" = 6f ... etc.

my $DEBUG = 0; # Is debugging on? Normally it's off (0).

my $pingExitCode = system(qq{ping -c 1 $baseIP > /dev/null}); # $? is a perl special variable
if ($pingExitCode != 0) {
    print STDERR colored("FAILED to ping the base station at \"$baseIP\" (exit code was non-zero).\nNote that the base station IP address will most likely start with either 10.0 or 192.168, depending on how your router works.\nCheck to make sure your base IP address is correct! Ping exit code was \"$pingExitCode\".", "white on_red");
    print STDERR "\n"; # <-- prevents the red background from continuing to the next line
}


sub quitWithUsageError($) { print($_[0] . "\n"); printUsage(); print($_[0] . "\n"); }
sub printUsageAndQuit() { printUsage(); exit(1); }
sub printUsage() {    print STDOUT <DATA>; exit(0); }

sub sleepFloat($) { # input: a decimal amount of seconds to sleep!
    Time::HiRes::usleep(int($_[0] * 1_000_000)); # <-- one million microseconds are in one second
}



my $MAX_LOGIN_ATTEMPTS = ($DEBUG) ? 1 : 10; # how many times to try to connect to the light before skipping it. Default: 10 attempts, unless we are debugging, in which case set it to 1 just to make things easier to check.
my $HUE_MAX = 65535; # Values of HUE can be [0-65535]
my $SAT_MAX = 255; # [0-255]
my $VAL_MAX = 255; # [0-255]

my $MIN_SECONDS_BETWEEN_UPDATES = 0; #(1.0/30.0);

my $lightIDStringCommaDelim = undef;
my $hue360 = undef;   # hue on a 0 to 360 scale. The Philips Hue light system actually takes a value from 0 to 65535.
my $sat100 = undef;
my $val100 = undef;

my $hue65535 = undef; # This is the number that the hue *actually* accepts. Defualt hue is 0, which is red.
my $sat255   = undef;
my $val255   = undef;


my $transitionInSeconds = undef;
my $off = 0;

my $doFlicker = 0;
my $doRainbow = 0;
my $doCarnival = 0;
my $doHorrid = 0;
my $doFireplace = 0;

GetOptions("help|?|man"       => sub { printUsageAndQuit(); }
	   , "n=s"                => \$lightIDStringCommaDelim
	   , "off"                => sub { $off = 1; }
	   , "hue|h=f"            => \$hue360   # accepts a floating point value for the hue, from 0.0 to 360.0
	   , "sat|saturation|s=f" => \$sat100 # 0 to 255
	   , "val|value|v=f"      => \$val100 # 0 to 255

	   , "longhue|hh=i"       => \$hue65535 # accepts an INTEGER for the value, from 0 to 65535.
	   , "longsat|ss=i"       => \$sat255 # 0 to 255. Accepts an INTEGER only.
	   , "longval|vv=i"       => \$val255 # 0 to 255. Accepts an INTEGER only.

	   , "time|t=f"           => \$transitionInSeconds
	   , "rainbow!"            => \$doRainbow
	   , "carnival!"           => \$doCarnival
	   , "horrid!"             => \$doHorrid

	   , "flicker!"             => \$doFlicker

	   , "fire|fireplace!"             => \$doFireplace

	   , "fast"    => sub { $transitionInSeconds = 0.15; }

	   , "red"         => sub { $hue360 = 0; $sat100 = 100; $val100 = 100; }
	   , "grapefruit|gf"  => sub { $hue360 = 30; $sat100 = 100; $val100 = 100; }
	   , "orange"  => sub { $hue360 = 50; $sat100 = 100; $val100 = 100; }
	   , "amber"   => sub { $hue360 = 90; $sat100 = 100; $val100 = 100; }
	   , "yellow"  => sub { $hue360 = 120; $sat100 = 100; $val100 = 100; }
	   , "green"   => sub { $hue360 = 145; $sat100 = 100; $val100 = 100; }
	   , "blue"    => sub { $hue360 = 255; $sat100 = 100; $val100 = 100; }
	   , "purple"  => sub { $hue360 = 260; $sat100 = 100; $val100 = 100; }
	   , "violet"  => sub { $hue360 = 270; $sat100 = 100; $val100 = 100; }
	   , "magenta" => sub { $hue360 = 310; $sat100 = 100; $val100 = 100; }
	   , "pink"    => sub { $hue360 = 320; $sat100 = 80; $val100 = 100; }
	   , "white"   => sub { $hue360 = 0; $sat100 = 0; $val100 = 100; }
	   , "gray|grey"    => sub { $hue360 = 0; $sat100 = 0; $val100 = 50; }

    ) or printUsageAndQuit();

my @lights = ();

if (!defined($lightIDStringCommaDelim)) {
    @lights = (1,2,3); # all lights, assuming there are three <------------- HARD CODED HERE!
} else {
    @lights = split(",", $lightIDStringCommaDelim);
}

## Figure out the hue in terms that the base station understands!
(not(defined($hue65535) && defined($hue360))) or die "You can specify EITHER --hue (-h) OR --longhue (--hh), but not both at the same time!\n";
(not(defined($sat100) && defined($sat255)))   or die "You can specify EITHER --sat (-s) OR --longsat (--ss), but not both at the same time!\n";
(not(defined($val100) && defined($val255)))   or die "You can specify EITHER --val (-v) OR --longval (--vv), but not both at the same time!\n";

if (defined($hue360)) {    $hue65535 = int($hue360 / 360.0 * $HUE_MAX); } # Convert 0-360 to 0-65535
if (defined($sat100)) {    $sat255   = int($sat100 / 100.0 * $SAT_MAX); } # Convert 0-100 to 0-255
if (defined($val100)) {    $val255   = int($val100 / 100.0 * $VAL_MAX); } # Convert 0-100 to 0-255

(!defined($hue65535) or ($hue65535 >= 0 && $hue65535 <= $HUE_MAX)) or die "The invalid hue value \"$hue65535\" was out of range (Valid range is from 0 to $HUE_MAX).";
(!defined($sat255) or ($sat255 >= 0 && $sat255 <= $SAT_MAX))     or die "The invalid saturation value \"$sat255\" was out of range (Valid range is from 0 to $SAT_MAX).";
(!defined($val255) or ($val255 >= 0 && $val255 <= $VAL_MAX))     or die "The invalid value (brightness) value \"$val255\" was out of range (Valid range is from 0 to $VAL_MAX).";

sub logIntoBaseStation($$$) {
    my ($baseIP, $deviceType, $user) = @_;
    for (my $count = 0; $count < $MAX_LOGIN_ATTEMPTS; $count++) {
	print STDOUT "[$count/$MAX_LOGIN_ATTEMPTS]: Attempting to connect... please press the link button on the Hue base station.\n";
	my $curlResult = `curl -4 -s -XPOST http://$baseIP/api/ -d '{\"devicetype\":\"$deviceType\", \"username\":\"$user\" }'`;
	last if ($curlResult =~ /success/); # break out of the loop if we actually succeeded!
	sleepFloat(1.0); # Sleep for one second before looping!
    }
    print colored("Successful authentication!\n", 'green');
}

sub setLights(     $$$$$$$@) {
    my($baseIP, $deviceType, $user, $hue, $sat, $brightness, $ttime, @lightNumArray) = @_;
    for my $lnum (@lightNumArray) {
	setSingleLight($baseIP, $deviceType, $user, $hue, $sat, $brightness, $ttime, $lnum);
    }
}

sub setSingleLight($$$$$$$$) {
    my($baseIP, $deviceType, $user, $hue, $sat, $brightness, $ttime, $lightNum) = @_;

    my $hueStr = (defined($hue)) ? qq{\"hue\":$hue,} : qq{};
    my $satStr = (defined($sat)) ? qq{\"sat\":$sat,} : qq{};
    my $briStr = (defined($brightness)) ? qq{\"bri\":$brightness,} : qq{};

    if (!defined($ttime)) { $ttime = 1.0; }
    my $timeInTenthsOfSecondsString = (defined($ttime)) ? (", \"transitiontime\":" . int($ttime * 10)) : '';

    if ($DEBUG) {
	print "User: $user\n";
	print "hue: $hueStr\n";
	print "sat: $satStr\n";
	print "bri: $briStr\n";
	print "ttime-in-tenths: $timeInTenthsOfSecondsString\n";
    }

    my $isOnBool = (defined($brightness) && $brightness <= 0) ? "false" : "true"; # The light is only off if brightness is both defined AND ALSO <= 0.

    for (my $i = 0; $i < $MAX_LOGIN_ATTEMPTS; $i++) {
	#print "Attempting to contact $baseIP...";
	my $cmd = qq{curl -s --request PUT --data '{$hueStr $satStr $briStr \"on\":$isOnBool $timeInTenthsOfSecondsString}' http://$baseIP/api/$user/lights/$lightNum/state/};
	my $curlResult = `$cmd`;
	($DEBUG) && print STDERR colored("COMMAND: $cmd\n", "white on_blue"); # Print the error text
	last if ($curlResult !~ /error/); # if it worked, then break out of the loop.
	($DEBUG) && print STDERR colored("FAILED COMMAND: $cmd\n", "white on_red"); # Print the error text
	($DEBUG) && print STDERR colored("RESULT: $curlResult\n", "white on_red"); # Print the error text

	if ($curlResult =~ /unauthorized user/) {
	    logIntoBaseStation($baseIP, $deviceType, $user); 
	}
	sleepFloat(0.1); # Sleep for 1/10 of a second before trying to log in again.
    }
}

sub setLightsOff($$$$@) {
    my($baseIP, $deviceType, $user, $ttime, @lightNumArray)=@_;
    setLights($baseIP, $deviceType, $user, undef, undef, 0, 2, @lightNumArray);
}

# sanity check
#if ( $LightNum < 1 or $LightNum > 8 ) { &usage( "LightNum ($LightNum) must be [1..8]" ); }
#if ( $Hue < 0 or $Hue > 65535 ) { &usage( "hue ($Hue) must be 0-65535" ); }
#if ( $Saturation < 0 or $Saturation > 255 ) { &usage( "Saturation ($Saturation)must be [0..255]" ); }
#if ( $Brightness < 0 or $Brightness > 255 ) { &usage( "Brightness ($Brightness) must be [0..255]" ); }
#if ( $TransitionTime < 0 or $TransitionTime > 65535 ) { &usage( "Transitiontime ($TransitionTime) must be [0..65535] (1/10 seconds)"); }

# set the HueLight 

if ($off) {
    setLightsOff($baseIP, $deviceType, $user, 10, @lights);
    exit(0); # <-- important!
}


if ($doRainbow) {
    print "RAINBOW\n";
    my $rainbowHue = 0;
    my $rainbowTTime = defined($transitionInSeconds) ? $transitionInSeconds : 0.2; # default carnival transition time is 0.1 seconds, unless this is specifically overridden
    while (1) {
	# INFINTE LOOP!
	$rainbowHue += 1000;
	if ($rainbowHue > $HUE_MAX) { $rainbowHue = 0; }
	setLights($baseIP, $deviceType, $user, $rainbowHue, $sat255, $val255, $rainbowTTime, @lights);
	my $sleepTime = max($rainbowTTime, $MIN_SECONDS_BETWEEN_UPDATES);
	sleepFloat($sleepTime+0.01);
    }
}

if ($doCarnival) {
    print "CARNIVAL\n";
    my $hhh = 0;
    my $carnivalTransitionTime = defined($transitionInSeconds) ? $transitionInSeconds : 0.1; # default carnival transition time is 0.1 seconds, unless this is specifically overridden
    while (1) {
	# INFINTE LOOP!
	$hhh += int(($HUE_MAX / 8) + rand($HUE_MAX / 2)); # randomly(ish) change the hue, but ensure that it changes by at least (max/8) 
	if ($hhh > $HUE_MAX) { $hhh = min($hhh - $HUE_MAX, $HUE_MAX); }
	setLights($baseIP, $deviceType, $user, int($hhh), $sat255, $val255, $carnivalTransitionTime, @lights);
	my $sleepTime = max($carnivalTransitionTime, $MIN_SECONDS_BETWEEN_UPDATES);
	sleepFloat($sleepTime+0.01);
	#print $hhh . "\t";
    }
}

if ($doHorrid) {
    print "HORRID\n";
    my $hhh = 0;
    my $RED    = int(  0.0/360.0 * $HUE_MAX);
    my $BLUE   = int(255.0/360.0 * $HUE_MAX);
    my $horridTransitionTime = defined($transitionInSeconds) ? $transitionInSeconds : 0.1;
    for (my $i = 0; ; $i++) {
	# INFINTE LOOP!
	$hhh = ($i % 2 == 0) ? $RED : $BLUE;
	setLights($baseIP, $deviceType, $user, int($hhh), $sat255, $val255, $horridTransitionTime, @lights);
	my $sleepTime = max($horridTransitionTime, $MIN_SECONDS_BETWEEN_UPDATES);
	sleepFloat($sleepTime);
    }
}

if ($doFlicker) {
    print "FLICKER\n"; # like a flickering CRT monitor
    for (my $i = 0; ; $i++) { # <-------------- INFINITE LOOP
	my $sceneTransitionSpeed = 0.15 + rand(0.6); # The "flicker" speed is between 0.1 and 0.6 seconds
	my $sceneLength = (rand(8.0) + 0.5) * ($sceneTransitionSpeed*$sceneTransitionSpeed); # Each scene is at least 1 second, but could be shorter
	my $lengthSoFar = 0;
	for (my $lengthSoFar = 0; $lengthSoFar <= $sceneLength; ) {
	    
	    my $fh = int((210+rand(40))/360*$HUE_MAX); # acceptable hue range on the 360 scale is from 180 to 260 (blues)

	    if (rand(1.0) < 0.02) { $fh = 0; } # 2% chance of setting the hue to red!
	    elsif (rand(1.0) < 0.02) { $fh = 120; } # ~2% chance of setting the hue to yellow!

	    my $fs = int((rand(0.5)+0.5)*$SAT_MAX); # saturation can vary from 0.5 to 1.0.
	    my $fv = int((rand(0.4)+0.6)*$VAL_MAX); # value can vary from 0.4 to 1.0
	    setLights($baseIP, $deviceType, $user, $fh, $fs, $fv, $sceneTransitionSpeed, @lights);
	    my $sleepTime = max($sceneTransitionSpeed, $MIN_SECONDS_BETWEEN_UPDATES);
	    sleepFloat($sleepTime);
	    $lengthSoFar += $sleepTime;
	}
    }
}


if ($doFireplace) {
    print "FIREPLACE\n"; # like a crackling fire
    for (my $i = 0; ; $i++) { # <-------------- INFINITE LOOP
	my $sceneTransitionSpeed = 0.5 + rand(0.4); # The "flickering" speed
	my $sceneLength = 2 + rand(4.0); # Each scene length
	my $lengthSoFar = 0;
	for (my $lengthSoFar = 0; $lengthSoFar <= $sceneLength; ) {
	    my $fh = int((40+rand(60))/360*$HUE_MAX); # acceptable hue range.
	    my $fs = undef; #int((rand(0.4)+0.4)*$SAT_MAX); # saturation can vary from 0.4 to 0.8.
	    my $fv = undef; #int((rand(0.15)+0.85)*$VAL_MAX); # value can vary from 0.85 to 1.0
	    setLights($baseIP, $deviceType, $user, $fh, $fs, $fv, $sceneTransitionSpeed, @lights);
	    my $sleepTime = max($sceneTransitionSpeed, $MIN_SECONDS_BETWEEN_UPDATES);
	    sleepFloat($sleepTime);
	    $lengthSoFar += $sleepTime;
	}
    }
}

setLights($baseIP, $deviceType, $user, $hue65535, $sat255, $val255, $transitionInSeconds, @lights);

__DATA__

Hue Setter.

Example usage:
   hue.pl --red -n 1,3
  Set lights 1 and 3  to red.

  hue.pl --rainbow
  Set the rainbow color-changing mode.
  Requires that this program continue running!

Other special modes. These all require that the program continue running:

--rainbow: rainbow
--fireplace: red/yellow slow changes
--flicker: Like the flickering of a monitor or TV
--horrid: red/blue "emergency light" flashing
  
-n (comma delimited numeric list): which lights to change
  e.g.     -n 1,3,4 (no spaces)

  --off: Turn lights off

  
  
  --hue or -h (0 to 360)  (0 = red, 90 = yellow, etc...)
  --sat or -s (0 to 100)  (saturation -- 0 = white, 100 = saturated)
  --val or -s (0 to 100)  (brightness: 0 = darkes)

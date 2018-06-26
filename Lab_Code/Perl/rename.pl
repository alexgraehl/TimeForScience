#!/usr/bin/env perl
#
#  This script was developed by Robin Barker (Robin.Barker@npl.co.uk),
#  from Larry Wall's original script eg/rename from the perl source.
#
#  This script is free software; you can redistribute it and/or modify it
#  under the same terms as Perl itself.
#
# Larry(?)'s RCS header:
#  RCSfile: rename,v   Revision: 4.1   Date: 92/08/07 17:20:30
#
# $RCSfile: rename,v $$Revision: 1.5 $$Date: 1998/12/18 16:16:31 $
#
# $Log: rename,v $
# Revision 1.6: 2018: Color support, certain warnings, updates.
# Revision 1.5  1998/12/18 16:16:31  rmb1
# moved to perl/source
# changed man documentation to POD
#
# Revision 1.4  1997/02/27  17:19:26  rmb1
# corrected usage string
#
# Revision 1.3  1997/02/27  16:39:07  rmb1
# added -v
#
# Revision 1.2  1997/02/27  16:15:40  rmb1
# *** empty log message ***
#
# Revision 1.1  1997/02/27  15:48:51  rmb1
# Initial revision
#

use strict;  use warnings;  use diagnostics;
use List::Util qw(max min);
use Getopt::Long;
use Carp; # backtrace on errors. Has the "confess" function. Use this instead of "die" if you want useful information!

$| = 1; # Always flush text output IMMEDIATELY to the console, don't wait to buffer terminal output! Setting this to zero can cause STDERR and STDOUT to be interleaved in weird ways.

#no warnings 'numeric';
#use Scalar::Util;
#print Scalar::Util::looks_like_number($string), "\n";

sub tryToLoadModule($) {
	my $x = eval("require $_[0]");
	if ((defined($@) && $@)) {
		warn "We FAILED to load module $_[0]. Skipping this module, but continuing with the program.";
		return 0;	# FAILURE
	} else {
		$_[0]->import();
		return 1;	# SUCCESS
	}
}

my $SHOULD_USE_COLORS = tryToLoadModule("Term::ANSIColor");
if ($SHOULD_USE_COLORS) {
	use Term::ANSIColor;
}

sub warnPrint($) { chomp($_[0]); warn(safeColor("[WARNING]: " . $_[0] . "", "yellow on_black")); } # regarding "warn": if the string ends with a newline it WON'T print the line number!

sub safeColor($;$) {		# one required and one optional argument
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
	if (! -t STDERR) {
		$col = undef;
	}			# no coloration if this isn't to a terminal
	print STDERR safeColor($msg, $col);
}

sub printColorStdout($;$) {
	# prints color to STDOUT *UNLESS* it is re-directed to a file, in which case NO COLOR IS PRINTED.
	my ($msg, $col) = @_; # Only prints in color if STDOUT is to a terminal, NOT if it is redirected to an output file!
	if (! -t STDOUT) {
		$col = undef;
	}			# no coloration if this isn't to a terminal
	print STDOUT safeColor($msg, $col);
}


sub printOkStdout($) {
	my ($msg) = @_;
	printColorStdout($msg, "green"); # OR: green on_black
}

sub printNonFatalFailure($) {
	my ($msg) = @_;
	printColorStderr($msg, "red"); # OR: red on_black
	warn $msg;
}

sub dryNotify(;$) {		# one optional argument
	my ($msg) = @_;
	$msg = (defined($msg)) ? $msg : "This was only a dry run, so we skipped executing a command.";
	print STDERR safeColor("[DRY RUN]: $msg\n", "black on_yellow");
}

sub notify($) {			# one required argument
	my ($msg) = @_;
	warn safeColor("[DRY RUN]: $msg\n", "cyan on_blue");
}

sub main();
sub quitWithUsageError($) { print($_[0] . "\n"); printUsageAndQuit(); print($_[0] . "\n"); }
sub printUsageAndQuit() { printUsage(); exit(1); }

sub printUsage() {
	print STDOUT "\"rename\" requires input arguments.\nUsage: rename [-v] [-n] [-f] perl_regular_expression [filenames]\nExample: rename 's/TXT/new_processed_text/' *.TXT\n\n";
	print STDOUT <DATA>;
	exit(0);
}




my ($quiet, $no_act, $force, $renameRegexp);
Getopt::Long::Configure('bundling');
$Getopt::Long::passthrough = 1; # ignore arguments we don't recognize in GetOptions, and put them in @ARGV

GetOptions('q|quiet'   => \$quiet
	   ,'n|no-act|dry|dry-run|dryrun'  => \$no_act
	   ,'f|force'   => \$force
	   ,'help|?|h'  => sub { printUsage(); exit(0); }
	  ) or printUsageAndQuit();

($renameRegexp = shift) or printUsageAndQuit();
$renameRegexp =~ m/^s/ or die "[rename.pl] INCORRECT USAGE: The regular expression MUST start with 's'. Example: 's/8/9/g' to change all '8's to '9's.\n";

#$renameRegexp =~ m/^s\/(.*)\/([a-zA-Z]*)$/i;
#my $rexMain = $1; # The main body of the regular expression. Not distinguishing the "before" and "after" parts yet.
#my $rexPost = $2; # The letters at the end that modify the regexp. Like 'i' for case-insensitive
#print("The input regexp was:  $renameRegexp, which is split into $rexMain and $rexPost\n");

#printColorStdout("hey", "red");
#printColorStdout("hey", "red");
#printColorStdout("what", "red");
#printColorStderr("stderr what", "blue on_white");
#printColorStdout("is", "red on_blue");
#printColorStdout("this", "red on_green");
#printColorStderr("stderr what", "red on_white");
#print STDERR "hellow\n";
#print STDERR "Test color!\n";

if (!@ARGV) {
	print "Reading filenames from STDIN...\n" if !$quiet;
	@ARGV = <STDIN>;
	chop(@ARGV);
}

if (!$quiet) {
	print "About to examine " . scalar(@ARGV) . " files/folders that might be renamed if they match the input regexp pattern of " . safeColor($renameRegexp, "white on_blue") . " . . .\n";
}

my $longestFilenameLen = 0;
for my $filename (@ARGV) {
	if (length($filename) > $longestFilenameLen) {
		$longestFilenameLen = length($filename);
	}
}

for (@ARGV) {
	my $fromName = $_;
	#print $rexMain . " is the main...\n";
	#if ($fromName =~ m/$rexMain/) {
	#	print "match found\n";
	#}
	eval $renameRegexp;
	my $toName = $_;

	# Dubious way of figuring out which areas of the string have changed
	#my @clef = (); # changes LEFT edge
	#my @crig = (); # changes RIGHT edge
	#my $minLen = min(length($fromName), length($toName));
	#my $goodLeft = 0;
	#my $goodRight = $maxLen;
	#while (my $i = 0; $i < $maxLen; $i++) {
	#	if ($fromName[$i] ne $toName[i]) {
	#		$goodThroughLeft++;
	#	}
	#}
	#if ($fromName =~ $renameRegexp) {
	#	print "yep\n";
	#}
	
	#my $evalResult = eval $toName =~ $renameRegexp;     ## apply the renaming operation to "$_" by default --- this changes "$_" to the newly-renamed version, but does NOT rename the file on disk yet!
	die $@ if $@;		## ?????????????? why is perl so confusing

	my $numSpacesToPad = $longestFilenameLen - length($fromName);
	my $pad            = (' ' x $numSpacesToPad);

	if ($fromName eq $toName) {
		if (!$quiet) {
			print STDOUT safeColor("$fromName ", "blue")
			  . $pad . safeColor(" (Not renamed--no changes)", "blue") . "\n";
		}
		next;
	}

	if ($no_act) {
		if (!$quiet) {
			print STDOUT safeColor("Dry run: $fromName ", "green")
			  . $pad . safeColor("---DRY-RUN--->", "cyan")
			  . " " . $toName . "\n";
		}
		next; # NEXT!!!!! Skip the rename code for this file...
	}
	
	if (-e $toName and !$force) {
		printNonFatalFailure("$fromName was not renamed: $toName already exists");
		next; # NEXT!!!!! Skip renaming...
	}

	if (rename $fromName, $toName) { ## <-- now try to actually rename the file on disk!
		if (!$quiet) {
			print STDOUT safeColor("$fromName ", "green")
			  . $pad . safeColor("--RENAMED-->", "yellow")
			  . " " . $toName . "\n";
		}
		#print "$fromName  renamed to  $toName\n";
	} else {
		warn "Can't rename $fromName $toName: $!\n"; ## any error codes are stored in "$!" automatically
	}
}


END {
    # Runs after everything else.
    # Makes sure that the terminal text is back to its normal color.
    if ($SHOULD_USE_COLORS) { print STDERR Term::ANSIColor::color("reset"); print STDOUT Term::ANSIColor::color("reset"); }
}

exit(0);

__DATA__

=head1 NAME

rename - renames multiple files

=head1 SYNOPSIS

B<rename> S<[ B<-v> ]> S<[ B<-n> ]> S<[ B<-f> ]> I<perlexpr> S<[ I<files> ]>

=head1 DESCRIPTION

C<rename>
renames the filenames supplied according to the rule specified as the
first argument.
The I<perlexpr> 
argument is a Perl expression which is expected to modify the C<$_>
string in Perl for at least some of the filenames specified.
If a given filename is not modified by the expression, it will not be
renamed.
If no filenames are given on the command line, filenames will be read
via standard input.

For example, to rename all files matching C<*.bak> to strip the extension,
you might say

	rename 's/\.bak$//' *.bak

To translate uppercase names to lower, you'd use

	rename 'y/A-Z/a-z/' *

Backreference example:

rename.pl -n 's/old_thing(.txt|.gz)/new_thing$1/' *.txt *.gz
(You generally should not need to use backreferences, but if you do, here they are.)


=head1 OPTIONS

=over 8

=item B<-v>, B<--verbose>

Verbose: print names of files successfully renamed.

=item B<-n>, B<--no-act>

No Action: show what files would have been renamed.

=item B<-f>, B<--force>

Force: overwrite existing files.

=back

=head1 ENVIRONMENT

No environment variables are used.

=head1 AUTHOR

Larry Wall

=head1 SEE ALSO

mv(1), perl(1)

=head1 DIAGNOSTICS

If you give an invalid Perl expression you'll get a syntax error.

=head1 BUGS

The original C<rename> did not check for the existence of target filenames,
so had to be used with care.  I hope I've fixed that (Robin Barker).

Note: are you trying to rename files in a directory to match a LIST of new names?
      Although that is not part of this script, here is a shell one-liner that will do it:
      (Requires a 1-column file named "NEW_FILENAME_LIST.txt" of new filesnames, in the same
      order that "ls -1" uses. Probably breaks on special characters, so beware.)
      while read -u 9 SRC; do read -u 8 DEST; echo mv "$SRC" "$DEST"; done ; 9<<<"$(ls -1 * | grep -v NEW_FILENAME_LIST.txt)" 8< NEW_FILENAME_LIST.txt && echo "[DONE]"
      Note that the "9<<<" is in fact a real thing, not a typo! It matches the 'read -u 9' earlier.

=cut

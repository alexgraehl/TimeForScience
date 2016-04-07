#!/usr/bin/perl

#@COMMENT@ sort.pl can handle COMPRESSED (gzip/bzip2) files and can accept header line(s). It uses the fast UNIX sort internally. Frequency-of-use rating: 9/10.

# This is basically just a wrapper to GNU "sort" that adds the option of having a header line.
# It takes ONE filename (the very last argument) and any number of arguments to pass along to UNIX sort.
#use POSIX      qw(ceil floor); #use List::Util qw(max min); #use File::Basename;

use File::Temp;
use Term::ANSIColor;
use Getopt::Long;
use strict; use warnings; use diagnostics;

sub printUsage() {  print STDOUT <DATA>; }
sub printUsageAndQuit() {  printUsage();  exit(0); } # Exit with an "OK" code
sub quitWithUsageError($) {
	print color("red");  print "Error: " . $_[0];  print color("reset");
	print "\n";
	printUsage();
	print color("red");  print "Error: " . $_[0];  print color("reset");
	print "\n";
	exit(1); # Exit with ERROR code
}

sub agw_make_temp() { my ($temp_filehandle, $temp_filename) = File::Temp::tempfile(DIR=>"./", SUFFIX=>".sort.tmp"); return($temp_filename); }

sub main();

# ==1==
sub main() { # Main program
	my $delim = "\t";
	my $numHeaderLines = 0;
	$Getopt::Long::passthrough = 1; # ignore arguments we don't recognize in GetOptions, and put them in @ARGV
	GetOptions("help|?|man" => sub { printUsageAndQuit(); }
		   , "header|h=s" => \$numHeaderLines
		   , "d=s" => \$delim
		  ) or printUsageAndQuit();
	#quitWithUsageError("1 == 0? Something is wrong!");

	my $numUnprocessedArgs = scalar(@ARGV);
	if ($numUnprocessedArgs == 0) {
		# I guess we're reading from STDIN.
		#quitWithUsageError("Error in arguments! You must send exactly one filename to this program. Note that currently you CANNOT pipe data to it through <STDIN>. That should be fixed.\n");
	}

	my $numValidFilesWeFound = 0;
	foreach my $arg (@ARGV) { # Go through the unprocessed args and look for an erroneous header argument
		if (-f $numValidFilesWeFound) { $numValidFilesWeFound++; }
		if ($arg =~ m/^[-]h\d+$/) {
			quitWithUsageError("Hey, you shouldn't use the syntax -h1 or -h2, use -h 1 or -h=1! This isn't 'head' or 'tail', where that kind of weird thing is allowed! And if you want to use the '-h' option in sort, it doesn't take any arguments!");
		}
	}

	if ($numValidFilesWeFound >= 2) {
		print color("yellow");  print "WARNING: It looks like you passed more than one file to sort.pl, which is NOT supported! You should double check this...\n";  print color("reset");
		#($numUnprocessedArgs == 1) or quitWithUsageError("ERROR: You apparently passed in MORE THAN ONE filename to sort.pl! Just pass in one!");
	}

	my $sortArgsString; # additional arguments to sort! Note that we still only take one filename

	my $STDIN_FILENAME = agw_make_temp();

	# Oh huh, we're reading from STDIN! We'll have to output a TEMP file, because we're going to be using the first line of this file twice.
	my $filename = (scalar(@ARGV) > 0) ? $ARGV[-1] : undef; # Last argument is the file. Previous ones could be arguments to 'sort'
	if (!defined($filename) or (not -f $filename)) {
		# Looks like the user did NOT actually pass in a file after all.
		# They must have passed in some more arguments to sort, I gues...
		$filename = $STDIN_FILENAME;
		$sortArgsString = join(" ", @ARGV);
		#if ($numHeaderLines > 0) {
		# I guess we need to have a temp file if we have header lines AND are reading from STDIN
		open (OUT, "> ${STDIN_FILENAME}"); # We are reading from stdin...
		while (<STDIN>) { print OUT $_; }
		close(OUT);
		#}
	} else {
		$sortArgsString = join(" ", @ARGV[0..($#ARGV-1)] ); # the remaining arguments on the command line!
	}

	my $HEADER_FILENAME = agw_make_temp();
	
	print STDERR color("yellow");
	print STDERR "STDERR: Sorting begins here: ===============================\n";
	print STDERR color("reset");

	my $sortCommand = "sort -t \'$delim\' $sortArgsString";

	# =============== TRANSPARENTLY READ A GZIP / BZIP2 / ZIP / UNCOMPRESSED FILE =======
	my $catter; # <-- the thing that does the 'cat' command
	#if (!defined($filename) or $filename eq $STDIN_FILENAME) { $catter = ''; }
	if    ($filename =~ /[.](gz|gzip)$/i)       { $catter = "gzip --decompress --stdout"; }     # Un-gzip a file and send it to STDOUT.
	elsif ($filename =~ /[.](bz2|bzip2)$/i)     { $catter = "bzip2 --decompress --stdout"; }    # Un-bz2 a file and send it to STDOUT
	elsif ($filename =~ /[.](zip)$/i)           { $catter = "unzip -p"; } # Un-regular-zip a file and send it to STDOUT with "-p": which is DIFFERENT from -c (-c is NOT what you want here). See 'man unzip'
	else                                        { $catter = "cat"; }  # Default: just read a file normally

	if ($filename =~ /[.](tar|rar|7z|Z)$/i) { quitWithUsageError("Unfortunately, we handle only gzip, bzip2, and zip files. We can't currently handle tar files or .Z ('compress') files, or rar files."); }

	if ($numHeaderLines > 0) {
		my $nHeaderPlus1 = $numHeaderLines+1;
		system(qq{$catter "$filename" | head -n$numHeaderLines > "$HEADER_FILENAME"});
		system(qq{$catter "$filename" | tail -n +$nHeaderPlus1 | $sortCommand | cat "$HEADER_FILENAME" - });  # <-- this generates (to STDOUT) the final output that the user sees
		unlink($HEADER_FILENAME); # delete the temp files
		unlink($STDIN_FILENAME); # delete the temp files
	} else {
		system(qq{$catter "$filename" | $sortCommand});
	}
	print STDERR color("green");
	print STDERR "STDERR: Sorting ends here: ===============================\n";
	my $headerStr = "";
	if ($numHeaderLines == 1) { $headerStr = "(with one header line) "; }
	if ($numHeaderLines > 1) { $headerStr = "(with $numHeaderLines header lines) "; }
	print STDERR "STDERR: Finished sorting \"$filename\" " . $headerStr . "! Sort command was: $sortCommand\n";
	print STDERR color("reset");
} # end main()


main();


END {
  # Runs after everything else.
  # Makes sure that the terminal text is back to its normal color.
    print STDERR color("reset");
}

exit(0);
# ====

__DATA__

sort.pl  [OPTIONS] [OPTIONS_TO_GNU_SORT]   FILENAME
  * Can handle compressed files (bz2 / gzip)
  * Can handle header lines.
  * [OPTIONS]: Options for this program (see below)
  * [OPTIONS_TO_GNU_SORT]: Options that are honored by the default unix "sort". "man sort" to see what they are.
  * FILENAME must come *last*.
  
  * Multiple filenames are not necessarily correctly supported.

by Alex Williams, 2009-2015

This is basically just a wrapper to GNU "sort" that adds the option of having a header line.

It does have one limitation: it can only sort ONE file while dealing with the header properly (whereas GNU sort cats all the command line arguments together).

See the examples below for more information.

CAVEATS:

* The very LAST argument must be the filename to sort.

* Only properly handles ONE file. Beware of passing in multiple files!

* The default delimiter is a TAB. (Unlike UNIX sort, which uses any whitespace.)

* sort -k 1,1 is DIFFERENT from sort -k 1   . (That is also true in UNIX sort.)

OPTIONS:

  --delim or -d = DELIMITER   (Default: tab)
     Sets the file delimiter to DELIMITER.

  Other options are just passed through to GNU sort.

EXAMPLES:

sort.pl somefile

cat somefile.gz | sort.pl -d 'Z' --rev

cat somefile.bz2 | tail -n 10 | sort.pl --rev

sort_uniq.pl -h 1 myfile.gz > out_with_header_still_at_top.txt

KNOWN BUGS:

  None known

TO DO:

  Nothing yet.


--------------

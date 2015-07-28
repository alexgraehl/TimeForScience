#!/usr/bin/perl

#@COMMENT@ flatten.pl 'flattens' a list of 'key' items in the first column and 'value' items in subsequent columns into a two-column edge list. Example: A B C D would turn into four lines: A B, A C, and A D. The opposite of flatten.pl is "expand.pl". Frequency-of-use rating: 8/10.
#@COMMENT@

use strict;

my @files;
my @cols     = ();
my $rev      = 0;
my $delim    = "\t";
my $outDelim = undef;
my $inc      = 1;
my $blanks   = 1;
my $all_cols = 0;
my $printOriginal = 0;

my $allowSingletons = 1; ## By default, we allow "singletons" with no tab-delimiter after them. If we put '-s' though, then anything without a delimiter afterwards gets kicked out

while(@ARGV)
{
    my $arg = shift @ARGV;

    if($arg eq '--help') {
	print STDOUT <DATA>;
	exit(0);
    }
    elsif($arg eq '-k') {
	$arg = shift @ARGV;

	if($arg eq 'all') {
	    $all_cols = 1;
	}
	else {
	    push(@cols, int(shift(@ARGV))-1);
	}
    }
    elsif($arg eq '-all') {
	$all_cols = 1;
    }
    elsif(($arg eq '-rev') or ($arg eq '--rev')) {
	$rev = 1;
    }
    elsif($arg eq '-d') {
	$delim = shift @ARGV;
    }
    elsif($arg eq '-do' or $arg eq '--do') {
	$outDelim = shift @ARGV;
    }
    elsif($arg eq '-i') {
	$inc = shift @ARGV;
    } elsif($arg eq '--ns') {
	$allowSingletons = 0;

    } elsif($arg eq '--orig') {
	$printOriginal = 1;
	
    } elsif($arg eq '-blanks' or $arg eq '-b') {
	$blanks = 0;
    }
    elsif((-f $arg) or (-l $arg) or ($arg eq '-')) {
	push(@files, $arg);
    }
    else {
	die("Bad argument '$arg' given.");
    }
}

if(scalar(@cols) == 0) {
    push(@cols, 0);
}

if (!defined($outDelim)) {
    $outDelim = $delim;
}

if(scalar(@files) == 0) { push(@files, '-'); }

foreach my $file (@files) {
    my $fin;

    open($fin, $file) or die("Could not open file '$file'");

    while(my $line = <$fin>) {
	next if ($line !~ /\S/); ## Always skip lines that are completely blank...
	my @tuple = split($delim, $line);
	chomp($tuple[$#tuple]);

	my $theOrigToPrint = undef;
	{
	    my $lineCopy = $line;
	    chomp($lineCopy);
	    $theOrigToPrint = ($printOriginal) ? (${lineCopy} . "\t") : '';
	}

	for(my $i = 0; $i < ($all_cols ? scalar(@tuple) : scalar(@cols)); $i++) {
	    my $col = ($all_cols) ? $i : $cols[$i];
	    my $key = splice(@tuple, $col, 1);

	    if (scalar(@tuple) == 0 && $allowSingletons) {
		## Ok, this line had only an item and then *no* tab after it. So there was only one element for this line!
		## i.e., the line was just "SOMETHING" and there was no delimiter afterward
		print $theOrigToPrint . $key . $outDelim . "\n"; ## we print this item, which is just a key but *no* delimiter, with a delimiter, even though it did not have one originally.
	    }
	    
	    for(my $i = 0; $i < scalar(@tuple); $i += $inc) {
		my $output = $key;
		for(my $j = 0; $j < $inc; $j++) {
		    my $data = $tuple[$i + $j];
		    if($blanks or $data =~ /\S/) {
			if(not($rev)) {
			    $output = $output . $outDelim . $data; ## REGULAR ORDER!!
			} else {
			    $output = $data . $outDelim . $output; ## REVERSED ORDER!
			}
		    }
		}
		print $theOrigToPrint . $output . "\n";
	    }
	    
	}
    }
    close($fin);
}

exit(0);

__DATA__
syntax: flatten.pl [OPTIONS] [FILE1 | < FILE1] [FILE2...]

Prints each tab-delimited data item in each row with its key seperately.
The first text field of data is assumed to be the key.

In other words, for each row like this:
"Alpha   Beta   Gamma  Delta  Beta"
It prints out:
   Alpha  Beta
   Alpha  Gamma
   Alpha  Delta
   Alpha  Beta

The "opposite" of flatten.pl is expand.pl .

OPTIONS are:

-k COL: Set the key column to COL (default is 1).  If COL equals 'all' then
        every column is treated as the key so that all pairwise combinations
        are printed.

-d DELIM: Set the input field delimiter to DELIM (default is <tab>).

-do OUTPUT_DELIM: When outputting, use this delimiter instead of the input one. Defaults to -d DELIM otherwise.

--ns: Disallow "singletons" (cases in which there is an item on a line, but no delimiter afterward)
    By default, we *do* print these items. But with --ns, we do not.

-i INC: increment through the members by INC (default is 1).

--rev: Print in reverse order.

--orig: Also print the original line (without ending whitespace/linebreaks), at the beginning of each line

-b: Do *not* print blank elements (default prints them). Totally blank lines are *always* skipped.

-all: Same as '-k all'



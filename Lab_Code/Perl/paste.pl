#!/usr/bin/env perl

use strict;

my $verbose   = 1;
my $delim     = "\t";
my $delim_in  = undef;
my $delim_out = undef;
my @is_file;
my @columns;
my @files;
my $num_files = 0;

# This is a bit of a mess but it's basically "getopt"
while (@ARGV) {
	my $arg = shift @ARGV;
	if ($arg eq '--help') {
		print STDOUT <DATA>;
		exit(0);
	} elsif ((-f $arg) or (-l $arg) or ($arg eq '-')) {
		$num_files++;
		push(@files, $arg);
		push(@is_file, $num_files);
		push(@columns, $arg);
	} elsif ($arg eq '-q') {
		$verbose = 0;
	} elsif ($arg eq '-d') {
		$delim = shift(@ARGV);
	} elsif ($arg eq '-di' or $arg eq '--di') {
		$delim_in = shift(@ARGV);
	} elsif ($arg eq '-do' or $arg eq '--do') {
		$delim_out = shift(@ARGV);
	} elsif ($arg eq '-nd') {
		# No delimiters
		$delim_in = '';
		$delim_out = '';
	} else {
		push(@is_file, 0);
		$arg =~ s/([^\\])\\t/$1\t/g;
		$arg =~ s/^\\t/\t/;
		$arg =~ s/([^\\])\\n/$1\n/g;
		$arg =~ s/^\\n/\n/;
		$arg =~ s/\\\\/\\/g;
		push(@columns, $arg);
	}
}

$delim_in  = defined($delim_in)  ? $delim_in  : $delim;
$delim_out = defined($delim_out) ? $delim_out : $delim;

if(scalar(@files) == 0) {
	push(@files,'-'); # Assume STDIN if there is no file
}

my @fins;
my @blanks;
for(my $f = 0; $f <= $#files; $f++) {
    open($fins[$f], $files[$f]) or die("Can't read file '$files[$f]'");
    $blanks[$f] = [];
}
my $done = 0;
while(not($done)) {
    my @tokens;
    my $num_file_tokens = 0;
    $done = 1;
    for(my $c = 0; $c < @columns; $c++)
    {
        my $f = $is_file[$c];
        if($f)
        {
           $f--;
           my $fin = $fins[$f];
	   # Alex: there was a bug here: if you have the line below
	   # check for "not(eof($fin))", then it actually CLIPS OFF
	   # the last line when pasting.
           if((my $line = <$fin>)) { # and not(eof($fin)))
               my @tuple = split($delim_in, $line);
               chomp($tuple[$#tuple]);
               if(scalar(@tuple) > scalar(@{$blanks[$f]})) {
		       my @blankList = '' x scalar(@tuple); # make blank list of size '@tuple'
		       $blanks[$f] = \@blankList; # pass in a reference
	       }
               push(@tokens, \@tuple);
               $done = 0; # keep going!
               $num_file_tokens++;
           } else {
		   push(@tokens, $blanks[$f]);
           }
        }
        else {
           push(@tokens, [$columns[$c]]);
        }
    }
    if (scalar(@files) == 0 or $num_file_tokens > 0) {
	    printMultiLists($delim_in, $delim_out, @tokens); # this function is defined right below here
    }
}
for(my $f = 0; $f < scalar(@files); $f++) {
    close($fins[$f]);
}

exit(0);

sub printMultiLists {
   my ($delim_within, $delim_between, @lists) = @_;
   for(my $i = 0; $i < @lists; $i++) {
	   print STDOUT ($i > 0 ? $delim_between : "") . join($delim_within, @{$lists[$i]});
   }
   print STDOUT "\n";
}

__DATA__
syntax: paste.pl [OPTIONS]  input1 [input2] [input3]...

'Pastes' together files (like UNIX paste) or literal text. Prints to standard out.

EXAMPLES:

paste.pl -do ''  '>>>>' INFILE  > OUTFILE
  Adds >>> to the BEGINNING of every line in INFILE, printing this to OUTFILE.
  Note the part with -do '', which suppresss the normal (tab) delimiting between inputs.

paste.pl INFILE 'UNKNOWN_SOURCE' > OUTFILE
  Puts a tab (the default delimiter) and then "UNKNOWN_SOURCE" at the END of every line.

paste.pl 'LEFT' 'MORE' FILE1 FILE2 'RIGHT'
  Pastes a bunch of files and text together, printing to the console.

WARNING: You cannot paste text that is ALSO a filename!! (Files override literal text.)
 Therefore:  paste.pl FILE "FILE"
 will paste TWO copies of the FILE, not the literal text "FILE", assuming "FILE" is an actual filename.
 Workaround:
 Instead of:   paste.pl  FILE  "FILE"
 you can do:   paste.pl  FILE  "__FILE__" (assuming "__FILE__" is not also a valid filename)

OPTIONS:

-q: Quiet mode (default is verbose)

-di: Set the input delimiter. I am not sure why there even IS
  a "delim_in" since it does not seem like it should matter.
  Paste.pl normally should not care about input delimters... I think.

-do: Set the output delimiter. Default is tab.
     Set this to '' in order to output nothing at
     all between the pasted value and the next item.

-d: Set the delimiter for both input AND output simultaneously.
    Same as:   -di SOMETHING  -do SOMETHING

-nd: NO delimiter: input & output delims are both set to ''
     so as to have no tabs separating items. Same as -d ''.



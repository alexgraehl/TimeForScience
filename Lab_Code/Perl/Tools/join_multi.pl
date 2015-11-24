

#@COMMENT@ join_multi.pl allows you to merge multiple files in different orders into one big matrix file. For example, if you have per-gene counts for 10 experiments, but sometimes an experiment has a different set of genes, you would run join_multi.pl EXP1.txt EXP2.txt ... EXPN.txt > out.matrix.txt

use lib "$ENV{MYPERLDIR}/lib"; use lib "$ENV{TIME_FOR_SCIENCE_DIR}/Lab_Code/Perl/LabLibraries"; require "libfile.pl";
use lib "$ENV{MYPERLDIR}/lib"; use lib "$ENV{TIME_FOR_SCIENCE_DIR}/Lab_Code/Perl/LabLibraries"; require "liblist.pl";
use lib "$ENV{MYPERLDIR}/lib"; use lib "$ENV{TIME_FOR_SCIENCE_DIR}/Lab_Code/Perl/LabLibraries"; require "libstd.pl";
use lib "$ENV{MYPERLDIR}/lib"; use lib "$ENV{TIME_FOR_SCIENCE_DIR}/Lab_Code/Perl/LabLibraries"; require "libstring.pl";

use strict;
use warnings;
use strict; use warnings; use Getopt::Long;

use List::MoreUtils;
  
sub quitWithUsageError($) { print($_[0] . "\n"); printUsageAndQuit(); print($_[0] . "\n"); }
sub printUsageAndQuit() { printUsageAndContinue(); exit(1); }
sub printUsageAndContinue() {    print STDOUT <DATA>; }
sub debugPrint($) {   ($isDebugging) and print STDERR $_[0]; }

sub closeSmartFilehandle($) { my($handle)=@_; if ($handle ne *STDIN) { close $handle; } }# Don't close STDIN, but close anything else!
sub openSmartAndGetFilehandle($) {
    # returns a FILEHANDLE. Can be standard in, if the 'filename' was specified as '-'
    # Transparently handles ".gz" and ".bz2" files.
    # This is the MARCH 6, 2013 version of this function.
    my ($filename) = @_;
    if ($filename eq '-') { return(*STDIN); } # <-- RETURN!!!
    my $reader;
    if    ($filename =~ /[.]gz$/)  { $reader = "gzip --stdout $filename |"; }     # Un-gzip a file and send it to STDOUT.
    elsif ($filename =~ /[.]bz2$/) { $reader = "bzcat $filename |"; }    # Un-bz2 a file and send it to STDOUT
    elsif ($filename =~ /[.]zip$/) { $reader = "unzip -p $filename |"; } # Un-regular-zip a file and send it to STDOUT with "-p": which is DIFFERENT from -c (-c is NOT what you want here). See 'man unzip'
    else                           { $reader = "$filename"; }  # Default: just read a file normally
    my $fh;
    open($fh, "$reader") or die("Couldn't read from <$filename>: $!");
    return $fh;
}

# Flush output to STDOUT immediately.
$| = 1;

my $addFilenamesToHeader = 0;
my $headers = 1; # default: 1 header line

$Getopt::Long::passthrough = 1; # ignore arguments we don't recognize in GetOptions, and put them in @ARGV
GetOptions("help|?|man" => sub { printUsageAndQuit(); }
	   , "q" => sub { $verbose = 0; }
#	   , "k=i" => \$key_col
	   , "d=s" => \$delim
	   , "o=s"  => \$blank
	   , "m=i"  => \$min
	   , "h=i"  => \$headers
	   , "fnames!"  => \$addFilenamesToHeader
	   #, "i|ignore-case!" => \$shouldIgnoreCase
#	   , "debug!" => \$isDebugging
    ) or printUsageAndQuit();

my $key_col        = 1 #$args{'-f'};
my $blank_placeholder = '___@@@_BLANK_VALUE_@@@___';

my $numUnprocessedArgs = scalar(@ARGV);
($numUnprocessedArgs == 2) or quitWithUsageError("[Error] in arguments! You must send at least TWO filenames to join_multi.pl , or else there is nothing to join in the first place.");

my @files = @ARGV; # ok, these are files
my $num_files = scalar(@files);

my @all_headers = ();

# The number of files the key appears in:
my %count;
for (my $j = 0; $j < scalar(@files); $j++) {
	$verbose and print STDERR "(", $j+1, "/" . scalar(@files) . "). ", "Collecting every key from file '$files[$j]'...";
	my $file = openSmartAndGetFilehandle($files[$j]);
	my $lineNum = 0;
	while (my $line = <$file>) {
		$lineNum++;
		my @tuple = split($delim, $line);
		chomp($tuple[$#tuple]);
		my @key_cols = (0);
		my @value_cols;
		for (my $i = 0; i < scalar(@tuple); i++) {
			my $is_a_key_column = List::MoreUtils(any($_ == $i), @key_cols);
			if (not $is_a_key_column) {
				push @value_cols, $i;
				push @all_headers, 
			}
		}
	}
	closeSmartFilehandle($file);
	foreach my $key (@file_keys) {
		if (!exists($count{$key})) {
			$count{$key} = 1;
			push(@keys_in_order, $key);
		} else {
			$count{$key}++;
		}
	}
	$verbose and print STDERR " done.\n";
}

my $num_keys_total = scalar(keys(%count));
my %row;
my @data;
my $num_keys_kept = 0;

foreach my $key (@keys_in_order) {
   if($count{$key} >= $min) {
      $data[$num_keys_kept][0] = $key;
      $row{$key} = $num_keys_kept;
      $num_keys_kept++;
   }
}

$verbose and print STDERR "$num_keys_kept keys present in $min or more files kept (out of $num_keys_total total found).\n";

my @blanks;

#my @headers = ();

for(my $j = 0; $j < scalar(@files); $j++) {
   $verbose and print STDERR "(", $j+1, "/" . scalar(@files) . "). ",
                             "Collecting data from file '$files[$j]'...";
   my $file = &openFile($files[$j]); # OPEN THE FILE HERE
   my %keys_not_seen = %count;
   $blanks[$j] = &duplicate($max_tuple[$j], $blank_placeholder);
   my $key_cols  = defined($cols[$j]) ? $cols[$j] : \@global_cols;
   my @sorted_key_cols = sort { $a <=> $b; } @{$key_cols};

   my $line_no = 0;
   while(my $line = <$file>) {
      $line_no++;
      # my @tuple = split($delim, $line);
      my $tuple = &mySplit($delim, $line);
      my $last = scalar(@{$tuple}) - 1;
      chomp($$tuple[$last]);
      my $key = &extractKey($tuple, $key_cols, \@sorted_key_cols);
      if(exists($row{$key})) {
         if(exists($keys_not_seen{$key})) {
            my $i = $row{$key};
            push(@{$data[$i]}, &pad($tuple, $blanks[$j]));
            delete($keys_not_seen{$key});
         }
      }
   }

   foreach my $key (keys(%keys_not_seen)) {
      if(exists($row{$key})) {
         my $i = $row{$key};
         push(@{$data[$i]}, $blanks[$j]);
      }
   }
   close($file); # close the file here!
   $verbose and print STDERR " done.\n";
}

$verbose and print STDERR "Printing out the combined data...\n";

for(my $i = 0; $i < $num_keys_kept; $i++) {
   my $key = defined($data[$i][0]) ? $data[$i][0] : $blank_placeholder;
   my $printable = $key;
   for(my $j = 1; $j <= scalar(@files); $j++) {
      my $val_list = defined($data[$i][$j]) ? $data[$i][$j] : $blanks[$j-1];
      my $values   = join($delim, @{$val_list});
      $printable .= $delim . $values;
   }
   $printable =~ s/$blank_placeholder/$blank/g; # replace the placeholder with the blank value again
   print STDOUT $printable, "\n";
}
$verbose and print STDERR " done.\n";

exit(0);

sub parseColsFromArgs {
   my ($args) = @_;
   my @cols;
   my $fileno = undef;
   while(@{$args}) {
       my $arg = shift @{$args};
       if ($arg =~ /^-(\d+)/) {
	   $fileno = int($1) - 1;
	   $arg = shift @{$args};
	   my @col_list = split(',', $arg);
	   foreach my $col (@col_list) {    push(@{$cols[$fileno]}, $col - 1); }
       }
   }
   return \@cols;
}

sub extractKey {
   my ($tuple, $cols, $ordered_cols) = @_;
   my $key = '';
   for(my $i = 0; $i < scalar(@{$cols}); $i++) {
      $key .= (($i > 0) ? "\t" : "") . $$tuple[$$cols[$i]];
   }
   for(my $i = scalar(@{$cols}) - 1; $i >= 0; $i--) {
      splice(@{$tuple}, $$ordered_cols[$i], 1);
   }
   return $key;
}

__DATA__
syntax: join_multi.pl [OPTIONS] FILE1 FILE2 [FILE3 FILE4 ...]

Joins multiple files together, saving you from having to call
join.pl over and over. Note that this is a specific kind of join,
and the results from join_multi.pl will most likely be *DIFFERENT*
from the results of regular join.pl .

Unclear how it handles headers. Seems to work, though.

If you want to merge files together, for example, to collect all
the runs of a various experiments into one file, you can use
join_multi.pl like this:

Example:

You have a file:

Experiment1  0.11
Experiment2  2.22

And a second file:

Experiment1  1.111
Experiment5  5.55
Experiment6  6.6

If you run join_multi.pl FILE1 FILE2 on them, you will get:

Experiment1  0.11   1.111
Experiment2  2.22
Experiment5         5.55
Experiment6         6.6

(Note that join_multi.pl may also insert an extra tab between
 the key and the first value, so you may need to check for this
 and use "cut" appropriately.)

CAVEATS:

BUG REPORT: IT DOES NOT HANDLE HEADERS PROPERLY!!!!!!!! The second file has a header that just
  goes down somewhere randomly in the file, instead of at the top.

join_multi.pl has weird behavior when you use a totally empty file.
It will usually add two tabs between fields instead of one. To fix this,
echo '' >> PREVIOUSLY_EMPTY_FILE before you join_multi.pl it.

OPTIONS are:

-q: Quiet mode (default is verbose)

-k COL: The key column MUST be column 1 right now, so this does not work.
        Otherwise, it would set the key column to COL (default is 1).

-d DELIM: Set the field delimiter to DELIM (default is tab).

-m MIN: Set the minimum number of files that an entry has to 
        exist in to MIN.  The default is 1 which means the entry
        has to appear in at least one file.

-o BLANK_OUTPUT: Set the blank character to BLANK (default is empty).
                 This is whatever gets printed when there is NO MATCH in a file.
                 You could set it to (for example) -o "NA".

--addfilenames or -a: Add filenames to the header line.



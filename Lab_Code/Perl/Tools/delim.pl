#!/usr/bin/perl

##############################################################################
##############################################################################
##
## del.pl
##
##############################################################################
##############################################################################
##
## Written by Josh Stuart.
##
##  Email address: jstuart@ucsc.edu
##          Phone: (650) 725-7612
##
## Postal address: Department of Bioengineering, UCSC
##                 1156 High Street
##                 Santa Cruz, CA 95060
##
##       Web site: http://www.soe.ucsc.edu/~jstuart/
##
##############################################################################
##############################################################################

require "$ENV{HOME}/develop/perl/lib/libfile.pl";

use strict;
use warnings;

# Flush output to STDOUT immediately.
$| = 1;


my @flags   = (
                  [    '-q', 'scalar',     0,     1]
                , [    '-k', 'scalar',     1, undef]
                , [    '-d', 'scalar',  "\t", undef]
                , ['--file', 'scalar',   '-', undef]
              );

my %args = %{&parseArgs(\@ARGV, \@flags)};

if(exists($args{'--help'}))
{
   print STDOUT <DATA>;
   exit(0);
}

my $verbose = not($args{'-q'});
my $col     = int($args{'-k'}) - 1;
my $delim   = $args{'-d'};
my $file    = $args{'--file'};

my $filep;
open($filep, $file) or die("Could not open file '$file' for reading");
while(<$filep>)
{
   my @x = split($delim, $_);
   chomp($x[$#x]);
   my $item  = $x[$col];
}
close($filep);

exit(0);

__DATA__
syntax: del.pl [OPTIONS]

OPTIONS are:

-q: Quiet mode (default is verbose)

-e: Empty the trash

-t TMPDIR: Set the temporary directory to TMPDIR (default is /tmp).


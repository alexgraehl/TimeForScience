#!/usr/bin/perl

use strict;

use lib "$ENV{MYPERLDIR}/lib"; use lib "$ENV{TIME_FOR_SCIENCE_DIR}/Lab_Code/Perl/LabLibraries"; require "libmap.pl";

my $runs_mak = &getMapDir('Templates') . '/Make/runs.mak';
my $makefile = './Makefile';

my $cmd = "ln -s $runs_mak $makefile";

print STDERR "$cmd";
system("$cmd");
print STDERR "\n";


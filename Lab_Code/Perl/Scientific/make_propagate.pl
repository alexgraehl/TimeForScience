#!/usr/bin/perl

use strict;

use lib "$ENV{MYPERLDIR}/lib"; use lib "$ENV{TIME_FOR_SCIENCE_DIR}/Lab_Code/Perl/LabLibraries"; require "libmap.pl";

my $propagate_mak = &getMapDir('Templates') . '/Make/propagate.mak';
my $makefile = './Makefile';

my $cmd = "ln -s $propagate_mak $makefile";
print STDERR "$cmd";
system("$cmd");
print STDERR "\n";


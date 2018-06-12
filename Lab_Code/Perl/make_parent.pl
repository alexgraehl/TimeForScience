#!/usr/bin/env perl

use lib "$ENV{MYPERLDIR}/lib"; use lib "$ENV{TIME_FOR_SCIENCE_DIR}/Lab_Code/Perl/LabLibraries"; require "libmap.pl";

my $makefile_template = &getMapDir('Templates') . '/Make/parent.mak';
my $pwd = "$ENV{PWD}";

system("ln -s $makefile_template $pwd/Makefile");


#!/usr/bin/perl -w

# This program removes ONE-LINE comments. It won't remove comments that are at the end of a line.
# It only removes lines that are ENTIRELY comments and nothing else.

# It also doesn't understand multi-line comments.


use strict;

use File::Basename;
use Getopt::Long;

my $comment_symbol = '#'; # default comment symbol is the hash mark
my $strip_blank_lines = 0;

GetOptions("help|?|man" => sub { print STDOUT <DATA>; exit(0); }
	   , "s|symbol=s" => \$comment_symbol
	   , "semicolon|lisp" => sub { $comment_symbol = ';'; }
	   , "hash|perl" => sub { $comment_symbol = '#'; }
	   , "c-style" => sub { $comment_symbol = '//'; }
	   , "b|strip_blank_lines" => sub { $strip_blank_lines = 1; }
	   );

#$symbol = s/\//\\\//g; # replace '/' with '\/' so it works with egrep (actually not necessary, it turns out)

my $cmd = '';

if ($strip_blank_lines) {
    $cmd .= "egrep -v '^([ |\t]*)\$' | ";
}

$cmd .= "egrep -v '^([ |\t]*)${comment_symbol}.*'";

system($cmd);
print $cmd;

#while (<>) { # read from STDIN
#    if ($_ =~ "^\s*${comment_symbol}" || ($strip_blank_lines && ($_ =~ "^\s*\$")) ) {
#	# don't print comment lines... skip to the next line
#	# also, if the user has set the program to strip blank lines, we don't print those either
#    } else {
#	print STDOUT $_;
#    }
#}

exit(0);

__DATA__

syntax: strip_comments.pl [OPTIONS] < STDIN

description: This program strips away one-line comment lines.
    It reads from STD, and outputs to STDOUT.
    The default comment character is the hash ('#').
    The main part is not perl--it uses egrep, so it should be fast.

-s=SYMBOL  (or --symbol=SYMBOL)
    Set the comment symbol to SYMBOL. Default is '#'.
    Can be multiple characters (like // in C).
       Examples:  -s=#  -s=;  -s=//
    
    Remember that this program will NOT handle multi-line comments (like /* */).

--perl
--lisp
--c-style
    These set the comment character to #, ;, or //, respectively.

-b: also strip out blank lines


--help: prints this help screen

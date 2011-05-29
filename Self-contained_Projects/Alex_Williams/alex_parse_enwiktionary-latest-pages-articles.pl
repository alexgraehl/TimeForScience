#!/usr/bin/perl


## Expecting to read from the file named:

## enwiktionary-latest-pages-articles.xml

## usage: ./programname.pl enwiktionary-latest-pages-articles.xml 

## expected tab-delimited output with three columns:
# line number //// type_of_thing (word or etymology) //// content
# 2376    WORD_11 portmanteau
# 2401    ETYMOLOGY1      From {{etyl|fr}} {{term|portemanteau|lang=fr}}, literally {{term|porte||carry}} + {{term|manteau|lang=fr||coat}}
# 2437    ETYMOLOGY2      Coined by [[w:Lewis Carrol|Lewis Carrol]] in [[s:Through The Looking Glass (And What Alice Found There)|Through The Looking Glass]] to describe the words he coined in [[w:Jabberwocky|Jabberwocky]].


use strict; use warnings;

my $next = 0;

my $nextLineIsAnEtymology = 0;

my $nextLineIsAPageTitleWord = 0;

my $lastWord = '';

my $lineNum = 0;
my $numEtymologiesFoundSinceLastWord = 0;
my $numWords = 0;
while (my $line = <>) {
    chomp($line);
    $lineNum++;
    if ($line =~ m{<title>(.*)</title>}) {
	$numWords++;
	$lastWord = $1;
	print "\n";
	print $lineNum . "\tWORD_$numWords\t" . $lastWord . "\n";
	$numEtymologiesFoundSinceLastWord = 0;
    }

    if ($nextLineIsAnEtymology) { # && ($numEtymologiesFoundSinceLastWord == 1)) {
	print $lineNum . "\tETYMOLOGY${numEtymologiesFoundSinceLastWord}\t";
	print $line . "\n";
	$nextLineIsAnEtymology = 0;
    }


    if ($line =~ /^===Etymology/) {
	$nextLineIsAnEtymology=1;
	$numEtymologiesFoundSinceLastWord++;
	#print "$lineNum: Expecting an etymology number $numEtymologiesFoundSinceLastWord on the next line for \"$lastWord\"...\n";

    }

}

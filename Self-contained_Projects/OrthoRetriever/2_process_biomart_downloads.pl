#!/usr/bin/perl -w

use strict;
use warnings;

use Term::ANSIColor;

sub fileHasContents($) {
    my ($filename) = @_;
    return (-f $filename and ((-s $filename)>0));
}

sub agwSystem($) {
    print colored("Executing this command: $_[0]\n", "black on_green");
    my $result = system($_[0]);
    print colored("Result code was: $result\n", "green");
}

my @mapFiles = split('\n', `ls -1 DL_BIOMART*map*`);

print "MAP FILES:\n";
print (join(", ", @mapFiles) . "\n");

my @masterFiles = split('\n', `ls -1 DL_BIOMART_PRIMARY*`);

print "PRIMARY / MASTER LIST FILES:\n";
print (join(", ", @masterFiles) . "\n");


#my $MAP_FROM_COL = 1; # column with the Ensembl name of the FROM species (indexed from ONE-- 1 is the first column, NOT zero!)
#my $MAP_TO_COL   = 4; # column with the Ensembl name of the TO species (indexed from ONE-- 1 is the first column, NOT zero!)
#my $MAP_TYPE_COL = 5; # "many-to-one" or "one-to-one" or whatever

foreach my $ff (@mapFiles) {
    my $annotFinal = "ANNOT_FINAL_${ff}";
    my $ffWithoutSpeciesName = "${ff}.fixed.header.tmp";
    my $AAA = "aaa.tmp";
    my $BBB = "bbb.tmp";
    my $CCC = "ccc.tmp";

    $ff =~ /DL_BIOMART_([^_]+)_map_to_([^_]+)/i;
    my $from = $1;
    my $to   = $2;

    print "The file <$ff> maps orthologs $from --> $to\n"; 
    
    my $fromPrimary = "DL_BIOMART_PRIMARY_${from}";
    my $toPrimary = "DL_BIOMART_PRIMARY_${to}";
    (-f $fromPrimary and -f $toPrimary) or die "From primary and to-primary files did not exist! Make sure <$fromPrimary> and <$toPrimary> are both valid files.";
    
    agwSystem(qq{cat ${ff} | perl -p -e 's/[^\t]+[ ](Ensembl Gene ID)/\$1/g' > ${ffWithoutSpeciesName} }); # remove the species-specific name!
    agwSystem(qq{cat $ffWithoutSpeciesName | select.pl -k 4 -lne '' > $AAA}); # only lines with NON-BLANK map-to-the-other-species entries!
    agwSystem(qq{join.pl -1 4 -2 1 -o "NO_ANNOTATION--THIS_SHOULD_NOT_HAPPEN" $AAA $toPrimary > $BBB }); # replace the species-specific Ensembl Gene ID field with the literal text "Ensembl Gene ID"
    
    my $uniqLineLenCount = `count_items_per_line.pl $BBB | uniq | wc -l`;    chomp($uniqLineLenCount);
    print STDERR "<$uniqLineLenCount>\n";
    ($uniqLineLenCount eq '1') or die "Problem: some of the output lines in <$BBB> had different lengths (number of tab-delimited elements) after joining them!";
    
    agwSystem(qq{join.pl -1 2 -2 1 -o "NO_ANNOTATION--THIS_SHOULD_NOT_HAPPEN" $BBB $fromPrimary > $CCC});

    agwSystem(qq{cat $CCC | perl -p -e 's/ [(]\\d+ of \\d+[)]//g' > $annotFinal}); # delete the " (2 of 3)" kind of thing in some gene names. This ONLY appears in the common names for the homology, not in the "real" data!

    my $annotUniqLineLenCount = `count_items_per_line.pl $annotFinal | uniq | wc -l`;    chomp($annotUniqLineLenCount);
    print STDERR "<$uniqLineLenCount>\n";
    ($annotUniqLineLenCount eq '1') or die "Problem: some of the output lines in <$annotFinal> had different lengths (number of tab-delimited elements) after joining them!";
    #die "ok done testing";
}

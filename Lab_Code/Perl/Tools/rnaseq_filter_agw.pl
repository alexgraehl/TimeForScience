#!/usr/bin/perl

## This script handles a few basic processing steps for an aligned SAM or BAM file.
## by Alex Williams, May 2012




use strict;  use warnings;

use Getopt::Long;
use File::Basename;

my $assumeFileIsAlreadySorted = 0; ## By default, assume the input SAM/BAM file will still require sorting.
my $shouldSort = 1;
my $shouldFilterMappedOnly = 1;
my $isDryRun = 1;


sub printUsageAndQuit() { ## Function for printing the "you used the wrong arguments to this program" error message
    print STDOUT <DATA>;
    exit(0);
}

sub datePrint($) { my $d = `date`; chomp($d); print $d . ":\t" . $_[0]; } ## Timestamping function that is like "print" except with a timestamp at the beginning of the line. Use it instead of "print."

sub alexSystemCall($) {
    my ($cmd) = @_;
    if ($isDryRun) {
	datePrint("Dry run of command: $cmd\n");
    } else {
	datePrint("Executing command: $cmd\n");
    }
    if (!$isDryRun) { system($cmd); }
}



my $SAMTOOLS_PATH = `which samtools`;  chomp($SAMTOOLS_PATH);
my $PICARD_PATH   = `which SortSam.jar`;  chomp($PICARD_PATH);
my $GIGABYTES_FOR_PICARD = 2;

if (length($SAMTOOLS_PATH) <= 1) {
    die "Could not find <samtools> in the \$PATH. Make sure the <samtools> executable is somewhere in your \$PATH. You should be able to type \"samtools\" and run the command; if that is not the case, then this program won't run either. You can install samtools using apt-get.\n";
}

if (length($PICARD_PATH) <= 1) {
    die "Could not find <SortSam.jar> (part of the Picard suite) in the \$PATH. Make sure the <SortSam.jar> java applet is somewhere in your \$PATH. You should be able to type \"SortSam.jar\" and get a weird java error message; if that is not the case, then this program won't run either. For installation instructions, check http://picard.sourceforge.net/.\n";
}



GetOptions("help|?|man"        => sub { printUsageAndQuit(); }
	   , "sort!" => \$shouldSort ## specify "--nosort" to avoid sorting
	   , "mapfilter!" => \$shouldFilterMappedOnly ## specify "--nomapfilter" to avoid filtering out the unmappable reads
	   , "dry|dryrun|dryRun|dry_run|dry-run" => sub { $isDryRun = 1; }
    ) or printUsageAndQuit();

if (scalar(@ARGV) < 1) { die "ARGUMENT ERROR: This script requires at least one filename---you need to give it a BAM or SAM file to operate on.\n"; }

print STDOUT "The following files will be processed by rnaseq_filter_agw:\n";
foreach my $tmp (@ARGV) {
    print STDOUT "  - ${tmp}\n";
}

my $numFilesSuccessfullyProcessed = 0;
my $numFilesNotOK = 0;

#use diagnostics;

foreach my $file (@ARGV) {
    my $fileOK = 1; # Assume the file is something we can read, unless we hear otherwise...
    my $isSam = ($file =~ m/[.]sam$/i);
    my $isBam = ($file =~ m/[.]bam$/i);
    if (not -r $file) { print "File <$file> is not readable by this user.\n"; $fileOK = 0; }
    if (-z $file) { print "File <$file> is of zero length.\n"; $fileOK = 0; }
    if (not -e $file) { print "File <$file> does not exist.\n"; $fileOK = 0; }
    if (!$isSam and !$isBam) { print "File <$file> doesn't end with .bam or .sam. Input filenames must be <something.sam> or <something.bam>."; $fileOK = 0; }
    if (!$fileOK) {
	datePrint("We could not read valid SAM/BAM data from <$file>  (Please double-check this file path.)\n");
	$numFilesNotOK++;
	next; # Skip to next iteration of the loop
    }

    datePrint("Now processing the file <$file>...\n");
    my $filteredFile = $file . ".filtered.tmp";
    my $sortedFile = $file . ".sorted.tmp";
    
    my $fileTypeOption = ($isSam ? " -S " : " "); ## SAM files require the "-S" option in order to be viewed in samtools. By default, a BAM file (no option) is expected.
    my $countAllReadsCmd = qq{${SAMTOOLS_PATH} view $fileTypeOption $file | wc -l };

    my $sortCmd = (qq{java -Xmx${GIGABYTES_FOR_PICARD}g -jar ${PICARD_PATH} }
		   . qq{ INPUT=${file} } ## Picard SortSam.jar accepts both SAM and BAM files as input!
		   . qq{ SORT_ORDER=coordinate }
		   . qq{ OUTPUT=${sortedFile} });

    my $filterMappedOnlyCmd = (qq{samtools view ${fileTypeOption} -h -F 4 $file }
			       . qq{ > $filteredFile });
    

    alexSystemCall($countAllReadsCmd);
    alexSystemCall($sortCmd);
    alexSystemCall($filterMappedOnlyCmd);
    
    #system(qq{samtools faidx $genomeFastaFile}); ## generate an index file
    #print STDOUT "Done. Wrote <$faiFile> to the filesystem.";    

    datePrint("[Done] with <$file>.\n\n");
    $numFilesSuccessfullyProcessed++;
}

# open FILE, ">>", $browserTrackDescriptionFile or die $!; ## APPEND TO THE FILE!!!
# print FILE "\n";
# print FILE browserTrackString("bam", ${sortBamFullFilename});
# print FILE browserTrackString("bigWig", ${bigWigOutFile});
# print FILE "\n";
# close(FILE);

datePrint("[DONE]\n\n");


__DATA__

convert_SAM_or_BAM_for_Genome_Browser.pl --fasta=<genome_fasta_file> <INPUT BAM / SAM FILE>

Processes a SAM or BAM alignment file to generate files for the UCSC Genome Browser.

OPTIONS:
   --nosort : Do not sort the input file. Not usually a useful option.
   --nomapfilter: Do not remove the unmapped reads.
            (By default, unmapped reads are removed.)

   --dryrun : Print what *would* be done, but do not actually perform any filtering.

This is a script that will ****************
************************
************************
************************

Requres the following additional software:
    * <samtools> must be installed (to filter the bam files)
        - You can install it with `apt-get install samtools`
    * <Picard> must be installed in order to properly sort the sam/bam files. samtools is not ideal because it does not set the sort flag at the top of the file.
        - Check at http://picard.sourceforge.net/

If you want a wiggle track, which you most likely do (if you do not, you can specify `--nowig`):
    * <genomeCoverageBed> must be installed
    * <wigToBigWig> must be installed

Inputs:
   *************

Output: 
 **********

These files can be viewed in the genome browser.


EXAMPLES:

 rnaseq_filter_agw *******************
 ??

Files that are generated from the input YOURFILE.sam:
 ** ?

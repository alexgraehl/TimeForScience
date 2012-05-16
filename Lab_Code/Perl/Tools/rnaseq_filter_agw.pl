#!/usr/bin/perl

## This script handles a few basic processing steps for an aligned SAM or BAM file.
## by Alex Williams, May 2012




use strict;  use warnings;

use Getopt::Long;
use File::Basename;

my $isDryRun = 0;


my $shouldSort = 1;
my $shouldFilterMappedOnly = 1;
my $shouldRemoveDupes = 1;
my $shouldCalculateSummary = 1;
my $shouldGenerateIndex = 1;
#my $genomeFastaFile = undef;

my $keepAllReads = 0;

my $SAMTOOLS_PATH = `which samtools`;  chomp($SAMTOOLS_PATH);
my $SORTSAM_PATH   = `which SortSam.jar`;  chomp($SORTSAM_PATH);
my $GIGABYTES_FOR_PICARD = 2;
my $MARKDUPLICATES_PATH = `which MarkDuplicates.jar`; chomp($MARKDUPLICATES_PATH);
my $ULIMIT_RESULT = 1024; ## result of running the shell command ulimit -n. Since this is a shell built-in, it can, for some reason, not be run like a real command, so backticks don't work. `ulimit -n`; chomp($ULIMIT_RESULT);
my $latest           = "latest_file_to_operate_on.tmp.bam"; ## This is a temporary symlink that gets updated

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

# sub verifyFastaIndexFileOrDie() {
#     if ($shouldGenerateIndex) {
# 	if (!defined($genomeFastaFile)) {
# 	    die "ARGUMENT ERROR: If you want to generate an index (.bai) file, you need to specify a fasta file on the command line using the --fasta=FILENAME option.\n"
# 		. "Or you can OMIT generating a fasta file by specifiying the option --noindex .\n"
# 		. "You did not specify a fasta file.\n";
# 	}
# 	if ((not -e $genomeFastaFile) or (not -r $genomeFastaFile)) {
# 	    die "ARGUMENT / FILE READING ERROR: We need a valid genome fasta file in order to generate an index index (.bai) file.\n"
# 		. "The file you specified cannot be read or does not exist. Here is the file you specified as the fasta file:\n"
# 		. "   >  $genomeFastaFile \n"
# 		. "Please double-check to make sure that file actually exists and is readable!\n"
# 		. "Alternatively, you can decide not to make an index at all by specifiying --noindex.\n";
# 	}
#     }
# }

sub verifyToolPathsOrDie() {
    # Below: Check for the location of some required binaries.
    # We depend on samtools & the Picard suite being installed already.
    if (length($SAMTOOLS_PATH) <= 1) { die "Could not find <samtools> in the \$PATH. Make sure the <samtools> executable is somewhere in your \$PATH. You should be able to type \"samtools\" and run the command; if that is not the case, then this program won't run either. You can install samtools using apt-get.\n"; }
    if (length($SORTSAM_PATH) <= 1) {  die "Could not find <SortSam.jar> (part of the Picard suite) in the \$PATH. Make sure the <SortSam.jar> java applet is somewhere in your \$PATH. You should be able to type \"SortSam.jar\" and get a weird java error message; if that is not the case, then this program won't run either. For installation instructions, check http://picard.sourceforge.net/.\n"; }
    if (length($MARKDUPLICATES_PATH) <= 1) {  die "Could not find <Markduplicates.jar> (part of the Picard suite) in the \$PATH. Make sure the <MarkDuplicates.jar> java applet is somewhere in your \$PATH. You should be able to type \"MarkDuplicates.jar\" and get a weird java error message; if that is not the case, then this program won't run either. For installation instructions, check http://picard.sourceforge.net/.\n"; }
}



GetOptions("help|?|man"        => sub { printUsageAndQuit(); }
	   , "sort!" => \$shouldSort ## specify "--nosort" to avoid sorting
	   , "mapfilter!" => \$shouldFilterMappedOnly ## specify "--nomapfilter" to avoid filtering out the unmappable reads
	   , "keepdupes" => sub { $shouldRemoveDupes = 0; }
#	   , "f|fasta" => \$genomeFastaFile ## used for generating the index
	   , "noindex" => sub { $shouldGenerateIndex = 0; }
	   , "keep|keepall!" => \$keepAllReads
	   , "dry|dryrun|dryRun|dry_run|dry-run" => \$isDryRun
    ) or printUsageAndQuit();

if (scalar(@ARGV) < 1) { die "ARGUMENT ERROR: This script requires at least one filename---you need to give it a BAM or SAM file to operate on.\n"; }

print STDOUT "The following files will be processed by rnaseq_filter_agw:\n";
foreach (@ARGV) {
    print STDOUT "  - $_\n";
}

if ($keepAllReads) {
    print STDOUT "We are running in the MOST CONSERVATIVE mode, where we DO NOT REMOVE ANY READS.\n";
    print STDOUT "   * Unmapped reads will be RETAINED.\n";
    print STDOUT "   * Duplicate-location reads will be RETAINED.\n";
    $shouldRemoveDupes = 0;
    $shouldFilterMappedOnly = 0;
}

verifyToolPathsOrDie();

my $numFilesSuccessfullyProcessed = 0;
my $numFilesNotOK = 0;

#use diagnostics;

my @outFiles = (); # No files generated yet.

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
    my $filteredFile     = "$file.filtered_mapped_only.tmp.bam";
    my $sortedFile       = "$file.sorted_by_coord.tmp.bam";
    my $dedupFile        = "$file.no_duplicates.tmp.bam";
    my $dedupExtraMetricsFile = "$file.no_duplicates.extra.txt";
    my $summaryStatsFile = "$file.summary.stats.txt";
    
    ## Note: Picard sorting takes both SAM *and* BAM files, and outputs to BAM. So from here on out, we will be operating on BAM files only.
    my $sortCmd = (qq{java -Xmx${GIGABYTES_FOR_PICARD}g -jar ${SORTSAM_PATH} }
		   . qq{ INPUT=${latest} } ## Picard SortSam.jar accepts both SAM and BAM files as input!
		   . qq{ SORT_ORDER=coordinate }
		   . qq{ OUTPUT=${sortedFile} });
    
    my $countAllReadsCmd = qq{${SAMTOOLS_PATH} view ${latest} | wc -l };
    
    my $filterMappedOnlyCmd = (qq{samtools view -h -F 4 $latest } ## <-- only include the mapped (mapping flag is "4" apparently) reads
			       . qq{ | samtools view -bS - } ## <-- Convert back to BAM
			       . qq{ > $filteredFile });   ## <-- output file location
    my $maxFileHandles = int(0.8 * $ULIMIT_RESULT); ## This is for the MAX_FILE_HANDLES_FOR_READ_ENDS_MAP parameter for MarkDuplicates: From the Picard docs: "Maximum number of file handles to keep open when spilling read ends to disk. Set this number a little lower than the per-process maximum number of file that may be open. This number can be found by executing the 'ulimit -n' command on a Unix system. Default value: 8000."
    my $dedupCmd = (qq{java -Xmx${GIGABYTES_FOR_PICARD}g -jar ${MARKDUPLICATES_PATH} }
		    . qq{ INPUT=${latest} } ## Picard SortSam.jar accepts both SAM and BAM files as input!
		    . qq{ REMOVE_DUPLICATES=TRUE }
		    . qq{ MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=$maxFileHandles }
		    . qq{ OPTICAL_DUPLICATE_PIXEL_DISTANCE=100 } ## <-- 100 is default
		    . qq{ METRICS_FILE=${dedupExtraMetricsFile} }
		    . qq{ OUTPUT=${dedupFile} });



    sub updateSymlinkToPointTo($) {
	## We always want each subsequent command to operate on the MOST RECENT file.
	## Sometimes, the user doesn't want to run every step in the pipeline. So steps are skipped.
	## We use this symlink to keep track of which file to work on.
	## In retrospect, we could probably do this with a variable instead, it would probably be safer than using the filesystem.
	my ($target) = @_;
	my $RM_CMD = "/bin/rm --preserve-root --force $latest";
	alexSystemCall("${RM_CMD} && ln -s ${target} $latest");
    }
    
    ### Now to actually RUN the various things ###
    if (!$shouldSort) {
	if (!$isBam) {
	    die "If you specified to NOT sort, then you have to input a BAM file. SAM files are not allowed when there is no sorting. Sorry. We could fix this by adding some kind of samtools convert-to-bam here.\n";
	}
    }

    updateSymlinkToPointTo($file);
    
    if ($shouldSort) {
	alexSystemCall($sortCmd);	
	updateSymlinkToPointTo($sortedFile);
    }

    my $numReadsBeforeWeFiddledWithTheFile = "UNDEFINED";
    if ($shouldCalculateSummary) { 
	# Count the number of reads in the file BEFORE we filter
	# but AFTER we run the sorting command, to make sure the file is a BAM file.
	$numReadsBeforeWeFiddledWithTheFile = `$countAllReadsCmd`;
    }

    if ($shouldFilterMappedOnly) {
	alexSystemCall($filterMappedOnlyCmd);
	updateSymlinkToPointTo($filteredFile);
    }

    if ($shouldRemoveDupes) {
	alexSystemCall($dedupCmd);
	updateSymlinkToPointTo($dedupFile);
    }

    my $fileToMakeIndexFrom = readlink "$latest";
    my $prefix = $file;
    $prefix =~ s/.[bs]am$//i; ## remove the original .bam or .sam prefix!
    my $finalBAM = "${prefix}.processed.bam";
    my $finalIndex = "${prefix}.processed.bai";

    if ($shouldGenerateIndex) {
	## This should always be the LAST step
	alexSystemCall(qq{samtools index ${fileToMakeIndexFrom} ${finalIndex} });
	alexSystemCall(qq{ mv -f $fileToMakeIndexFrom $finalBAM });
    }

    if ($shouldCalculateSummary) { 
	open FILE, ">", $summaryStatsFile or die $!; { ## Note: OVERWRITING this summary file.
	    print FILE "Number of reads found in <$file> before any filtering: " . $numReadsBeforeWeFiddledWithTheFile . "\n"; 
	    print FILE "\n";
	} close(FILE);
    }

    push(@outFiles, $finalBAM);
    if ($shouldGenerateIndex) {
	datePrint("[Done] with <$file>. Generated the output files <${finalBAM}> and <${finalIndex}>\n\n");
	push(@outFiles, $finalIndex);
    } else {
	datePrint("[Done] with <$file>. Generated the output file <${finalBAM}>\n\n");
    }
    $numFilesSuccessfullyProcessed++;
}

datePrint("*************************************************************\n");
datePrint("Successfully generated the following files:\n");
foreach (@outFiles) {
    datePrint("  - $_\n");
}
datePrint("[DONE]\n\n");

__DATA__

rnaseq_filter_agw.pl --fasta=<genome_fasta_file> <INPUT BAM / SAM FILE>

Processes a SAM or BAM alignment file to generate SORTED, INDEXED, AND FILTERED BAM files.

The output is always a BAM file, regardless of the input.

OPTIONS:
   --keep or --keepall: Do NOT remove any reads. Just sort & index.

   --dryrun : Print what *would* be done, but do not actually perform any filtering.

FINE-TUNING OPTIONS:
   --noindex : Do not generate a bam index (Default: generate a ".bai" index file.)

   --nosort : Do not sort the input file (Default: output is sorted.)

   --nomapfilter: Do not remove the unmapped reads (Default: unmapped reads are removed.)



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

#!/usr/bin/perl

## rnaseq_filter_agw.pl
## by Alex Williams, May 2012

## This script handles a few basic processing steps for an aligned SAM or BAM file.

## =========================================================================
## DETAILED README IS AT THE BOTTOM OF THIS FILE --> SCROLL DOWN to __DATA__
## (Or run this script with "--help" to see the full description.)
## =========================================================================

use strict;  use warnings;  use English '-no_match_vars';
use Getopt::Long; use File::Basename;
my $VERSION = "0.1.1"; ## Which version of the program is this. Change this if you update the program significantly.

my $SUCCESS_STATUS = 0; ## <-- this should be ZERO

my $isDryRun = 0;
my $shouldSort = 1;
my $shouldFilterMappedOnly = 1;
my $shouldRemoveDupes = 1;
my $shouldQC = 1;
#my $genomeFastaFile = undef;
my $keepAllReads = 0;

my $SAMTOOLS_PATH = `which samtools`;  chomp($SAMTOOLS_PATH); ## Tries to find "samtools" on the user's $PATH
my $SORTSAM_PATH   = `which SortSam.jar`;  chomp($SORTSAM_PATH); ## Tries to find "SortSam.jar" from the Picard suite in the user's $PATH
my $GIGABYTES_FOR_PICARD = 8;
my $MARKDUPLICATES_PATH = `which MarkDuplicates.jar`; chomp($MARKDUPLICATES_PATH); ## Tries to find "MarkDuplicates.jar" from the Picard suite in the user's $PATH
my $ULIMIT_RESULT = 1024; ## result of running the shell command ulimit -n. Since this is a shell built-in, it can, for some reason, not be run like a real command, so backticks don't work. `ulimit -n`; chomp($ULIMIT_RESULT);

my @successfullyMadeFiles = (); # No files generated yet.
my @failedFiles = (); # No files generated yet.

sub printUsageAndQuit() { ## Function for printing the "you used the wrong arguments to this program" error message
    print STDOUT <DATA>;
    exit(0);
}

sub printVersionAndQuit() { print STDOUT "rnaseq_filter_agw.pl: Version $VERSION\n"; exit(0); }

sub datePrint($) { my $d = `date`; chomp($d); print $d . ":\t" . $_[0]; } ## Timestamping function that is like "print" except with a timestamp at the beginning of the line. Use it instead of "print."

sub alexSystemCall($) {
    my ($cmd) = @_;
    if ($isDryRun) {
	datePrint("Dry run of command: $cmd\n");
    } else {
	datePrint("Executing command: $cmd\n");
    }
    if (!$isDryRun) { return(system($cmd)); }
    else { return $SUCCESS_STATUS; } ## assume that
}

sub appendToSummaryFile($$) {
    my ($whatToWrite, $filename) = @_;
    open FILE, ">>", $filename or die $!; { ## Note: APPENDING TO this summary file.
	print FILE $whatToWrite;
	print STDERR $whatToWrite;
    } close(FILE);
}

sub reportCommandFailure($$) {
    my ($cmd, $fileThatFailed) = @_;
    datePrint("[ERROR]: In the processing for <$fileThatFailed>, this command failed: $cmd\n" . "We are skipping to the next file now.\n");
    push(@failedFiles, $fileThatFailed);
}

sub verifyToolPathsOrDie() {
    # Below: Check for the location of some required binaries.
    # We depend on samtools & the Picard suite being installed already.
    if (length($SAMTOOLS_PATH) <= 1) { die "Could not find <samtools> in the \$PATH. Make sure the <samtools> executable is somewhere in your \$PATH. You should be able to type \"samtools\" and run the command; if that is not the case, then this program won't run either.\nYou should be able to install samtools with the command `sudo apt-get samtools`.\n"; }
    if (length($SORTSAM_PATH) <= 1) {  die "Could not find <SortSam.jar> (part of the Picard suite) in the \$PATH. Make sure the <SortSam.jar> java applet is somewhere in your \$PATH. You should be able to type \"SortSam.jar\" and get a weird java error message; if that is not the case, then this program won't run either. For installation instructions, check http://picard.sourceforge.net/.\nNote that SortSam.jar is PROBABLY not in your default path; try running `locate -i SortSam.jar` to find it. Then add it to somewhere that's in your path. \n"; }
    if (length($MARKDUPLICATES_PATH) <= 1) {  die "Could not find <Markduplicates.jar> (part of the Picard suite) in the \$PATH. Make sure the <MarkDuplicates.jar> java applet is somewhere in your \$PATH. You should be able to type \"MarkDuplicates.jar\" and get a weird java error message; if that is not the case, then this program won't run either. For installation instructions, check http://picard.sourceforge.net/.\nNote that MarkDuplicates.jar is PROBABLY not in your default path; try running `locate -i MarkDuplicates.jar` to find it. Then add it to somewhere that's in your path.\n"; }
}

sub generateIndex($$) {
    my ($bamFile, $indexName) = @_;
    my $indexCmd = qq{samtools index ${bamFile} ${indexName} }; ## Generates the file $indexName
    my $result = alexSystemCall($indexCmd);
    my $COMMAND_SUCCESSFUL_STATUS = 0;
    if ($COMMAND_SUCCESSFUL_STATUS != $result) {
	reportCommandFailure($indexCmd, $bamFile);
    }
    return ($result == $COMMAND_SUCCESSFUL_STATUS); ## Returns 1 if the command ran with no problems, 0 otherwise.
}

GetOptions("help|?|man"        => sub { printUsageAndQuit(); }
	   , "sort!" => \$shouldSort ## specify "--nosort" to avoid sorting
	   , "mapfilter!" => \$shouldFilterMappedOnly ## specify "--nomapfilter" to avoid filtering out the unmappable reads
	   , "keepdupes" => sub { $shouldRemoveDupes = 0; }
#	   , "f|fasta" => \$genomeFastaFile ## used for generating the index
	   , "qc!" => \$shouldQC
	   , "keep|keepall!" => \$keepAllReads
	   , "v|V|version|Version" => sub { printVersionAndQuit(); }
	   , "dry|dryrun|dryRun|dry_run|dry-run" => \$isDryRun
    ) or printUsageAndQuit();

if (scalar(@ARGV) < 1) { die "ARGUMENT ERROR: This script requires at least one input filename---you need to give it a BAM or SAM file to operate on.\nRun it with --help in order to see the options.\n"; }

print STDOUT "The following files will be processed by rnaseq_filter_agw:\n";
foreach (@ARGV) {
    print STDOUT "  - $_\n";
}

if ($keepAllReads) {
    print STDOUT "We are running in the MOST CONSERVATIVE mode, where we DO NOT REMOVE ANY READS.\n";
    print STDOUT "   * Unmapped reads will be RETAINED.\n";
    print STDOUT "   * Duplicate-location reads will be RETAINED.\n";
    print STDOUT "If you specified 'nomapfilter', that option is now OVERRIDDEN.\n";
    $shouldRemoveDupes = 0;
    $shouldFilterMappedOnly = 0;
}

verifyToolPathsOrDie();

my $numFilesSuccessfullyProcessed = 0;
my $numFilesNotOK = 0;

foreach my $file (@ARGV) {
    my $prefix = $file; $prefix =~ s/.[bs]am$//i; ## remove the original .bam or .sam prefix!
    my $latest = undef; ## keep track of which file we should be operating on. This is important in case we skip steps, which can happen if the user specifies certain flags to this program
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
	next; # Go to the next file...
    }

    datePrint("Now processing the file <$file>...\n");
    my $filteredFile          = "$prefix.filtered.tmp.bam";
    my $sortedFile            = "$prefix.sorted.tmp.bam";
    my $dedupFile             = "$prefix.no_duplicates.tmp.bam";
    my $dedupExtraMetricsFile = "$prefix.no_duplicates.extra.txt";
    my $summaryStatsFile      = "$prefix.summary.stats.txt";
    
    open FILE, ">", $summaryStatsFile or die $!; { ## <-- Note: NOT APPENDING -- clearing out the file completely!!!
	print FILE "Summary file for <$file>:\n";
    } close(FILE);
    
    $latest = $file; ## "$latest" is a variable that holds the filename of the MOST RECENT file. It needs to constantly be updated as we change the filenames after processing files. To start, it's the base $file name.
    
    ### Now to actually RUN the various things ###
    if (!$isBam) {
	appendToSummaryFile("Converting SAM --> BAM for file <$latest>...\n", $summaryStatsFile);
	my $bamName = $latest; $bamName =~ s/[.]sam$/.bam/;
	my $convertToBamCmd = qq{samtools view -bS ${latest} > ${bamName} };

	if (-e $bamName && (-s $bamName > 0)) { 
	    appendToSummaryFile("The BAM file version ($bamName) of the sam file ($latest) that you specified on the command line ALREADY EXISTS. We will not be re-generating it, but will instead be using the same-named BAM file ($bamName) as input. You must MANUALLY REMOVE $bamName if you want to re-generate it.\n", $summaryStatsFile);
	    print STDERR 
	} else {
	    my $result = alexSystemCall($convertToBamCmd);
	    if ($result != $SUCCESS_STATUS) { reportCommandFailure($convertToBamCmd, $latest); next; } ## go to the next file...
	}
	$latest = $bamName;
	appendToSummaryFile("Converted SAM --> BAM.\n  Command was: $convertToBamCmd\n", $summaryStatsFile);
    }


    if (!$shouldSort) {
	## If the user said NOT to sort the files, then let's double-check the header to make sure it says 'SO:coordinate' or else MarkDuplicates.jar
	## will crash stupidly
	my $sortOrder = `samtools view -H ${latest} | grep '\@HD.*SO:'`; chomp($sortOrder);
	if (not($sortOrder =~ m/SO:coordinate/)) {
	    print STDERR "ERROR IN <$latest>: Input file warning: When --nosort is specified, the sort order flag in the BAM file MUST be 'coordinate' (SO:coordinate).\nAnything else---including just plain 'sorted'---will cause 'MarkDuplicates.jar' to crash (this is a bug in MarkDuplicates.jar that CANNOT be solved even by forcing ASSUME_SORTED=true)! Therefore, you must run this program again WITHOUT the --nosort directive!\n\nIn your specific case, the line in the BAM header of <$latest> that SHOULD say 'SO:coordinate' actually says the following:\n  This --> $sortOrder\n\n";
	    die "Critical problem in input files; don't specify '--nosort' on these files, as the header is screwed up! Run this command again WITHOUT --nosort!\n";
	}
    }
    
    my $numOriginalReads = `samtools view -c $latest`; chomp($numOriginalReads);
    appendToSummaryFile(("Number of reads in <$latest> before filtering:\t" . $numOriginalReads . "\n"), $summaryStatsFile);
    
    ## Note: Picard sorting takes both SAM *and* BAM files, and outputs to BAM. So from here on out, we will be operating on BAM files only.
    if ($shouldSort) {
	appendToSummaryFile("Sorting <$latest>...\n", $summaryStatsFile);
	my $sortCmd = (qq{java -Xmx${GIGABYTES_FOR_PICARD}g -jar ${SORTSAM_PATH} }
		       . qq{ INPUT=${latest}  SORT_ORDER=coordinate  OUTPUT=${sortedFile} });
	if (alexSystemCall($sortCmd) != $SUCCESS_STATUS) { reportCommandFailure($sortCmd, $latest); datePrint("DEBUGGING MESSAGE FROM ALEX: Maybe the input file <$latest> didn't have a *header* line? The header is REQUIRED for sorting---check your input file ($latest) and make sure it has header lines. If it doesn't, then you'll need to re-header the file with samtools (`samtools reheader <in.header.sam> <in.bam>`). (That command assumes that you already have a properly-headered SAM file somewhere.)\n"); next; } ## go to the next file...
	$latest = $sortedFile; ## update the "latest" file to be this newly-sorted one
	appendToSummaryFile("Sorted <$latest> with Picard SortSam.jar.\n  Command was: $sortCmd\n", $summaryStatsFile);
    }
    
    ## COUNTING TIME!
    ## lower case -f means ONLY INCLUDE THESE READS. upper case -F means EXCLUDE THESE READS.
    my $numFailedToMap       = `samtools view -c -f 0x004 $latest`; chomp($numFailedToMap);
    my $numNonPrimary        = `samtools view -c -f 0x100 $latest`; chomp($numNonPrimary);
    my $numFailedBasicQC     = `samtools view -c -f 0x200 $latest`; chomp($numFailedBasicQC);
    my $numOpticalDuplicates = `samtools view -c -f 0x400 $latest`; chomp($numOpticalDuplicates);
    my $numMateDoesNotMap    = `samtools view -c -f 0x008 $latest`; chomp($numMateDoesNotMap); ## <-- this one only occurs in paired-end reads

    appendToSummaryFile(('-'x100 . "\n"), $summaryStatsFile);
    appendToSummaryFile(("  > Original number of reads BEFORE any filtering                 :\t" . $numOriginalReads . "\n"), $summaryStatsFile);
    appendToSummaryFile(("  * Number of reads that failed to map              (SAM flag 0x4):\t" . $numFailedToMap . "\n"), $summaryStatsFile);
    appendToSummaryFile(("  * Number of non-primary reads                   (SAM flag 0x100):\t" . $numNonPrimary . "\n"), $summaryStatsFile);
    appendToSummaryFile(("  * Number of reads that failed basic QC          (SAM flag 0x200):\t" . $numFailedBasicQC . "\n"), $summaryStatsFile);
    appendToSummaryFile(("  * Number of reads that were optical duplicates  (SAM flag 0x400):\t" . $numOpticalDuplicates . "\n"), $summaryStatsFile);
    appendToSummaryFile(("  * Number of reads where the mate pair did not map (SAM flag 0x8):\t" . $numMateDoesNotMap . "\n"), $summaryStatsFile);

    my $totalReadsFilteredBySamtools = ($numFailedToMap + $numNonPrimary + $numFailedBasicQC + $numOpticalDuplicates + $numMateDoesNotMap);
    appendToSummaryFile(('-'x100 . "\n"), $summaryStatsFile);
    appendToSummaryFile(("Total number of reads that <samtools -F> removes:\t" . $totalReadsFilteredBySamtools . "\n"), $summaryStatsFile);
    appendToSummaryFile(("              Out of the initial number of reads:\t" . $numOriginalReads . "\n"), $summaryStatsFile);
    appendToSummaryFile(("           Percentage of reads that were removed:\t" . sprintf("%.2f", ($totalReadsFilteredBySamtools/$numOriginalReads*100)) . "\n"), $summaryStatsFile);

    appendToSummaryFile("\n\n", $summaryStatsFile);

    ## -F means "skip flag"  (-f means "require flag")
    ## -F 0x0004 means "skip flag ID #4", which is "the query sequence itself is unmapped"
    ## -F 0x0008 means skip pairs where the MATE is unmapped . We are NOT skipping these for now.
    ## -F 0x0100 means "skip the 'not primary alignment' flag" (No secondary alignments)
    ## -F 0x0200 means "skip reads that fail vendor/platform quality control"
    ## -F 0x0400 means "skip optical / PCR duplicates"
    ## Note that the "-u" flag means "send the data to the next step in UNCOMPRESSED form."
    ## The "-u" flag should be used for all these cases EXCEPT for the very last one, where we want the data compressed again.
    my $filterMappedOnlyCmd = (qq{     samtools view -hbu -F 0x004  $latest } ## <-- only include the mapped (mapping flag is "4" apparently) reads, only include PRIMARY alignments, and skip any failing-QC reads
			       . qq{ | samtools view -hbu -F 0x100 - } ## remove "alignment is not primary"
			       . qq{ | samtools view -hbu -F 0x200 - } ## <-- remove "read fails QC"
			       . qq{ | samtools view -hbu -F 0x400 - } ## <-- remove "read is PCR or optical duplicate"
			       . qq{ | samtools view -hb  -F 0x008 - } ## <-- remove "mate pair did not map" reads
			       . qq{ > $filteredFile });   ## <-- output file location
    if ($shouldFilterMappedOnly) {
	my $result = alexSystemCall($filterMappedOnlyCmd);
	if ($result != $SUCCESS_STATUS) { reportCommandFailure($filterMappedOnlyCmd, $latest); next; } ## go to the next file...
	
	$latest = $filteredFile;
	appendToSummaryFile(qq{samtools: removed reads with certain flags.\n  Command was: $filterMappedOnlyCmd\n\n}, $summaryStatsFile);
	appendToSummaryFile(qq{  * Removed reads with flag 0x4 (\"sequence is unmapped\")\n}, $summaryStatsFile);
	appendToSummaryFile(qq{  * Removed reads with flag 0x100 (\"alignment is not primary\")\n}, $summaryStatsFile);
	appendToSummaryFile(qq{  * Removed reads with flag 0x200 (\"read fails platform/vendor quality checks\")\n}, $summaryStatsFile);
	appendToSummaryFile(qq{  * Removed reads with flag 0x400 (\"read is either a PCR or an optical duplicate\")\n}, $summaryStatsFile);
	appendToSummaryFile(qq{  * Removed reads with flag 0x8 (\"in paired-end data, the mate pair is unmapped\")\n}, $summaryStatsFile);
    }
    
    appendToSummaryFile("\n\n", $summaryStatsFile);
    
    my $maxFileHandles = int(0.8 * $ULIMIT_RESULT); ## This is for the MAX_FILE_HANDLES_FOR_READ_ENDS_MAP parameter for MarkDuplicates: From the Picard docs: "Maximum number of file handles to keep open when spilling read ends to disk. Set this number a little lower than the per-process maximum number of file that may be open. This number can be found by executing the 'ulimit -n' command on a Unix system. Default value: 8000."
    my $dedupCmd = (qq{java -Xmx${GIGABYTES_FOR_PICARD}g -jar ${MARKDUPLICATES_PATH} }
		    . qq{ INPUT=$latest } ## Picard SortSam.jar accepts both SAM and BAM files as input!
		    . qq{ REMOVE_DUPLICATES=TRUE }
		    #. ((!$shouldSort) ? qq{ ASSUME_SORTED=TRUE } : qq { }) ## <-- if we use --nosort, then ASSUME they are sorted no matter what!
		    . qq{ MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=$maxFileHandles }
		    . qq{ OPTICAL_DUPLICATE_PIXEL_DISTANCE=100 } ## <-- 100 is default
		    . qq{ METRICS_FILE=$dedupExtraMetricsFile }
		    . qq{ OUTPUT=$dedupFile }
	);
    if ($shouldRemoveDupes) {
	my $result = alexSystemCall($dedupCmd);
	if ($result != $SUCCESS_STATUS) { reportCommandFailure($dedupCmd, $latest); next; } ## go to the next file...
	$latest = $dedupFile;
	appendToSummaryFile("Removed duplicate reads.\n  Command was: $dedupCmd\n", $summaryStatsFile);
	
	appendToSummaryFile(qq{\n\n}, $summaryStatsFile);
	appendToSummaryFile(qq{MarkDuplicates.jar status report:\n}, $summaryStatsFile);
	my $appendMarkDuplicatesMetricsCmd = (qq{grep -A 2 'METRICS CLASS' $dedupExtraMetricsFile >> $summaryStatsFile }
					      . qq{ && /bin/rm $dedupExtraMetricsFile});
	if (alexSystemCall($appendMarkDuplicatesMetricsCmd) != $SUCCESS_STATUS) { reportCommandFailure($dedupCmd, $latest); next; } ## go to the next file...
	appendToSummaryFile(qq{\n\n}, $summaryStatsFile);
    }

    my $numReadsAfterProcessing = `samtools view -c $latest`; chomp($numReadsAfterProcessing);
    my $percentRemaining = sprintf("%.2f", $numReadsAfterProcessing/$numOriginalReads*100);

    appendToSummaryFile(("Statistics about <$latest> after all filtering steps:\n"), $summaryStatsFile);
    appendToSummaryFile(("  *      Number of reads remaining:\t" . $numReadsAfterProcessing . "\n"), $summaryStatsFile);
    appendToSummaryFile(("  *  ...out of this initial number:\t" . $numOriginalReads . "\n"), $summaryStatsFile);
    appendToSummaryFile(("  *              % reads remaining:\t" . $percentRemaining . "\n"), $summaryStatsFile);
    appendToSummaryFile(("  *               ...and % removed:\t" . (100.0 - $percentRemaining) . "\n"), $summaryStatsFile);
    
    my $finalBAM = "${prefix}.processed.bam";
    my $finalIndex = "${prefix}.processed.bai";
    
    ## Index generation should always be the LAST processing step
    my $indexWasGenerated = generateIndex($latest, $finalIndex);
    if (!$indexWasGenerated) { next; } # skip to the next file if we failed this step.
    
    alexSystemCall(qq{ mv -f ${latest} $finalBAM }); ## Move the file to the done-processing-it location.
    
    push(@successfullyMadeFiles, $finalBAM);
    datePrint("[Done] with <$file>. Generated the output files <${finalBAM}> and <${finalIndex}>\n\n");
    push(@successfullyMadeFiles, $finalIndex);
    $numFilesSuccessfullyProcessed++;
}

if (scalar(@successfullyMadeFiles) > 0) {
    datePrint("*************************************************************\n");
    datePrint("Successfully generated the following files:\n");
    foreach (@successfullyMadeFiles) {
	datePrint("  - $_\n");
    }
}

if (scalar(@failedFiles) > 0) {
    ## If any files FAILED to be processed, then let's report them at the very end here.
    datePrint("*************************************************************\n");
    datePrint("FAILED TO PROCESS THE FOLLOWING FILES (see console log for more details:\n");
    foreach (@failedFiles) {
	datePrint("FAILED:  - $_\n");
    }
}

datePrint("[DONE]\n\n");

__DATA__

rnaseq_filter_agw.pl   [any number of bam/sam files may be specified here]

Processes input SAM/BAM alignment files, and:
   1. Removes any UNMAPPED reads (unless you say --keepall)
   2. SORTS the files in COORDINATE order (unless you say --nosort)
   3. Generates an output BAM file and a BAI index file.

The output is always a BAM file, regardless of the input.

Example:   rnaseq_filter_agw.pl  test.bam samtest2.sam test3.bam

Output: For the example above, the output filenames would be:
   BAM: "test.processed.bam" and "samtest2.processed.bam" and "test3.processed.bam"
   BAI: "test.processed.bai" and "samtest2.processed.bai" and "test3.processed.bai"

Note that this is a DESTRUCTIVE operation that removes reads:
   * If you want to keep EVERY read, specify "--keepall" on the command line.
   * Otherwise, unmapped reads are removed and duplicates are collapsed down to one.

OPTIONS:
   --keep or --keepall: Do NOT remove any reads. Just sort & index.
   --dryrun:            Print what *would* be done, without doing it.

FINE-TUNING OPTIONS:
   --nosort:      Do not sort the input file (Default: output is sorted)
   --nomapfilter: Do not remove the unmapped reads (Default: unmappable reads are removed)
   --noqc:        Do not run the QC metrics. Runs faster.

Requires the following additional programs to already be installed AND in your $PATH:
    * These programs must be in your $PATH (even the Picard JARs must be on the path!)
    * <samtools> must be installed (this is used to remove the unmapped reads)
        - You can install it with `apt-get install samtools`
    * <Picard> must be installed in order to properly sort the sam/bam files.
        - Specifically, we use SortSam.jar and MarkDuplicates.jar.
        - (samtools is not ideal because it does not set the sort flag at the top of the file)
        - Check at http://picard.sourceforge.net/
    * If you cannot find SortSam.jar or MarkDuplicates.jar, try "locate -i SortSam.jar" .


#!/usr/bin/perl

## This script converts a BAM or SAM file into the required parts for the genome browser.

## by Alex Williams, Feb. 2011

## Convert SAM or BAM to Genome Browser format

use strict;  use warnings;  use diagnostics;
use List::Util qw/max min/;
use Getopt::Long;
use File::Basename;

sub printUsageAndQuit() {
    print STDOUT <DATA>;
    exit(0);
}

sub datePrint($) {
    my $d = `date`;
    chomp($d);
    print $d . ":\t" . $_[0];
}

sub cleanedUpFilename($) {
    my ($name) = @_;
    $name =~ s/[.]bam$//i; $name =~ s/[.]bw$//i; ## remove recognized file extensions
    $name =~ s/[-._]accepted_hits//i; $name =~ s/[-._]sort[e]?[d]?//i; $name =~ s/Browser[-._]//i;
    return($name);
}

my %groupsColHash = (); # Key = group name. Value = color string.

sub guessColorFromFilename($) {
    my ($filename) = @_;
    # Hopefully guesses a reasonable color based on the filename(s). Tries to cluster samples into groups that are the same color.
    my @colors = (  "255,0,0" # red
		    , "0,0,255" # blue
		    , "0,128,0" # green
		    , "128,0,128" # purple
		    , "128,128,128" # gray
		    , "255,140,0" # orange
		    , "0,0,0"     # black
		    , "128,128,0" # olive green
		    , "0,255,0"   # lime green
		    , "0,190,255" # sky blue
		    , "255,0,255" # magenta pink
		    , "0,255,255" # cyan blue
		    , "0,128,128" # teal blue-green
		    , "0,0,128"   # navy blue
		    , "230,220,0" # dandelion yellow
		    , "128,0,0"   # dark red
		    , "100,100,150"); # dark blue-gray
    my $thisGroup = cleanedUpFilename($filename);
    $thisGroup = s/[-._:]*//i; # Delete everything after the first hyphen/dot/underscore/other spacer character. We assume whatever was at the BEGINNING was the group name.
    ($thisGroup ne '') or print STDERR "[WARNING] We could not auto-guess the group name for filename '$filename'. Normally, we assume the grouop is whatever appears BEFORE the first period/hyphen/underscore, and the rest is a replicate number or other ignore-able ID. As a result, colors might not be auto-set by groups!\n";
    if (!exists($groupsColHash{$thisGroup})) { # Haven't seen this group yet. Assign a color index to this group!
	my $numGroupsSeenBeforeThis = scalar(keys(%groupsColHash));
	$groupsColHash{$thisGroup} = $colors[ $numGroupsSeenBeforeThis % scalar(@colors) ]; # Get a color for this group!
    }
    my $thisColor = $groupsColHash{$thisGroup};
    return($thisColor);
}

sub browserTrackString($$) {
    ## Takes "bigWig" or "bam" as its first argument
    ## and a filename as its second argument.
    ## Output: a NEW string that is what the genome browser track description would look like for this file.
    my ($type, $filename) = @_; ## input arguments: type must be bigWig or bam
    if ($type ne "bigWig" and $type ne "bam") { die "Sorry, we only know how to write BIGWIG and BAM files. Problem!\n"; }
    my $redColor = "255,0,0";  ## R, G, B, from 0 to 255
    my $blueColor = "0,0,255"; ## R, G, B, from 0 to 255
    my $color;
    my $parenthetical = "NA"; ## <-- a user-visible additional comment that goes at the end of the track name
    if ($type eq "bigWig") {
	my $colorByGroup = guessColorFromFilename($filename);
	$color = qq{color="$colorByGroup"}; $parenthetical = qq{};
    }
    if ($type eq "bam") { $color = qq{colorByStrand="$redColor $blueColor"}; $parenthetical = qq{ (reads)}; }
    
    my $visStatus = "full";
    if ($type eq "bam") { $visStatus = qq{dense}; } ## bam tracks are DENSE (almost totally hidden) by default
    if ($type eq "bigWig") { $visStatus = qq{full}; } ## bigwig tracks are FULL (shown) by default
    my $cleanedUpName = cleanedUpName($filename);
    my $url = "https://gb.ucsf.edu/bio/browser/YOURLOCATION/${filename}";
    my $str = qq{track type=${type} name="${cleanedUpName}${parenthetical}" description="${cleanedUpName}${parenthetical}" bigDataUrl="${url}" visibility=${visStatus} ${color}\n};
    return($str);
}

my $makeWig = 1; ## By default, generate a bigwig track too. Specify "nowig" to avoid this.
my $shouldSort = 1; ## By default, assume the input SAM/BAM file will still require sorting.
my $shouldScale = 1; ## By default, scale the wiggle tracks to a common height. If we SHOULD scale, then we divide all the values by the (number of reads in the BIGGEST file / number of reads in THIS file)

my $shouldKeepTempFiles = 0;

# ==============================================================================
GetOptions("help|?|man"        => sub { printUsageAndQuit(); }
	   , "wig!" => \$makeWig ## specify "--nowig" to avoid making a wiggle track
	   , "sort!" => \$shouldSort ## "--nosort" avoids the slow sorting step
	   , "scale!" => \$shouldScale ## "--noscale" doesn't scale the output wiggle files. BAM files are never scaled. Default: DO scale
	   , "keeptemp!" => \$shouldKeepTempFiles
    ) or printUsageAndQuit();

if (scalar(@ARGV) == 0) { die "ARGUMENT ERROR: This script requires ONE OR MORE SAM/BAM filenames as arguments!\nIt cannot read from STDIN---sorry!\nExample: convert_SAM_or_BAM_for_Genome_Browser.pl  mySamFile.sam\n[Quitting now.]\n"; }

my @INPUT_FILES = @ARGV;
chomp(@INPUT_FILES);
# ==============================================================================

my $random_number = rand();
my $chrLenTempFile = "Browser.tmp.TEMPORARY.chrlen_" . rand() . "-" . time() . "_convert_SAM_or_BAM_for_Genome_Browser.tmp";

print STDOUT "Files to be processed by convert_SAM_or_BAM_for_Genome_Browser:\n";
foreach (@INPUT_FILES) {
    my $fname = $_;
    print STDOUT "  - ${fname}\n";
    if (!(-e ${fname}) or (!(-s $fname))) { die "Curses, one of the input files ($fname) appears to be missing or zero-length!\n"; }
}

# ==============================================================================
for my $in (@INPUT_FILES) {
    my $isSAM = ($in =~ m/\.sam$/i); # result: a true/false value
    my $isBAM = ($in =~ m/\.bam$/i); # result: a true/false value
    ($isSAM || $isBAM) or die "Alas, the input file <$in> was apparently neither SAM nor BAM, but those are the only file types we support! We can't continue! Fix the input arguments and only include BAM/SAM files. The file names must also end with '.sam' or '.bam'.";
}

# ==============================================================================
my %numReadsHash = ();

if ($shouldScale) {
    # Figure out which bam/sam file is the LARGEST. We will need this data later.
    print "We will be re-scaling all the wiggle tracks to make them comparable between input files with different sequencing depths.\n";
    for my $inFile (@INPUT_FILES) {
	my $isSAM = ($inFile =~ m/\.sam$/i); # result: a true/false value
	my $isBAM = ($inFile =~ m/\.bam$/i); # result: a true/false value
	my $displayCommand = "samtools view" . (($isSAM) ? " -S " : " "); # sam files require the '-S' flag to view with samtools. Do NOT just use 'cat', as it displays the header as well, which we do not want.
	my $numLinesThisFile = `$displayCommand $inFile | wc -l`;
	chomp($numLinesThisFile); # remove the newline!
	$numReadsHash{$inFile} = $numLinesThisFile;
	print " * Number of reads in <$inFile>: $numLinesThisFile\n";
    }
    print "Highest number of reads in a single input file: " . List::Util::max(values(%numReadsHash)) . "\n";
} else {
    print "We will NOT be re-scaling the wiggle tracks. You may want to specify scaling with --scale. Otherwise, beware when comparing wiggle tracks between files---a file with higher sequencing depth will tend to have greater values in the wiggle track, even if a gene had lower expression in reality!\n";
}
# ==============================================================================

if ($shouldSort) {
    print "*"x80 . "\n";
    print "NOTE: We are COORDINATE SORTING each .bam file here.\n";
    print "      If the .bam files are ALREADY sorted by coordinate,\n";
    print "      (which is the case if you used Tophat, but not necessarily with Bowtie)\n";
    print "      you should cancel this command and re-run it with\n";
    print "      --nosort, and it will complete 20 times faster.\n";
    print "*"x80 . "\n";
}

# ==============================================================================

for my $originalInputFilename (@INPUT_FILES) {
    ## Go through each file...
    ## wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bigWigToBedGraph
    # Everything about BAM -> browser is described here http://genome.ucsc.edu/goldenPath/help/bam.html
    my $thingToTypeForBrowser = qq{track type=bam name="Control_EB" color=0,128,255 bigDataUrl="http://lighthouse.ucsf.edu/public_files_no_password/browser_custom_bed/sorted_ctrl.bam"};
    (-e $originalInputFilename) or die "ARGUMENT ERROR: The argument to this program must be a single BAM/SAM file that already exists.\n";
    (-r $originalInputFilename) or die "FILE ERROR: The argument to this program must be a single file that already exists AND is also readable. We couldn't read the file <$originalInputFilename> though!\n";
    (not (-z $originalInputFilename)) or die "FILE ERROR: The file <$originalInputFilename> was 0-length. It's not a real BAM/SAM file!\n";

    my $isFileSamFormat = ($originalInputFilename =~ m/\.sam$/i); # result: a true/false value
    my $isFileBamBinaryFormat = ($originalInputFilename =~ m/\.bam$/i); # result: a true/false value
    my $bamFilename = undef;

    if ($isFileSamFormat) {
	print "<$originalInputFilename> appears to be a SAM file.\n";
	## User specified a SAM file, so we're going to CONVERT it to a BAM file. Samtools must be installed for this to work!
	$bamFilename = $originalInputFilename;
	$bamFilename =~ s/\.sam$/\.bam/i;
	if (not(-e $bamFilename)) {
	    print "Time to make the bam file named <$bamFilename> from the input SAM file that was named <$originalInputFilename>... this will take a minute or two...\n";
	    my $exitCode = system("samtools view -S -b $originalInputFilename > $bamFilename"); ## generate a BAM file from the sam file
	    if (0 != $exitCode) { 
		unlink($bamFilename); # Definitely remove the screwed-up bam file first!! Otherwise we will think this is a valid input file.
		print "ERROR of type SAM99 (exit code $exitCode)! The attempted conversion from SAM -> BAM ($originalInputFilename -> $bamFilename) FAILED, probably because the input SAM file was either: 1) not a real SAM file or 2) did not have a header (which is mandatory).\n";
		my $checkHeaderCmd = "samtools view -S -H $originalInputFilename 2>&1 | tee"; # <-- "tee" redirects even samtools' bizarre output messages to a file! If you try just capturing STDERR by other methods, it is a huge mess.
		my $results = `$checkHeaderCmd`; # this also prints to the screen, by the way!
		if ($results =~ /no .*lines in .*header/) {
		    print "ADDITIONAL INFORMATION: Sure enough, the problem was that there was NO HEADER DATA in your sam file! Specifically, here is the error message from tunning the command <$checkHeaderCmd>:\n$results\n";
		} else {
		    print "ADDITIONAL INFORMATION: It was unclear whether the header was the problem, but you should check the input sam file yourself. View <$originalInputFilename> with 'less' or another tool, and see if you can figure out what is wrong with it.\n$results";
		}
		print "FAILURE: fix the SAM file and try again!\n";
		exit(1);
	    }
	} else {
	    print "Not remaking the BAM file named <$bamFilename>, because it already existed.\n";
	}
    } elsif ($isFileBamBinaryFormat) {
	print "<$originalInputFilename> appears to be a BAM (Binary SAM) file.\n";
	$bamFilename = $originalInputFilename; # that was easy, the user's specified file was ALREADY a bam file
    } else {
	die "INPUT FILE ERROR: input file <$originalInputFilename> must be sequence data in either SAM or BAM format!\n";
    }

    # At this point, we can assume the file is DEFINITELY a bam file--if it was SAM, we have already converted it!
    if ($makeWig) {
	system(qq{samtools view -H $bamFilename | grep "^\@SQ" | sed -e "s/SN://" -e "s/LN://" | cut -f 2,3 > ${chrLenTempFile} });
	## make a temp file with chromosome lengths in it!
    }

    my $bamPrefixWithoutFileExtension = $bamFilename; ## Don't just get the basename: get the WHOLE PATH that led up to this file. basename($bamFilename);
    $bamPrefixWithoutFileExtension =~ s/\.bam$//i;
    $bamPrefixWithoutFileExtension =~ s/[\/:;,]/_/g; ## slashes and ':;,' characters go to underscores

    my $browserTrackDescriptionFile = "Browser.Track.Descriptions.${bamPrefixWithoutFileExtension}.txt";
    my $sortBamFilePrefix           = "Browser.sort.${bamPrefixWithoutFileExtension}";
    my $sortBamFullFilename         = "${sortBamFilePrefix}.bam";
    my $bamIndexOutfile             = "${sortBamFilePrefix}.bam.bai";
    my $bedGraphIntermediateFile    = "Browser.tmp.${bamPrefixWithoutFileExtension}.bedgraph";
    
    if ((-e $sortBamFullFilename) && (-s $sortBamFullFilename) and (-e $bamIndexOutfile) && (-s $bamIndexOutfile)) {
	datePrint(": [Skipping] We are NOT continuing with the sorted-BAM-file generation, because such a file already exists. Remove it if you want to recompute it!\n");
    } else {
	if (!$shouldSort) {
	    datePrint(": Assuming that the file <$bamFilename> is already COORDINATE sorted---we will make a symlink to <$sortBamFullFilename> rather than re-sorting it with Samtools.\n");
	    system("ln -s $bamFilename $sortBamFullFilename");
	} else {
	    my $sortCmd = "samtools sort $bamFilename $sortBamFilePrefix"; ## <-- Note: this should NOT have the .bam extension, because samtools automatically adds the ".bam" extension
	    datePrint(": Running the SAMTOOLS sort-by-leftmost-coordinate command: ${sortCmd}...\n");
	    system($sortCmd);
	}
    }
    
    if ((-e $bamIndexOutfile) && (-s $bamIndexOutfile > 0)) {
	print "[Skipping] We are NOT re-generating an index file. One already exists ($bamIndexOutfile).\n";
    } else {
	my $indexCmd = "samtools index ${sortBamFullFilename}";
	datePrint(": Running the SAMTOOLS index command: ${indexCmd}...\n");
	system($indexCmd);
    }

    my $scaleFactor      = "ERROR_UNDEFINED";
    my $maxReads         = $shouldScale ? List::Util::max(values(%numReadsHash)) : "ERROR_UNDEFINED";
    my $numReadsThisFile = $shouldScale ? $numReadsHash{$originalInputFilename} : "ERROR_UNDEFINED";
    if ($shouldScale) {
	if (!defined($numReadsHash{$originalInputFilename})) { die "Programming error: the numReadsHash was not defined for <$originalInputFilename>! Something is wrong here."; }
	if ($numReadsThisFile == 0) {
	    print("************\nThat's weird, there were NO reads in this file, <$originalInputFilename>. That could indicate something being wrong.\n********\n");
	    $scaleFactor = 1.0; # this is a weird condition, but don't scale anything if there are no reads!
	} else {
	    $scaleFactor = $maxReads / $numReadsThisFile;
	    $scaleFactor = sprintf("%.3f", $scaleFactor); # No point in keeping a zillion degrees of useless decimal precision (which also makes the filesizes larger)
	    # Scales all the values!
	}
	($scaleFactor >= 0.99999) or die "BUG REPORT 87X: how is the scale factor LESS than 1? That should be impossible. This is a programming error! We are supposed to scale the smaller files UP to the largest file, so no files should be scaled down at all!";
    }

    my $scaleString = ($shouldScale) ? ".scaled_by_${scaleFactor}" : ""; # Change the title of the files, if the user specified scaling.

    my $bedGraphScaledFile = "Browser.tmp.${bamPrefixWithoutFileExtension}${scaleString}.bedgraph";
    my $bigWigOutFile      = "Browser.${bamPrefixWithoutFileExtension}${scaleString}.bw";

    if ($makeWig) {
	if ((-e $bigWigOutFile) && (-s $bigWigOutFile > 0)) {
	    datePrint(": [Skipping] We are NOT continuing with the generation of a bigwig file, because the bigwig file $bigWigOutFile already exists. Remove it if you want to recompute it!\n");
	} else {
	    datePrint(": Now generating a 'bigWig' browser wiggle track named $bedGraphIntermediateFile (this is slow, and can take up to an hour per input RNASeq file!)...\n");
	    ## Note: genomeCoverageBed REQUIRES that the input bam file be sorted by position. Luckily, we just sorted it!
	    ## The "-g" genome file needs the chromosome sizes. Luckily, this is information we can find in the $chrLenTempFile that we just made!
	    my $wigCmd1 = (qq{genomeCoverageBed -split -bg -ibam $sortBamFullFilename -g $chrLenTempFile > $bedGraphIntermediateFile});
	    datePrint(": Now running this command: $wigCmd1\n");
	    system($wigCmd1);
	    my $whichFileToWigify = undef;
	    if ($shouldScale) {
		if ($scaleFactor == 1.0) {
		    datePrint(": Skipping scaling the largest bedGraph file by a factor of $scaleFactor, since that would have no effect.\n");
		    # No need to "scale" this file by a factor of 1.0 (since that will change nothing)! Instead we'll just use the original bedGraph file
		    $whichFileToWigify = $bedGraphIntermediateFile; # <-- just use the OLD intermediate file, don't waste time recomputing a new scaled one.
		} else {
		    # Now we will write out the SCALED wiggle track
		    datePrint(": About to scale the bedGraph file by a factor of $scaleFactor ($maxReads / $numReadsThisFile)...\n");
		    open SCALED, "> $bedGraphScaledFile" or die $!;
		    open BWFILE, $bedGraphIntermediateFile or die $!;
		    while (my $line = <BWFILE>) {
			chomp($line);
			my @arr = split(/\t/, $line);
			my $scaledArr3 = ($arr[3]) * $scaleFactor; # scale the column with the actual bedgraph values in it. # 4th column in a bedgraph file, so index (counting from 0) is 3
			print SCALED (join("\t", $arr[0], $arr[1], $arr[2], $scaledArr3) . "\n");
		    }
		    close BWFILE or die $!;
		    close SCALED or die $!;
		    $whichFileToWigify = $bedGraphScaledFile;
		}
	    } else {
		$whichFileToWigify = $bedGraphIntermediateFile; # there isn't a scaled file, so don't look for one!
	    }
	    ## -clip means "allow errant entries off the end of the chromosome". This is important, because otherwise wigToBigWig quits with errors like "something went off the end of chr12_random"
	    my @wigCmd2 = ("wigToBigWig", "-clip", $whichFileToWigify, $chrLenTempFile, $bigWigOutFile);
	    datePrint("Now running this command: @wigCmd2\n");
	    system(@wigCmd2); ## makes a BIGWIG file from the WIG file
	    if (!$shouldKeepTempFiles) {
		if (-e $chrLenTempFile)           { unlink($chrLenTempFile) }; ## <-- delete the $chrLenTempFile now! We don't need it anymore.
		if (-e $bedGraphScaledFile)       { unlink($bedGraphScaledFile); }
		if (-e $bedGraphIntermediateFile) { unlink($bedGraphIntermediateFile); }
	    }
	}
    } else {
	print qq{[Skipping] the generation of a bigWig file, because "--nowig" was specified on the command line.\n};
    }

    open FILE, ">>", $browserTrackDescriptionFile or die $!; ## APPEND TO THE FILE!!!
    print FILE "\n";
    if ($makeWig) { print FILE browserTrackString("bigWig", ${bigWigOutFile}); }
    print FILE browserTrackString("bam", ${sortBamFullFilename});
    print FILE "\n";
    close(FILE);

    datePrint("[DONE]\n\n");
    print "Here are the final output files that you will probably want to put on your server:\n"
	. " - ${sortBamFullFilename}\n"
	. " - ${bamIndexOutfile}\n"
	. " - ${bigWigOutFile}\n";
    print "\n";
    print "1. You will probably want to move these files to your web server now.\n"
	. "   For my server, the command: scp sorted* lighthouse.ucsf.edu:\n"
	. "   copies everything from this machine to my home directory on that machine.\n";
    print "\n";
    print "2. Make sure that the files are READABLE by everyone, and that the directories they\n"
	. "   are in are also readable and executable by everyone. Double-check that you can access\n"
	. "   your files in a regular web browser.\n";
    print "\n";

    print "3. Then you can go to the genome browser gateway page (http://genome.ucsc.edu/cgi-bin/hgGateway)\n"
	. "   and click the \"add custom tracks\" button near the top. Then, in \"Paste URLs or data\", type\n"
	. "   " . $thingToTypeForBrowser . "\n"
	. "   (replacing sorted_ctrl.bam with your actual .bam file---note that the .bai files aren't mentioned here, but MUST be in the same directory.\n";
    print "\n";
    print "4. Note that if you want to save this custom track or send it to others, you will need\n"
	. "   a custom genome browser SESSION---just copying / saving the URL is not sufficient!\n";
    print "\n";
    print "[DONE] The file $originalInputFilename appears to have been processed successfully!\n";
#print "Assuming that your files are hosted on lighthouse, here are the lines to paste into the table browser to make things work:\n";
#print "track type=bam name="Control_EB" color=0,128,255 bigDataUrl=http://lighthouse.ucsf.edu/public_files_no_password/browser_custom_bed/sorted_ctrl.bam"
    print "\n";
}

__DATA__

convert_SAM_or_BAM_for_Genome_Browser.pl  <INPUT BAM / SAM FILES>

Example: convert_SAM_or_BAM_for_Genome_Browser.pl   myFile.bam   otherFile.bam

Processes SAM or BAM alignment files to generate files for the UCSC Genome Browser.
You can now specify multiple BAM input files at once.

This is a script that will generate UCSC-genome-browser-ready BAM files.
It will convert / produce the UCSC-ready files from any SAM or BAM file you want.
The BAM/SAM files can be HUGE (10+ GB), but can be hosted locally. Then you tell the UCSC Genome Browser
where your BAM files are, and it magically uses the BAM files from your local web server, without copying the entire file.

Requres the following additional software:
    * <samtools> must be installed (to sort the bam files)
        - You can install it with `apt-get install samtools`

If you want a wiggle track, which you most likely do (if you do not, you can specify `--nowig`):
    * <genomeCoverageBed> must be installed
    * <wigToBigWig> must be installed

Input to program:
    * Any number of SAM or BAM files (with aligned reads to a reference genome).

These files can be viewed in the genome browser.

OPTIONS:

   --noscale: Do NOT scale. Default is to rescale the wiggle files.
            --scale resizes the wiggle files so that they are on the same scale with
            respect to the total number of reads in each file.
            Does NOT affect the BAM files.
            Example: if you have file A with 100 reads and file B with 50 reads,
            the values in file B will be scaled by (100/50) = 2.0 and
            the values in file A will be scaled by (100/100) = 1.0 (unchanged).
            This way you can compare wiggle tracks between files with different
            sequencing depth.

   --nosort: Assume the input file is ALREADY sorted. Speeds things up.

   --nowig: Add this to NOT generate a wiggle track.
            Generating a wiggle track is SLOW.
            Also, if you specify --nowig, then you no longer need to specify
            a "fai" fasta index file.

   --keeptemp: Keeps the temp files (bedgraph files) instead of deleting them.

EXAMPLES:

convert_SAM_or_BAM_for_Genome_Browser.pl alignment.bam
convert_SAM_or_BAM_for_Genome_Browser.pl --nowig  alignment.bam

Files that are generated from the input YOURFILE.sam:
 1. Browser.sort.YOURFILE.bam (sorted version of the BAM/SAM file) (track type=bam)
 2. Browser.sort.YOURFILE.bam.bai (BAM index file -- required for the browser)
 3. Browser.tmp.YOURFILE.bedgraph (bedgraph intermediate file used to make the bigwig. Can be deleted!)
 4. Browser.YOURFILE.bw (bigWig file for Genome Browser)
 5. Browser.Track.Descriptions.YOURFILE.txt (track descriptions that you paste into the Genome Browser Custom Tracks)

Then you will just need to make a genome browser track description for those two files.

You can find a browser track description in file 5, the browser track description. Or make your own, like this sample one:

track type=bam name="Control_EB" color=0,128,255 bigDataUrl="http://lighthouse.ucsf.edu/public_files_no_password/browser_custom_bed/sorted_ctrl.bam"

  (Note that there is no need to mention the .bam.bai file in the track description,
   *but* it is implicitly used by the Genome Browser and must be present!)

This script just runs the instructions found at:
      http://genome.ucsc.edu/goldenPath/help/bam.html

This script will probably take about 15 minutes to run on a single 5 GB input SAM file.

Generating the bigWig file is SLOW. You can prevent this by specifying --nowig.

This script generates a single 5 kilobyte chromosome length temp file in /tmp/, which it deletes
upon successful script completion.


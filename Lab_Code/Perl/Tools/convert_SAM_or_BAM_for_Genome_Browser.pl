#!/usr/bin/perl

## This script converts a BAM or SAM file into the required parts for the genome browser.

## by Alex Williams, Feb. 2011

## Convert SAM or BAM to Genome Browser format

use strict;  use warnings;  use diagnostics;

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
    if ($type eq "bigWig") { $color = qq{color="$redColor"}; $parenthetical = qq{}; }
    if ($type eq "bam") { $color = qq{colorByStrand="$redColor $blueColor"}; $parenthetical = qq{ (reads)}; }
    
    my $visStatus = "full";
    if ($type eq "bam") { $visStatus = qq{dense}; } ## bam tracks are DENSE (almost totally hidden) by default
    if ($type eq "bigWig") { $visStatus = qq{full}; } ## bigwig tracks are FULL (shown) by default

    my $cleanedUpName = $filename;
    $cleanedUpName =~ s/\.bam$//i; $cleanedUpName =~ s/\.bw$//i; ## remove file extensions
    $cleanedUpName =~ s/.accepted_hits//i; $cleanedUpName =~ s/.sort//i; $cleanedUpName =~ s/Browser.//i;
    
    my $url = "http://lighthouse.ucsf.edu/public_files_no_password/browser_custom_bed/YOUR_FILE_LOCATION/${filename}";
    my $str = qq{track type=${type} name="${cleanedUpName}${parenthetical}" description="${cleanedUpName}${parenthetical}" bigDataUrl="${url}" visibility=${visStatus} ${color}\n};
    return($str);
}



#my $genomeFastaFile = undef;
my $faiFile = undef;  ## <-- "fai" means genome "fasta index (fai)" file
my $makeWig = 1; ## By default, generate a bigwig track too. Specify "nowig" to avoid this.
my $shouldSort = 1; ## By default, assume the input SAM/BAM file will still require sorting.

GetOptions("help|?|man"        => sub { printUsageAndQuit(); }
	   , "wig!" => \$makeWig ## specify "--nowig" to avoid making a wiggle track
	   , "sort!" => \$shouldSort ## "--nosort" avoids the slow sorting step
#	   , "fasta=s" => \$genomeFastaFile
    ) or printUsageAndQuit();

if (scalar(@ARGV) == 0) { die "ARGUMENT ERROR: This script requires ONE OR MORE SAM/BAM filenames as arguments!\nIt cannot read from STDIN---sorry!\nExample: convert_SAM_or_BAM_for_Genome_Browser.pl  mySamFile.sam\n[Quitting now.]\n"; }

my $random_number = rand();
my $chrLenTempFile = "/tmp/chr_length_" . rand() . "-" . time() . "_from_convert_SAM_or_BAM_for_Genome_Browser.tmp";


print STDOUT "File to be processed by convert_SAM_or_BAM_for_Genome_Browser:\n";
foreach (@ARGV) {
    my $fname = $_;
    print STDOUT "  - ${fname}\n";
    if (!(-e ${fname}) or (!(-s $fname))) { die "Curses, one of the input files ($fname) appears to be missing or zero-length!\n"; }
}

for my $originalInputFilename (@ARGV) {
    ## Go through each file...
## wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bigWigToBedGraph

    my $thingToTypeForBrowser = qq{track type=bam name="Control_EB" color=0,128,255 bigDataUrl="http://lighthouse.ucsf.edu/public_files_no_password/browser_custom_bed/sorted_ctrl.bam"};

    unless (-e $originalInputFilename) { die "ARGUMENT ERROR: The argument to this program must be a single BAM/SAM file that already exists.\n"; }
    unless (-r $originalInputFilename) { die "FILE ERROR: The argument to this program must be a single file that already exists AND is also readable. We couldn't read the file <$originalInputFilename> though!\n"; }
    unless (not (-z $originalInputFilename)) { die "FILE ERROR: The file <$originalInputFilename> was 0-length. It's not a real BAM/SAM file!\n"; }

# Everything about BAM -> browser is described here http://genome.ucsc.edu/goldenPath/help/bam.html

    my $fileIsSAM = ($originalInputFilename =~ m/\.sam$/i);
    my $fileIsBAM = ($originalInputFilename =~ m/\.bam$/i);
    my $bamFilename = undef;

    if ($fileIsSAM) {
	print "<$originalInputFilename> appears to be a SAM file.\n";
	## User specified a SAM file, so we're going to CONVERT it to a BAM file. Samtools must be installed for this to work!
	$bamFilename = $originalInputFilename;
	$bamFilename =~ s/\.sam$/\.bam/i;
	if (not(-e $bamFilename)) {
	    print "Time to make the bam file named <$bamFilename> from the input SAM file that was named <$originalInputFilename>... this will take a minute or two...\n";
	    system("samtools view -S -b $originalInputFilename > $bamFilename"); ## generate a BAM file from the sam file
	} else {
	    print "Not remaking the BAM file named <$bamFilename>, because it already existed.\n";
	}
    } elsif ($fileIsBAM) {
	print "<$originalInputFilename> appears to be a BAM (Binary SAM) file.\n";
	$bamFilename = $originalInputFilename; # that was easy, the user's specified file was ALREADY a bam file
    } else {
	die "INPUT FILE ERROR: input file <$originalInputFilename> must be sequence data in either SAM or BAM format!\n";
    }


## At this point, we can assume the file is DEFINITELY a bam file--if it was SAM, we have already converted it!
    if ($makeWig) {
	system(qq{samtools view -H $bamFilename | grep "^\@SQ" | sed -e "s/SN://" -e "s/LN://" | cut -f 2,3 > ${chrLenTempFile} });
	## make a temp file with chromosome lengths in it!
    }

    my $bamPrefixWithoutFileExtension = $bamFilename; ## Don't just get the basename: get the WHOLE PATH that led up to this file. basename($bamFilename);
    $bamPrefixWithoutFileExtension =~ s/\.bam$//i;
    $bamPrefixWithoutFileExtension =~ s/[\/:;,]/_/g; ## slashes and ':;,' characters go to underscores

    my $browserTrackDescriptionFile = "Browser.Track.Descriptions.${bamPrefixWithoutFileExtension}.txt";
    my $sortBamFilePrefix   = "Browser.sort.${bamPrefixWithoutFileExtension}";
    my $sortBamFullFilename = "${sortBamFilePrefix}.bam";
    my $bamIndexOutfile     = "${sortBamFilePrefix}.bam.bai";
    my $bedGraphIntermediateFile = "Browser.tmp.${bamPrefixWithoutFileExtension}.bedgraph";
    my $bigWigOutFile       = "Browser.${bamPrefixWithoutFileExtension}.bw";


    if ((-e $sortBamFullFilename) && (-s $sortBamFullFilename) and (-e $bamIndexOutfile) && (-s $bamIndexOutfile)) {
	print STDERR "[Skipping] We are NOT continuing with the sorted-BAM-file generation, because such a file already exists. Remove it if you want to recompute it!\n"
    } else {
	if (!$shouldSort) {
	    print(">> Assuming that the file <$bamFilename> is already sorted---we will make a symlink to <$sortBamFullFilename> rather than re-sorting it.\n");
	    system("ln -s $bamFilename $sortBamFullFilename");
	} else {
	    my $sortCmd = "samtools sort $bamFilename $sortBamFilePrefix"; ## <-- Note: this should NOT have the .bam extension, because samtools automatically adds the ".bam" extension
	    print(">> Running the SAMTOOLS sort command: $sortCmd\n...\n");
	    system($sortCmd);
	}
}

    if ((-e $bamIndexOutfile) && (-s $bamIndexOutfile > 0)) {
	print STDERR "[Skipping] We are NOT re-generating an index file. One already exists ($bamIndexOutfile).\n";
    } else {
	my $indexCmd = "samtools index ${sortBamFullFilename}";
	print(">> Running the SAMTOOLS index command: $indexCmd\n...\n");
	system($indexCmd);
    }

    if ($makeWig) {
	if ((-e $bigWigOutFile) && (-s $bigWigOutFile > 0)) {
	    print STDERR "[Skipping] We are NOT continuing with the generation of a bigwig file, because the bigwig file $bigWigOutFile already exists. Remove it if you want to recompute it!\n";
	} else {
	    print "Now generating a 'bigWig' browser wiggle track named $bedGraphIntermediateFile (this is slow, and can take up to an hour per input RNASeq file!)...\n";
	    
	    if (-e $bedGraphIntermediateFile && (-s $bedGraphIntermediateFile > 0)) {
		print STDERR "[Skipping] We are NOT creating the wiggle temp file <$bedGraphIntermediateFile>, because it already exists and has non-zero size. Remove it if you want to recompute it.\n";
	    } else {
		## Note: genomeCoverageBed REQUIRES that the input bam file be sorted by position. Luckily, we just sorted it!
		## The "-g" genome file needs the chromosome sizes. Luckily, this is information we can find in the $chrLenTempFile that we just made!
		my $wigCmd1 = (qq{genomeCoverageBed -split -bg -ibam $sortBamFullFilename -g $chrLenTempFile > $bedGraphIntermediateFile});
		datePrint(": Now running this command:\n  $wigCmd1\n"); 
		system($wigCmd1);
	    }
	    ## -clip means "allow errant entries off the end of the chromosome". This is important, because otherwise wigToBigWig quits with errors like "something went off the end of chr12_random"
	    my $wigCmd2 = (qq{wigToBigWig -clip $bedGraphIntermediateFile $chrLenTempFile $bigWigOutFile});
	    datePrint("Now running this command:\n  $wigCmd2\n");
	    system($wigCmd2); ## makes a BIGWIG file from the WIG file
	}
	unlink($chrLenTempFile); ## <-- delete the $chrLenTempFile now! We don't need it anymore.
    } else {
	print STDERR qq{[Skipping] the generation of a bigWig file, because "--nowig" was specified on the command line.\n};
    }

    open FILE, ">>", $browserTrackDescriptionFile or die $!; ## APPEND TO THE FILE!!!
    print FILE "\n";
    print FILE browserTrackString("bam", ${sortBamFullFilename});
    print FILE browserTrackString("bigWig", ${bigWigOutFile});
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

convert_SAM_or_BAM_for_Genome_Browser.pl <ONE SINGLE INPUT BAM / SAM FILE>

Processes a SAM or BAM alignment file to generate files for the UCSC Genome Browser.

To handle MULTIPLE files in the Bash shell, try this: for f in `ls */accepted_hits.bam`; do convert_SAM_or_BAM_for_Genome_Browser.pl $f; done

Options:
    * --nowig: If you do not want to generate a wiggle track, you can say `--nowig` and omit the --fasta file.
    * --alreadysorted: If your file is already coordinate-sorted, then you can say `--alreadysorted` to avoid re-sorting it.

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

Inputs:
   1. A SAM or BAM file, with aligned reads to a reference genome.

Note: if you do not want the bigwig file, you can omit generating it with --nowig 

Output: 2 files:
    1. Generates a big pileup track (bigBed format).
    2. A bigWig (bw) format file for the UCSC genome browser.

These files can be viewed in the genome browser.

OPTIONS:

   --nosort: Assume the input file is ALREADY sorted. Speeds things up.

   --nowig: Add this to NOT generate a wiggle track.
            Generating a wiggle track is SLOW.
            Also, if you specify --nowig, then you no longer need to specify
            a "fai" fasta index file.

EXAMPLES:

convert_SAM_or_BAM_for_Genome_Browser.pl alignment.bam
convert_SAM_or_BAM_for_Genome_Browser.pl --nowig  alignment.bam

Files that are generated from the input YOURFILE.sam:
 1. Browser.sort.YOURFILE.bam (sorted version of the BAM/SAM file) (track type=bam)
 2. Browser.sort.YOURFILE.bam.bai (BAM index file)
 3. Browser.tmp.YOURFILE.bedgraph (bedgraph intermediate file used to make the bigwig. Can be deleted!)
 4. Browser.YOURFILE.bw (bigWig file for Genome Browser)
 5. Browser.Track.Descriptions.YOURFILE.txt (track descriptions that you paste into the Genome Browser Custom Tracks)

Then you will just need to make a genome browser track description for those two files.

You can find one in file 5, the browser track description. Or make your own, like this sample one:

track type=bam name="Control_EB" color=0,128,255 bigDataUrl="http://lighthouse.ucsf.edu/public_files_no_password/browser_custom_bed/sorted_ctrl.bam"

  (Note that there is no need to mention the .bam.bai file in the track description,
   but it is implicitly used by the Genome Browser and must be present!)

This script just runs the instructions found at:
      http://genome.ucsc.edu/goldenPath/help/bam.html

This script will probably take about 15 minutes to run on a single 5 GB input SAM file.

Generating the bigWig file is SLOW. You can prevent this by specifying --nowig.

This script generates a single 5 kilobyte chromosome length temp file in /tmp/, which it deletes
upon successful script completion.


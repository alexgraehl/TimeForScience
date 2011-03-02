#!/usr/bin/perl

## This script converts a BAM or SAM file into the required parts for the genome browser.

## by Alex Williams, Feb. 2011

use strict;  use warnings;  use diagnostics;

## Convert SAM or BAM to Genome Browser format

print STDERR "This is a script that will generate UCSC-genome-browser-ready BAM files.\n"
    . "It will convert / produce the UCSC-ready files from any SAM or BAM file you want.\n"
    . "The BAM/SAM files can be HUGE (10+ GB), but can be hosted locally. Then you tell the UCSC Genome Browser\n"
    . "where your BAM files are, and it magically uses the BAM files from your local web server, without copying the entire file.\n"
    . "\n"
    . "Files that are generated from the input YOURFILE.sam:\n"
    . " 1. YOURFILE.bam\n"
    . " 1. sorted_YOURFILE.bam (sorted version of the previous file)\n"
    . " 1. sorted_YOURFILE.bam.bai (BAM index file)\n"
    . "\n"
    . "This script just runs the instructions found at:\n"
    . "      http://genome.ucsc.edu/goldenPath/help/bam.html\n"
    . "\n"
    . "Note also that <samtools> must be installed. You can install that with apt-get install samtools, if it isn't already installed.\n"
    . "\nThis script will probably take about 10 minutes to run on its single input file.\n";


unless (scalar(@ARGV) == 1) { die "ARGUMENT ERROR: This script takes exactly one argument: one single file name of a BAM or SAM file!\n[Quitting now.]\n"; }


my $file = $ARGV[0];


unless (-e $file) { die "ARGUMENT ERROR: The argument to this program must be a single BAM/SAM file that already exists.\n"; }

unless (-r $file) { die "FILE ERROR: The argument to this program must be a single file that already exists AND is also readable. We couldn't read the file <$file> though!\n"; }

unless (not (-z $file)) { die "FILE ERROR: The file <$file> was 0-length. It's not a real BAM/SAM file!\n"; }



# Everything about BAM -> browser is described here http://genome.ucsc.edu/goldenPath/help/bam.html

my $fileIsSAM = ($file =~ m/\.sam$/i);
if ($fileIsSAM) { print "<$file> appears to be a SAM file.\n"; }

my $fileIsBAM = ($file =~ m/\.bam$/i);
if ($fileIsBAM) { print "<$file> appears to be a BAM (Binary SAM) file.\n"; }


if (!$fileIsBAM && !$fileIsSAM) { die "FILE ERROR: File must be either SAM or BAM!\n"; }

my $bamFilename;

if ($fileIsSAM) {
    ## User specified a SAM file
    $bamFilename = $file;
    $bamFilename =~ s/\.sam$/\.bam/i;
    if (not(-e $bamFilename)) {
	print "Time to make the bam file named <$bamFilename> from the input SAM file that was named <$file>... this will take a minute or two...\n";
	system("samtools view -S -b $file > $bamFilename");
    } else {
	print "Not remaking the BAM file named <$bamFilename>, because it already existed.\n";
    }
}


if ($fileIsBAM) {
    $bamFilename = $file; # that was easy, the user's specified file was ALREADY a bam file
}

my $bamPrefixWithoutFileExtension = $bamFilename;
$bamPrefixWithoutFileExtension =~ s/\.bam$//i;

my $sortedOutNameWithoutFileExtension = "sorted_${bamPrefixWithoutFileExtension}";
my $sortCmd = "samtools sort $bamFilename $sortedOutNameWithoutFileExtension"; ## <-- automatically adds the ".bam" extension for some weird reason
print("Ok, we're going to run the samtools sort command now...\n>> Running the command: $sortCmd\n...\n");
system($sortCmd);

my $indexCmd = "samtools index ${sortedOutNameWithoutFileExtension}.bam";
print("Ok, we're going to run the samtools index command now...\n>> Running the command: $indexCmd\n...\n");
system($indexCmd);

print "DONE!\n\n";
print "Here are the final output files that you will probably want to put on your server:\n"
    . " - ${sortedOutNameWithoutFileExtension}.bam\n"
    . " - ${sortedOutNameWithoutFileExtension}.bam.bai\n";
print "\n";
print "1. You will probably want to move these files to your web server now.\n"
    . "   For my server, the command: scp sorted* lighthouse.ucsf.edu:\n"
    . "   copies everything from this machine to my home directory on that machine.\n";
print "\n";
print "2. Make sure that the files are READABLE by everyone, and that the directories they\n"
    . "   are in are also readable and executable by everyone. Double-check that you can access\n"
    . "   your files in a regular web browser.\n";
print "\n";

my $thingToTypeForBrowser = qq{track type=bam name="Control_EB" bigDataUrl=http://lighthouse.ucsf.edu/public_files_no_password/browser_custom_bed/sorted_ctrl.bam"};
print "3. Then you can go to the genome browser gateway page (http://genome.ucsc.edu/cgi-bin/hgGateway)\n"
    . "   and click the \"add custom tracks\" button near the top. Then, in \"Paste URLs or data\", type\n"
    . "   " . $thingToTypeForBrowser . "\n"
    . "   (replacing sorted_ctrl.bam with your actual .bam file---note that the .bai files aren't mentioned here, but MUST be in the same directory.\n";
print "\n";
print "4. Note that if you want to save this custom track or send it to others, you will need\n"
    . "   a custom genome browser SESSION---just copying / saving the URL is not sufficient!\n";
#print "Assuming that your files are hosted on lighthouse, here are the lines to paste into the table browser to make things work:\n";
#print "track type=bam name="Control_EB" bigDataUrl=http://lighthouse.ucsf.edu/public_files_no_password/browser_custom_bed/sorted_ctrl.bam"
print "\n";





    

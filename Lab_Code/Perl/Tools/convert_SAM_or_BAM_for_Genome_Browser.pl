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

my $genomeFastaFile = undef;
my $faiFile = undef;  ## <-- "fai" means genome "fasta index (fai)" file
my $makeWig = 1; ## By default, generate a bigwig track too. Specify "nowig" to avoid this.

GetOptions("help|?|man"        => sub { printUsageAndQuit(); }
	   , "wig!" => \$makeWig ## specify "--nowig" to avoid making a wiggle track
	   , "fasta=s" => \$genomeFastaFile
    ) or printUsageAndQuit();

if (scalar(@ARGV) != 1) { die "ARGUMENT ERROR: This script takes the name of a genome fasta file with --fasta=/location/to/species.fa and then exactly one non-named argument: one single file name of a BAM or SAM file!\nExample: convert_SAM_or_BAM_for_Genome_Browser.pl --fasta=/path/to/fasta/hg19.fa mySamFile.sam\n[Quitting now.]\n"; }



if ($makeWig) {
    if (!defined($genomeFastaFile)) { die "ARGUMENT ERROR: You must supply an input genome fasta file (with --fasta=/location/to/species.fa) for computing the wiggle tracks, OR you can specify --nowig on the command line to disable wiggle track generation."; }
    ## If we want to make a wiggle track, we will need a FAI fasta index file.
    if (!(-r $genomeFastaFile) or !(-e $genomeFastaFile) or (-z $genomeFastaFile)) {
	die "We could not read the specified --fasta file (genome sequence file) at:\n   $genomeFastaFile\n   Please check to make sure this is a valid and readable fasta file!\n";
    }

    $faiFile = ($genomeFastaFile . ".fai");  ## <-- "fai" means genome "fasta index (fai)" file
    
    if (!(-e $faiFile) && (-s $faiFile)) { ## -s means "return size of file" (to check for zero-length files)
	print STDOUT "There wasn't already a \".fai\" file where the .fasta file was located, so we are now automatically genearting the fasta index file <$faiFile> from the specified genome fasta file <$genomeFastaFile>...\n";
	system(qq{samtools faidx $genomeFastaFile}); ## generate an index file
	print STDOUT "Done. Wrote <$faiFile> to the filesystem.";
    }
}

print STDOUT "File to be processed by convert_SAM_or_BAM_for_Genome_Browser:\n";
foreach (@ARGV) {
    print STDOUT "  - $_\n";
}

my $originalInputFilename = $ARGV[0];

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

my $bamPrefixWithoutFileExtension = $bamFilename; ## Don't just get the basename: get the WHOLE PATH that led up to this file. basename($bamFilename);
$bamPrefixWithoutFileExtension =~ s/\.bam$//i;
$bamPrefixWithoutFileExtension =~ s/[\/:;,]/_/g; ## slashes and ':;,' characters go to underscores

my $browserTrackDescriptionFile = "Browser_Track_Descriptions.${bamPrefixWithoutFileExtension}.txt";
my $sortBamFilePrefix   = "Browser_sorted.${bamPrefixWithoutFileExtension}";
my $sortBamFullFilename = "${sortBamFilePrefix}.bam";
my $bamIndexOutfile     = "${sortBamFilePrefix}.bam.bai";
my $wigIntermediateFile = "Browser_tmp.${bamPrefixWithoutFileExtension}.wig";
my $bigWigOutFile       = "Browser.${bamPrefixWithoutFileExtension}.bigwig.bw";

if ((-e $sortBamFullFilename) && (-s $sortBamFullFilename) and (-e $bamIndexOutfile) && (-s $bamIndexOutfile)) {
    print STDOUT "[Skipping] We are NOT continuing with the sorted-BAM-file generation, because such a file already exists. Remove it if you want to recompute it!\n"
} else {
    my $sortCmd = "samtools sort $bamFilename $sortBamFilePrefix"; ## <-- Note: this should NOT have the .bam extension, because samtools automatically adds the ".bam" extension
    print(">> Running the SAMTOOLS sort command: $sortCmd\n...\n");
    system($sortCmd);
    
    my $indexCmd = "samtools index ${sortBamFullFilename}";
    print(">> Running the SAMTOOLS index command: $indexCmd\n...\n");
    system($indexCmd);
}

if ($makeWig) {
    if (-e $bigWigOutFile) {
	print STDOUT "[Skipping] We are NOT continuing with the generation of a bigwig file, because the bigwig file $bigWigOutFile already exists. Remove it if you want to recompute it!\n";
    } else {
	print "Now generating a BIG WIG browser wiggle track named $wigIntermediateFile (this is slow, and can take up to an hour per input RNASeq file!)...\n";
	
	if (-e $wigIntermediateFile && (-s $wigIntermediateFile > 0)) {
	    print STDOUT "[Skipping] We are NOT creating the wiggle temp file <$wigIntermediateFile>, because it already exists and has non-zero size. Remove it if you want to recompute it.\n";
	} else {
	    my $wigCmd1 = (qq{samtools mpileup -f $faiFile $bamFilename }
			   . qq(  | awk '{print \$1, \$2-1, \$2, \$4}' )
			   . qq{    > $wigIntermediateFile}); ## <-- generates a wig file
	    datePrint(": Now running this command:\n  $wigCmd1\n"); 
	    system($wigCmd1);
	}

	## -clip means "allow weird errant entries off the end of the chromosome, rather than exploding". This is important, because otherwise wigToBigWig will quit with errors like "something went off the end of chr12_random"
	my $wigCmd2 = (qq{wigToBigWig -clip $wigIntermediateFile $faiFile $bigWigOutFile});
	datePrint("Now running this command:\n  $wigCmd2\n");
	system($wigCmd2); ## makes a BIGWIG file from the WIG file
    }
} else {
    print STDOUT qq{[Skipping] the generation of a bigWig file, because "--nowig" was specified on the command line.\n};
}





sub browserTrackString($$) {
    my ($type, $filename) = @_; ## input arguments: type must be bigWig or bam
    if ($type ne "bigWig" and $type ne "bam") { die "Sorry, we only know how to write BIGWIG and BAM files. Problem!\n"; }
    my $redColor = "255,0,0";
    my $blueColor = "0,0,255";
    my $color;
    if ($type eq "bigWig") { $color = qq{color="$redColor"}; }
    if ($type eq "bam") { $color = qq{colorByStrand="$redColor $blueColor"}; }
    my $url = "http://lighthouse.ucsf.edu/public_files_no_password/browser_custom_bed/YOUR_FILE_LOCATION/${filename}";
    my $str = qq{track type=${type} name="${filename} (${type})" description="${filename} (${type})" bigDataUrl="${url}" visibility=2 ${color}\n};
    return($str);
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
#print "Assuming that your files are hosted on lighthouse, here are the lines to paste into the table browser to make things work:\n";
#print "track type=bam name="Control_EB" color=0,128,255 bigDataUrl=http://lighthouse.ucsf.edu/public_files_no_password/browser_custom_bed/sorted_ctrl.bam"
print "\n";


__DATA__

convert_SAM_or_BAM_for_Genome_Browser.pl --fasta=<genome_fasta_file> <INPUT_SAM_FILE>

Process a SAM or BAM alignment file to generate files for the UCSC Genome Browser.

This is a script that will generate UCSC-genome-browser-ready BAM files.
It will convert / produce the UCSC-ready files from any SAM or BAM file you want.
The BAM/SAM files can be HUGE (10+ GB), but can be hosted locally. Then you tell the UCSC Genome Browser
where your BAM files are, and it magically uses the BAM files from your local web server, without copying the entire file.

Inputs:
   1. --fasta=<genome.fasta.file.fa> . The genome of the species in question. Used to generate wig files.
   2. A SAM or BAM file, with aligned reads to a reference genome.

Note: if you do not have a genome fasta file, you can omit generating the wig files by saying --nowig instead
of specifying a fasta file.

Output: 2 files:
    1. Generates a big pileup track (bigBed format).
    2. A bigWig format file for the UCSC genome browser. The "fai" fasta index file is required in order to generate this bigWig file.

These files can be viewed in the genome browser.


OPTIONS:

   --fai=filename: Required UNLESS you specify "--nowig"
                   A fasta index file. Format is: some_species.fai or some_species.fa.fai

   --nowig: Add this to NOT generate a wiggle track.
            Generating a wiggle track is SLOW.
            Also, if you specify --nowig, then you no longer need to specify
            a "fai" fasta index file.

EXAMPLES:

 convert_SAM_or_BAM_for_Genome_Browser.pl --fasta=/Genomes/Mouse/mm9.fa  alignment.bam
 convert_SAM_or_BAM_for_Genome_Browser.pl --fasta=/Genomes/HS/hg19.fa  alignment.sam
 convert_SAM_or_BAM_for_Genome_Browser.pl --fasta=/work/Common/Data/Genome/hg19wholeGenomeFasta/hg19.fa  alignment.bam
 convert_SAM_or_BAM_for_Genome_Browser.pl --nowig  alignment.bam

Files that are generated from the input YOURFILE.sam:
 1. Browser_sorted_YOURFILE.bam (sorted version of the BAM/SAM file) (track type=bam)
 2. Browser_sorted_YOURFILE.bam.bai (BAM index file)
 3. Browser_tmp.YOURFILE.wiggle_track.bed (regular wiggle track -- BED file. Can be deleted.)
 4. Browser.YOURFILE.bigwig.bw (bigWig file for Genome Browser)
 5. Browser_Track_Descriptions._YOURFILE.txt (track descriptions that you paste into the Genome Browser Custom Tracks)

Then you will just need to make a genome browser track description for those two files.

You can find one in file 5, the browser track description. Or make your own, like this sample one:

track type=bam name="Control_EB" color=0,128,255 bigDataUrl="http://lighthouse.ucsf.edu/public_files_no_password/browser_custom_bed/sorted_ctrl.bam"

  (Note that there is no need to mention the .bam.bai file in the track description,
   but it is implicitly used by the Genome Browser and must be present!)

This script just runs the instructions found at:
      http://genome.ucsc.edu/goldenPath/help/bam.html

Note also that <samtools> must be installed. You can install that with apt-get install samtools, if it is not already installed.

This script will probably take about 30 minutes (!) to run on a single 5 GB input SAM file.

Generating the bigWig file is SLOW. You can prevent this by specifying --nowig.

#!/usr/bin/perl

use constant USE_COLORS_CONSTANT => 1; ## 1 = true, 0 = false
use POSIX      qw(ceil floor);
use List::Util qw(max min);
if (USE_COLORS_CONSTANT) { use Term::ANSIColor; }

use File::Basename;
use Getopt::Long;
#no warnings 'numeric';
#use Scalar::Util;
#print Scalar::Util::looks_like_number($string), "\n";

use strict;  use warnings;  use diagnostics;

sub main();

my $INTERSECT_BED_EXE = "/work/Apps/Bio/bedtools/bin/intersectBed";
my $globalAnnotDir = "./Z_AGNOMEN_DATA"; # default!

my $AGNOMEN_ENSEMBL_GRABBER_SCRIPT = "/work/Common/Code/ProjectCode/0_Analyze_0277_Agnomen_Annotate/agnomen-ensembl-hg19-grabber.pl";

my %GLOBAL_ANNOT_PATHS = (); ## New hash that will keep track of all the annotation files. Keys are the descriptions of the files, values are the actual paths of the annotation files.

sub systemBash {
    my @args = ( "bash", "-c", shift );
    system(@args);
}

sub agwSystemDieOnNonzero($) {
    my ($cmd) = @_;
    my $exitCode = system($cmd);
    if ($exitCode != 0) { die "QUITTING: Exit code of <$cmd> was non-zero (specifically, it was <$exitCode>).\n";
}

sub stderrPrint {
    print STDERR @_;
}

sub colorString($) {
    # Requires: use Term::AnsiColor at the top of your file. Note: you have to
    # PRINT the result of this function! Only sets the output color if the output
    # is a TERMINAL. If the output is NOT a terminal, then colorizing output
    # results in lots of garbage characters (the color control characters)
    # written to the screen.
    my ($theColor) = @_;
    if (-t STDOUT && USE_COLORS_CONSTANT) { # <-- checks to see if STDOUT goes directly to the terminal (instead of, say, outputting to a file with a redirect)
	return color($theColor);
    } else {
	return "";
    }
}

sub agwFileHasContents($) { my ($fname) = @_; return((-e $fname) && (-s $fname > 0)); }
sub quitWithUsageError($) { print($_[0] . "\n"); printUsage(); print($_[0] . "\n"); exit(1); }
sub printUsageAndQuit() { printUsage(); exit(1); }
sub printUsage() {  print STDOUT <DATA>; }



sub agnomenRemoteDownloadAnnot($$$$$$) {
    # Tries to generate the file $finalName, UNLESS that file already exists and is valid.
    ## Arguments:
    ## 0. The "short name" that gets used as the hash ID and also in the filename of the output
    ## 1. The URL to get
    ## 2. the genome build (e.g. hg19, hg18, mm9)
    ## 3. the name to save the RAW file into
    ## 4. postprocessingCmd: something to do with the file after downloading it. Set it to "undef" if it isn't defined
    ## 5. finalName: the final filename after processing
    my ($shortName, $url, $theBuild, $rawNameIncludingCompression, $postprocessingCmd, $finalName) = @_;

    my $localDir = "${globalAnnotDir}/${theBuild}"; #${species}/${theBuild}";
    my $rawFullPath = $localDir . "/" . $rawNameIncludingCompression;

    my $unzippedPath = $rawFullPath;
    $unzippedPath =~ s/[.](bz2|gz|zip)$//i;

    if ($rawFullPath =~ m/[ ,\t\s]/) { die "Uh oh, the local full path has a comma or whitespace in it, which is not legal for the directory name!!"; }
    my $curlCmd = 'curl --remote-name ' . '"' . $url . '"' ;

    system("mkdir -p $localDir");
    if ($rawNameIncludingCompression ne File::Basename::basename($url)) {
	stderrPrint("This is unusual; the local name of the file ($rawNameIncludingCompression) does not match the filename from the URL (the end of this URL: $url)!");
    }

    if (agwFileHasContents($unzippedPath)) {
	stderrPrint(colorString("cyan"));
	stderrPrint("[Skipping] re-download of <$rawFullPath> from <$url> -- the final downloaded file <$unzippedPath> already exists.\n");
	stderrPrint(colorString("reset"));
    } else {
	## Do we need to download the file?
	if (!agwFileHasContents($rawFullPath)) {
	    ## Better download it!
	    stderrPrint("Downloading $url --> $rawFullPath...\n");
	    agwSystemDieOnNonzero($curlCmd);
	} else {
	    stderrPrint("[Skipping] re-download of <$rawFullPath> from $url, as that file already exists.\n");
	}
	if ($rawFullPath =~ m/[.]gz$/) { agwSystemDieOnNonzero("gunzip --stdout $rawFullPath > $unzippedPath"); }
	if ($rawFullPath =~ m/[.]bz2$/) { agwSystemDieOnNonzero("bunzip2 --stdout $rawFullPath > $unzippedPath"); }
	if ($rawFullPath =~ m/[.]zip$/) { agwSystemDieOnNonzero("unzip $rawFullPath"); }
	stderrPrint("[Done] Successfully downloaded a file to the uncompressed location <$unzippedPath>\n");
    }

    ## Now the file should be downloaded, if necessary. We may also have to do a post-processing step, below.
    if (defined($postprocessingCmd) && $postprocessingCmd && length($postprocessingCmd) > 0) {
	if (agwFileHasContents($finalName)) {
	    stderrPrint(colorString("cyan"));
	    stderrPrint("[Skipping] re-generation of $finalName, as the output file already exists.\n");
	    stderrPrint(colorString("reset"));
	} else {
	    stderrPrint("Now running this post-processing command: $postprocessingCmd\n");
	    agwSystemDieOnNonzero($postprocessingCmd);
	    if (not -e $finalName) { die "Uh oh, we failed to generate the file <$finalName>! Maybe the shell command had a mis-typed final filename, or perhaps something else went wrong?\n"; }
	}
	
    } else {
	# No postprocessing to do, so we just want to use the regular file as an annotation file
	unless ($unzippedPath ne "$localDir/$finalName") die "Uh oh, the final file path was expected to be <$unzippedPath>, but it was manually specified as <$localDir/$finalName>.";
    }

    $GLOBAL_ANNOT_PATHS{$shortName} = "$localDir/$finalName"; ## <-- Add this file to the list of annotation files.
}

sub agnomenGenerateLocalAnnot($$$$$) {
    # Doesn't try to download anything; just runs '$cmd' UNLESS the $finalFilename already is a non-zero-length file.
    my ($shortName, $finalFilename, $theBuild, $cmd) = @_;
    my $localDir = "${globalAnnotDir}/${theBuild}";	
    my $finalFullPath = $localDir . "/" . $finalFilename;
    if (agwFileHasContents($finalFullPath)) {
	stderrPrint(colorString("cyan"));
	stderrPrint("[Skipping] re-generation of file <$finalFilename>, which already exists in <$finalFullPath>!\n");
	stderrPrint(colorString("reset"));
    } else {
	my $fullCommand = "cd $localDir && $cmd";
	stderrPrint("[Generating the file $finalFilename], by running the command:\n$fullCommand\n\n\nRunning now...\n");
	agwSystemDieOnNonzero($fullCommand); ## cd into the local
    }
    if (exists($GLOBAL_ANNOT_PATHS{$shortName})) { die "Programming error 19A: Uh oh, there is ALREADY an annotation file with the short name <$shortName>. Can't have two duplicate short names! Check the source code and fix it.\n"; }
    $GLOBAL_ANNOT_PATHS{$shortName} = $finalFullPath;
}






# ==1==
sub main() { # Main program
    my ($delim) = "\t";
    my ($genomeBuild) = undef;
    my ($shouldUpdate) = 0;
    my ($decimalPlaces) = 4; # How many decimal places to print, by default
    $Getopt::Long::passthrough = 1; # ignore arguments we don't recognize in GetOptions, and put them in @ARGV

    GetOptions("help|?|man" => sub { printUsageAndQuit(); }
	       , "delim|d=s" => \$delim
	       , "species|s=s" => \$genomeBuild
	       , "update!" => \$shouldUpdate
	       , "annotdir|a" => \$globalAnnotDir
	       , "dp=i" => \$decimalPlaces
	) or printUsageAndQuit();

    stderrPrint(colorString("green"));
    stderrPrint("-------------------------------\n");
    stderrPrint("Starting...!\n");
    stderrPrint("-------------------------------\n");
    stderrPrint(colorString("reset"));

    if (1 == 0) {
	quitWithUsageError("1 == 0? Something is wrong!");
    }

    if (not -e $INTERSECT_BED_EXE) { die "Error 95a: <intersectBed> was not found!\nWe REQUIRE the program <intersectBed> to be installed. It is part of BedTools. We expected it to be located in <${INTERSECT_BED_EXE}>, but there was no file there at all! Check to make sure this path exists.\n"; }
    if (not -x $INTERSECT_BED_EXE) { die "Error 99x: <intersectBed> was not executable by this user!\nWe REQUIRE <intersectBed>, located in <${INTERSECT_BED_EXE}>, to be executable. However, the file at that location was not executable by this user! You may have to use chmod to fix this.\n"; }

    my $numUnprocessedArgs = scalar(@ARGV);
    #if ($numUnprocessedArgs != 2) {
#	quitWithUsageError("Error in arguments! You must send TWO filenames to this program.\n");
#    }

    my @inputFiles = (); # Input files to annotate

    foreach (@ARGV) { # these were arguments that were not understood by GetOptions
	if (-e $_) {
	    push(@inputFiles, $_);
	} else {
	    die "Uh oh, the input argument <$_>, which we expected to be a file, apparently was not an actual file!";
	}
    }

    unless (-d $globalAnnotDir) { quitWithUsageError(">>> ERROR >>> The annotation directory ($globalAnnotDir) does not exist!") }
    unless (defined($genomeBuild) && (length($genomeBuild) > 0)) { quitWithUsageError(">>> ERROR >>> it is MANDATORY to specify a genome build (example: hg19 or mm9). You do this with a command like this:   --species=hg19  or  -s hg19 \nCheck the annotation directory for a list of valid species.\n"); }

    stderrPrint("* Searching for annotation files in: <$globalAnnotDir>\n");
    stderrPrint("* Using the specified genome build: <$genomeBuild>\n");
    stderrPrint("* Searching for BED annotation files in the following file path: <$globalAnnotDir/$genomeBuild/>\n");

    if (scalar(@inputFiles) == 0 && !$shouldUpdate) {
	## Need to have some files to annotate!
	quitWithUsageError(">>> ERROR >>> We expected AT LEAST ONE valid file to annotate to be specified on the command line. We didn't get any, however.\n");
    }

    
    if ($shouldUpdate) {
	my $DOLLAR_SIGN = '$';
	agnomenRemoteDownloadAnnot("Mouse_Ensembl_GTF", "ftp://ftp.ensembl.org/pub/release-68/gtf/mus_musculus/Mus_musculus.GRCm38.68.gtf.gz", "mm9", "Mus_musculus.GRCm38.68.gtf.gz", undef, "Mus_musculus.GRCm38.68.gtf");
	agnomenGenerateLocalAnnot("Human_Ensembl_BED", "human_ensembl_via_perl_api.tab", "hg19",
				  ( qq{ if [ ! -f human_ensembl_1_raw.tab.gz ]; then echo "\n>>>[NOTE] -- The ENSEMBL GRABBER takes about 18 hours to download data via Perl API! If you are reading this message, be prepared to wait QUITE A WHILE for the data to get downloaded!\n"; $AGNOMEN_ENSEMBL_GRABBER_SCRIPT | gzip > human_ensembl_1_raw.tab.gz ; fi; }
				    . qq{ zcat human_ensembl_1_raw.tab.gz | grep -v '^[#]' | grep -v '^$DOLLAR_SIGN' | sort -k 2,2 -k 3,3g | uniq > human_ensembl_2.tmp ; }
				    . qq{ sed 's/^/chr/' human_ensembl_2.tmp > human_ensembl_3.tmp ; } ## add 'chr' to the beginning of each line
				    . qq{ cat human_ensembl_2.tmp >> human_ensembl_3.tmp ; }
				    . qq{ sort human_ensembl_3.tmp > human_ensembl_via_perl_api.tab ; }
				    . qq{ /bin/rm human_ensembl_[123].tmp ; }
				  ));
    }
    
    #print "\t" . sprintf("%.${numDecimalPointsToPrint}f", $levFraction);
    #my @col1 = @{readFileColumn($filename1, 0, $delim)};
    
    #my $annot1 = "AnnotatedFeatures.gff"; #"hb.bed"; #"human_ens_manual.bed"; #"b.gtf"; #"AnnotatedFeatures.gff"; #"annot.bed";
    #stderrPrint("Warning: using a manually-selected annotation file here.\n");

    if (scalar(keys(%GLOBAL_ANNOT_PATHS)) == 0) {
	die "Uh oh, there were no annotation files, for some reason! Maybe you need to run 'agnomen-gene-annotate.pl --update' to re-download / re-generate them?\n";
    }
    
    ## For each file, we're going to annotate that file and write an output file. For each file, we're writing one output file per annotation type!
    ## Therefore, if you have (say) 3 input files to annotate, and 4 types of annotation, you will end up with (3*4 = 12) final output files.
    for my $input (@inputFiles) {
	if ($input !~ m/[.]bed$/i) {
	    ## Expecting the $input to be a BED file--- check to see if it has ".bed" as an extension?
	    stderrPrint("Warning: the input file $input was EXPECTED to be a .bed file, but it appears not to have the .bed filename extension. Double-check to make sure this is the right file! Continuing on anyway, assuming that it is a bed file...\n");
	}
	## 
	while (my($annotDescription, $annotFile) = each(%GLOBAL_ANNOT_PATHS)) {
	    ## annotFile better exist, or we're in trouble!!
#	    if ($annotFile =~ m/.gz$/) {
#		$annotFile = "<(zcat $annotFile)"; # bash subshell
#		print STDERR "SUBSHELL\n";
#	    }

	    my $outputFile = "ANNOT_" . $input . "--" . $annotDescription;
	    $outputFile =~ s/[;:,\/\\]/_/g; ## Remove potentially "unsafe" characters from the output filename
	    my $cmd = qq{ $INTERSECT_BED_EXE -wao -a ${input} -b ${annotFile} > $outputFile};

	    stderrPrint(colorString("yellow"));
	    stderrPrint(qq{[ANNOTATING]...\n    Running this command: $cmd\n});
	    stderrPrint(colorString("reset"));
	    my $intersectStatus = systemBash($cmd);
	    
	    if ($intersectStatus != 0) {
		stderrPrint("Uh oh, non-zero exit status!");
	    }
	}
    }

    stderrPrint(colorString("green"));
    stderrPrint("-------------------------------\n");
    stderrPrint("[Done!]\n");
    stderrPrint("-------------------------------\n");
    stderrPrint(colorString("reset"));
} # end main()


main();


END {
    # Runs after everything else.
    # Makes sure that the terminal text is back to its normal color.
    stderrPrint(colorString("reset"));
}

exit(0);
# ====

__DATA__

agnomen-gene-annotate.pl  [OPTIONS]  -s GENOME_BUILD   INPUT_FILE_TO_ANNOTATE.bed

by Alex Williams, 2012

GENOME_BUILD must be a UCSC-style genome build (examples: hg19, mm9)

The input file to annotation must be a valid tab-delimited BED-format file.

Example BED file:
    chr1   <tab>   100    <tab>   200
    chr2   <tab>  1150    <tab>  1300

If you want to learn more about file formats (GFF, etc), check: http://genome.ucsc.edu/FAQ/FAQformat.html

The local annotation is, by default, expected to be found in a directory called Z_AGNOMEN_DATA.

Agnomen is supposed to be able to manage the Z_AGNOMEN_DATA directory by itself, so you should in theory
not have to deal with this directory.

EXAMPLES:

agnomen-gene-annotate.pl  -s hg19  somefile.bed
  Annotates "somefile.bed" and outputs a new annotated file for each annotation type:
  Example output files that might be generated:
              ANNOT_somefile.bed--HUMAN_TRANSCRIPTS.txt
              ANNOT_somefile.bed--HUMAN_Structure_Pred.txt

agnomen-gene-annotate.pl  --update
  If necessary, re-downloads the annotation files.

CAVEATS:

NA

OPTIONS:

  --delim = DELIMITER  (Default: tab)
    or: -d DELIMITER
       Sets the input delimiter to DELIMITER.

  --species = STRING  (Default: undefined)
    or: -s STRING
       Sets which genome build should be used for annotation. Example: hg19 (a human genome), mm9 (a mouse genome).
       Currently supported: hg19, mm9.
       Default is to use output for ALL the different species.

  --annotdir = STRING (Default: Z_AGNOMEN_DATA)
    or: -a STRING
       A path on the filesystem where we keep the annotation data. The annotation data must be in species-specific
       sub-folders. Example filesystem:  MY_DATA/hg19/human_annot.bed
                                         MY_DATA/mm10/mouse_annot_1.bed
                                         MY_DATA/mm10/mouse_annot_extended.bed
       In the example above, you would specify --annotdir=MY_DATA

  --update
       Tells Agnomen to check for (and potentially update) newer annotation files. Does not fully work yet!
       Default is --noupdate.


KNOWN BUGS:

  None known.

TO DO:

  Add ???.

--------------

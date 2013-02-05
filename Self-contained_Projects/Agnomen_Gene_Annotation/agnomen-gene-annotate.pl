#!/usr/bin/perl

use constant USE_COLORS_CONSTANT => 1; ## 1 = true, 0 = false
use POSIX      qw(ceil floor);
use List::Util qw(max min);
use Cwd;
if (USE_COLORS_CONSTANT) { use Term::ANSIColor; }

use File::Basename;
use Getopt::Long;
#no warnings 'numeric';

use strict;  use warnings;  use diagnostics;

sub main();

my $INTERSECT_BED_EXE = "/work/Apps/Bio/bedtools/bin/intersectBed";
my $globalAnnotDir = "./Z_AGNOMEN_DATA"; # default!

my $DATABASE_INPUT_STRING_DELIM = ','; # comma-delimiting in the input string, example: --databases=human_ensembl,hg19browser,somethingelse
my $AGNOMEN_ENSEMBL_GRABBER_SCRIPT = "/home/alexgw/TimeForScience/Self-contained_Projects/Agnomen_Gene_Annotation/ag-ensembl-grabber-for-agnomen.pl";

my %GLOBAL_ANNOT_PATHS = (); ## New hash that will keep track of all the annotation files. Keys are the descriptions of the files, values are the actual paths of the annotation files.

sub systemBash {
    my @args = ( "bash", "-c", shift );
    system(@args);
}

sub agwSystemDieOnNonzero($) {
    my ($cmd) = @_;
    my $exitCode = system($cmd);
    if ($exitCode != 0) {
	die "QUITTING: Exit code of <$cmd> was non-zero (specifically, it was <$exitCode>).\n";
    }
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


sub addPathToGlobalAnnotationHash($$) {
    my ($shortName, $filePath) = @_;
    if (exists($GLOBAL_ANNOT_PATHS{$shortName})) { die "Programming error 19A: Uh oh, there is ALREADY an annotation file with the short name <$shortName>. Can't have two duplicate short names! Check the source code and fix it.\n"; }
    if (!agwFileHasContents($filePath)) {
	die "ERROR: The file at path <$filePath> was expected to exist, but it does not appear to exist at that path!\n";
    }
    $GLOBAL_ANNOT_PATHS{$shortName} = $filePath; ## <-- Add this file's FULL PATH to the list of annotation files.
}

sub agnomenGetAnnot($$$$$$$$) {
    # Tries to generate the file $finalName, UNLESS that file already exists and is valid.
    ## Arguments:
    ## 0. The "short name" that gets used as the hash ID and also in the filename of the output
    ## 1. The URL to get
    ## 2. the genome build (e.g. hg19, hg18, mm9) for this specific file. Set to "undef" to include this annotation file for ALL genomes.
    ## 3. the name to save the RAW file into
    ## 4. postprocessingCmd: something to do with the file after downloading it. Set it to "undef" if it isn't defined
    ## 5. finalName: the final filename after processing
    ## 6: allowedToUpdateFiles: whether we are allowed to actually modify the files on disk (to update them) or not
    ## 7: genomeBuildToUse: what genome build the user has specified for their input files. We will only add files that ALSO match this genome build.
    my ($shortName, $url, $theBuildForThisFile, $rawZipped, $postprocessingCmd, $finalName, $allowedToUpdateFiles, $genomeBuildToUse) = @_;

    my $origWorkingDir = cwd();

    if (defined($genomeBuildToUse) && defined($theBuildForThisFile) && ($theBuildForThisFile ne $genomeBuildToUse)) {
	stderrPrint(colorString("cyan"));
	stderrPrint("Not updating or adding the <$theBuildForThisFile> annotation file <$shortName> to the global annotation list, as it is for a different species.\n");
	stderrPrint(colorString("reset"));
	return;
    }

    stderrPrint(colorString("green"));
    if (defined($genomeBuildToUse)) { stderrPrint("Adding the annotation file <$shortName> to the global annotation list, as it is a valid annotation file for <$genomeBuildToUse>.\n"); }
    stderrPrint(colorString("reset"));

    my $localDir    = "${globalAnnotDir}/${theBuildForThisFile}";
    my $unzippedRaw = $rawZipped;
    $unzippedRaw    =~ s/[.](bz2|gz|zip)$//i;

    if (!$allowedToUpdateFiles) {
	stderrPrint(colorString("cyan"));
	stderrPrint("[Skipping] re-updating of the file of <$finalName>: the --update flag was not passed into this script.\n");
	stderrPrint(colorString("reset"));
    } else {
	system("mkdir -p $localDir");
	chdir($localDir);

	if (defined($url) && $url) {
	    my $curlCmd = 'curl --remote-name ' . '"' . $url . '"' ;
	    if ($rawZipped ne File::Basename::basename($url)) {
		stderrPrint("This is unusual; the local name of the file ($rawZipped) does not match the filename from the URL (the end of this URL: $url)!");
	    }
	    
	    if (agwFileHasContents($unzippedRaw)) {
		stderrPrint(colorString("cyan"));
		stderrPrint("[Skipping] re-download of <$localDir/$rawZipped> from <$url> -- the final downloaded file <$localDir/$unzippedRaw> already exists.\n");
		stderrPrint(colorString("reset"));
	    } else {
		## Do we need to download the file?
		if (!agwFileHasContents($rawZipped)) {
		    ## Better download it!
		    stderrPrint(colorString("green"));
		    stderrPrint("Since we didn't find a file at <$rawZipped>, we will start the download here...\n");
		    stderrPrint("Downloading $url --> $localDir/$rawZipped...\n");
		    stderrPrint(colorString("reset"));
		    agwSystemDieOnNonzero($curlCmd);
		} else {
		    stderrPrint("[Skipping] re-download of <$rawZipped> from $url, as that file already exists.\n");
		}
		if ($rawZipped =~ m/[.]gz$/) { agwSystemDieOnNonzero("gunzip --stdout $rawZipped > $unzippedRaw"); }
		if ($rawZipped =~ m/[.]bz2$/) { agwSystemDieOnNonzero("bunzip2 --stdout $rawZipped > $unzippedRaw"); }
		if ($rawZipped =~ m/[.]zip$/) { agwSystemDieOnNonzero("unzip $rawZipped"); }
		stderrPrint("[Done] Successfully downloaded a file to the uncompressed location <$localDir/$unzippedRaw>\n");
	    }
	}
    
	## Now the file should be downloaded, if necessary. We may also have to do a post-processing step, below.
	if (defined($postprocessingCmd) && $postprocessingCmd && length($postprocessingCmd) > 0) {
	    if (agwFileHasContents($finalName)) {
		stderrPrint(colorString("cyan"));
		stderrPrint("[Skipping] re-generation of $localDir/$finalName, as the output file already exists.\n");
		stderrPrint(colorString("reset"));
	    } else {
		stderrPrint(colorString("green"));
		stderrPrint("Now running this post-processing command:\n");
		stderrPrint(colorString("reset"));
		stderrPrint("$postprocessingCmd\n");
		agwSystemDieOnNonzero($postprocessingCmd);
		if (not -e $finalName) { die "Uh oh, we failed to generate the file <$finalName>! Maybe the shell command had a mis-typed final filename, or perhaps something else went wrong?\n"; }
	    }
	} else {
	    # No postprocessing to do, so we just want to use the regular file as an annotation file
	    if ($unzippedRaw ne $finalName) { die "Uh oh, the final file path was expected to be <$unzippedRaw>, but it was manually specified as <$finalName>, which is NOT identical."; }
	}
    }

    chdir($origWorkingDir);

    addPathToGlobalAnnotationHash($shortName, "$localDir/$finalName");
}


sub handleUserSpecifiedDatabases($$) {
    my ($dataString, $verbose) = @_;
    # Input arg: delimited string (example: "this,is,a,string")
    # Return value: a hash, with all the elements from the string split up.
    # and: $referenceToHash: a reference (\%hash would be passed in) to a hash where we store the results.
    my %newHash = ();
    if (defined($dataString)) {
	($verbose) && stderrPrint(colorString("green"));
	($verbose) && stderrPrint(qq{[INFO] The user specified that only the following annotation databases should be used: $dataString\n});
	($verbose) && stderrPrint(colorString("reset"));
	my @splitUpDb = split(/$DATABASE_INPUT_STRING_DELIM/, $dataString);
	for my $theItem (@splitUpDb) {
	    $newHash{$theItem} = 1; ## Let's remember that we SHOULD be using this database that the user specified on the command line!
	    $newHash{lc($theItem)} = 1; ## Let's also remember the LOWER CASE version of this name.
	}
	(scalar(keys(%newHash)) > 0) or die "The user specified some specific databases to use with the --databases (or -d) flag, but we couldn't interpret those as any actual valid databases! We ended up with ZERO databases. Note that this needs to be a comma-delimited string with NO SPACES in it.\n";
    } else {
	($verbose) && stderrPrint(colorString("green"));
	($verbose) && stderrPrint(qq{[INFO] Since no particular databases were specified, we will use ALL valid annotation databases.\n});
	($verbose) && stderrPrint(colorString("reset"));
    }
    return(%newHash);
}

# ==1==
sub main() { # Main program
    my $delim = "\t";
    my $genomeBuild = undef;
    my $shouldUpdate = 0;
    my $databaseStr = undef; # comma-delimited string
    my $longOutputNames = 0;
    $Getopt::Long::passthrough = 1; # ignore arguments we don't recognize in GetOptions, and put them in @ARGV

    GetOptions("help|?|man" => sub { printUsageAndQuit(); }
	       , "delim|d=s" => \$delim
	       , "species|s=s" => \$genomeBuild
	       , "update!" => \$shouldUpdate
	       , "longnames|longname!" => \$longOutputNames
	       , "annotdir|a=s" => \$globalAnnotDir
	       , "databases|database|db=s" => \$databaseStr # comma-delimited string
	) or printUsageAndQuit();

    stderrPrint(colorString("green"));
    stderrPrint("-------------------------------\n");
    stderrPrint("Starting...!\n");
    stderrPrint("-------------------------------\n");
    stderrPrint(colorString("reset"));


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

    unless (-d $globalAnnotDir) { quitWithUsageError(">>> ERROR >>> The annotation directory ($globalAnnotDir) does not exist! This directory must exist, as this is where we are putting all the annotation we download.\n") }
    if (not -e $INTERSECT_BED_EXE) { die "Error 95a: <intersectBed> was not found!\nWe REQUIRE the program <intersectBed> to be installed. It is part of BedTools. We expected it to be located in <${INTERSECT_BED_EXE}>, but there was no file there at all! Check to make sure this path exists.\n"; }
    if (not -x $INTERSECT_BED_EXE) { die "Error 99x: <intersectBed> was not executable by this user!\nWe REQUIRE <intersectBed>, located in <${INTERSECT_BED_EXE}>, to be executable. However, the file at that location was not executable by this user! You may have to use chmod to fix this.\n"; }
    
    my %databaseHash = handleUserSpecifiedDatabases($databaseStr, "verbose_messages");
    
    stderrPrint(colorString("green"));
    stderrPrint("* Searching for annotation files in: <$globalAnnotDir>\n");
    if (defined($genomeBuild)) { stderrPrint("* Using the specified genome build: <$genomeBuild>\n"); }
    if (defined($genomeBuild)) { stderrPrint("* Searching for BED annotation files in the following file path: <$globalAnnotDir/$genomeBuild/>\n"); }
    stderrPrint(colorString("reset"));

    if (scalar(@inputFiles) == 0 && !$shouldUpdate) {
	## Need to have some files to annotate!
	quitWithUsageError(">>> ERROR >>> We expected AT LEAST ONE valid file to annotate to be specified on the command line. We didn't get any, however.\n");
    }

    if ($shouldUpdate && scalar(@inputFiles) != 0) {
	quitWithUsageError(">>> ERROR >>> If you specify --update, you should NOT also specify any files to annotate!\n");
    }

    my $DOLLAR_SIGN = '$';
    agnomenGetAnnot("Mouse_Ensembl", "ftp://ftp.ensembl.org/pub/release-68/gtf/mus_musculus/Mus_musculus.GRCm38.68.gtf.gz"
		    , "mm9", "Mus_musculus.GRCm38.68.gtf.gz"
		    , undef
		    , "Mus_musculus.GRCm38.68.gtf", $shouldUpdate, $genomeBuild);
    agnomenGetAnnot("Human_Ensembl", undef, 
		    , "hg19", "human_ensembl_via_perl_api.tab"
		    , ( qq{ if [ ! -f human_ensembl_1_raw.tab.gz ]; then echo "\n>>>[NOTE] -- The ENSEMBL GRABBER takes about 18 hours to download data via Perl API! If you are reading this message, be prepared to wait QUITE A WHILE for the data to get downloaded!\n"; $AGNOMEN_ENSEMBL_GRABBER_SCRIPT | gzip > human_ensembl_1_raw.tab.gz ; fi; }
			. qq{ zcat human_ensembl_1_raw.tab.gz | grep -v '^[#]' | grep -v '^$DOLLAR_SIGN' | sort -k 2,2 -k 3,3g | uniq > human_ensembl_2.tmp ; }
			. qq{ sed 's/^/chr/' human_ensembl_2.tmp > human_ensembl_3.tmp ; } ## add 'chr' to the beginning of each line
			. qq{ cat human_ensembl_2.tmp >> human_ensembl_3.tmp ; }
			. qq{ sort human_ensembl_3.tmp > human_ensembl_via_perl_api.tab ; }
			. qq{ /bin/rm human_ensembl_[123].tmp ; }
		    ), "human_ensembl_via_perl_api.tab", $shouldUpdate, $genomeBuild);

    agnomenGetAnnot("Human_Ensembl_Mini_Test_First_10K_Rows", undef, 
		    , "hg19", "human_ensembl_via_perl_api.tab"
		    , ( qq{ head -n 10000 human_ensembl_via_perl_api.tab > minitest.tab }
		    ), "minitest.tab", $shouldUpdate, $genomeBuild);
    
    #my @col1 = @{readFileColumn($filename1, 0, $delim)};
    
    #my $annot1 = "AnnotatedFeatures.gff"; #"hb.bed"; #"human_ens_manual.bed"; #"b.gtf"; #"AnnotatedFeatures.gff"; #"annot.bed";
    #stderrPrint("Warning: using a manually-selected annotation file here.\n");

    if ($shouldUpdate) {
	    stderrPrint(colorString("green"));
	    stderrPrint("-------------------------------\n");
	    stderrPrint("[Finished updating!] Because we are UPDATING, we ignored any additional items specified on the command line.\n");
	    stderrPrint("-------------------------------\n");
	    stderrPrint(colorString("reset"));
    } else { 
	unless (defined($genomeBuild) && (length($genomeBuild) > 0)) { quitWithUsageError(">>> ERROR >>> it is MANDATORY to specify a genome build (example: hg19 or mm9). You do this with a command like this:   --species=hg19  or  -s hg19 \nCheck the annotation directory for a list of valid species.\n"); }

	if (scalar(keys(%GLOBAL_ANNOT_PATHS)) == 0) {
	    die "Uh oh, there were no annotation files, for some reason! Maybe you need to run 'agnomen-gene-annotate.pl --update' to re-download / re-generate them?\n";
	}
	
	my $numAnnotationsDone = 0;
	## For each file, we're going to annotate that file and write an output file. For each file, we're writing one output file per annotation type!
	## Therefore, if you have (say) 3 input files to annotate, and 4 types of annotation, you will end up with (3*4 = 12) final output files.
	for my $input (@inputFiles) {
	    if ($input !~ m/[.]bed$/i) {
		## Expecting the $input to be a BED file--- check to see if it has ".bed" as an extension?
		stderrPrint("Warning: the input file $input was EXPECTED to be a .bed file, but it appears not to have the .bed filename extension. Double-check to make sure this is the right file! Continuing on anyway, assuming that it is a bed file...\n");
	    }
	    ## 

	    while (my($annotShortName, $annotFile) = each(%GLOBAL_ANNOT_PATHS)) {
		
		## annotFile better exist, or we're in trouble!!
#	    if ($annotFile =~ m/.gz$/) {
#		$annotFile = "<(zcat $annotFile)"; # bash subshell
#		print STDERR "SUBSHELL\n";
#	    }
		if (!defined($databaseStr) || exists($databaseHash{lc($annotShortName)})) {
		    # Looks like we should include this annotation file!
		    my $tempFile = "1.agnomen.agw.in.progress." . int(rand(99999999)) . ".temp.tmp"; ## Temp file that is unlikely to already be in use!
		    
		    my $theOutputNameBase = ($longOutputNames) ? $input : basename($input); # see if we should use the FULL output name, with directory paths, or just the filename without the path. This is an option (--longnames)
		    my $agnomenOutputFile = qq{agnomen.$theOutputNameBase.$annotShortName.out.txt};
		    $agnomenOutputFile =~ s/[;:,\/\\]/_/g; ## Remove potentially "unsafe" characters from the output filename
		    my $cmd = qq{ $INTERSECT_BED_EXE -wao -a ${input} -b ${annotFile} > $tempFile } . qq{ && } . qq{ /bin/mv -f $tempFile $agnomenOutputFile};
		    stderrPrint(colorString("green"));
		    stderrPrint(qq{[ANNOTATING] using the data in $annotShortName\n});
		    stderrPrint(qq{[ANNOTATING]...\n    Running this command: $cmd\n});
		    stderrPrint(colorString("reset"));
		    my $intersectStatus = systemBash($cmd);
		    if ($intersectStatus != 0) {
			stderrPrint("Uh oh, non-zero exit status!");
		    } else {
			$numAnnotationsDone++;
		    }
		} else {
		    stderrPrint(colorString("cyan"));
		    stderrPrint(qq{[OMITTING] the annotation data in $annotShortName\n});
		    stderrPrint(colorString("reset"));
		}
	    }
	}

	if ($numAnnotationsDone == 0) {
	    stderrPrint(colorString("red"));
	    stderrPrint(qq{[FAILURE] It appears that we did not actually generate any annotation files. Check to see if perhaps the databases you have specified are invalid, or the species is not correct. Or maybe you need to run --update to generate the local databases.\n});
	    stderrPrint(colorString("reset"));
	} else {
	    stderrPrint(colorString("green"));
	    stderrPrint("-------------------------------\n");
	    stderrPrint("[Done!]\n");
	    stderrPrint("-------------------------------\n");
	    stderrPrint(colorString("reset"));
	}
    }

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

  --species = STRING  (Default: undefined. This is a REQUIRED parameter for annotation.)
    or: -s STRING
       Sets which genome build should be used for annotation. Example: hg19 (a human genome), mm9 (a mouse genome).
       Currently supported: hg19, mm9.
       Note that only ONE species may be specified here.

  --databases = STRING,OF,DATABASES (Default: undefined, meaning "use ALL databases for this species")
    or   --db = STRING,OF,DATABASES
       If this is specified, then we ONLY use certain databases for a species. Multiple databases for one species
       may be included here; the databases must be comma-delimited.
       You cannot specify databases for another species besides the one in --species (unless species is left blank, in which case
       *all* species are considered to be valid).
       You specify the database(s) by the SHORT NAME and not FILE NAME.
       SHORT NAME is defined in the perl code for agnomen-gene-annotate.pl; you cannot figure it out just by looking at a filename. 
       You have to check the code and look for the lines "agnomenGetAnnot(...)". The first argument to that function is the SHORT NAME.
       Example: specifying --databases=Hu_Ensembl,hg19Ribosome,hgExtra  is CORRECT
                       but --databases=human_ens.bed   <-------------   is WRONG (it uses a filename, not a short name)

  --longnames or --longname (Default: DISABLED)
       By default, the output names only include the filename: example: "agnomen.human.something.bed"
       If you want to include the entire path to the file, you can specify --longnames. Then you would end up with
       something like "agnomen._path_to_somewhere_on_filesystem.human.something.bed" . Much uglier, but may be useful if you have a
       lot of files with identical names in different directories.

  --annotdir = STRING (Default: Z_AGNOMEN_DATA)
    or: -a STRING
       A path on the filesystem where we keep the annotation data. The annotation data must be in species-specific
       sub-folders. Example filesystem:  MY_DATA/hg19/human_annot.bed
                                         MY_DATA/mm10/mouse_annot_1.bed
                                         MY_DATA/mm10/mouse_annot_extended.bed
       In the example above, you would specify --annotdir=MY_DATA

  --update
       Tells Agnomen to check for (and potentially update) newer annotation files. Does not fully work yet!
       By default, Agnomen will NOT update databases.
       * If you specify --update and NO species, we will attempt to download databases for ALL species.
       * If you specify --update and a particular species, we will only download databases for that one species.


KNOWN BUGS:

  None known.

TO DO:

  Add ???.

--------------

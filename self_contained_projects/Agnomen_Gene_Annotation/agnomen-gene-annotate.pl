#!/usr/bin/perl

use constant USE_COLORS_CONSTANT => 1; ## 1 = true, 0 = false
use POSIX      qw(ceil floor);
use List::Util qw(max min);
use Cwd;
if (USE_COLORS_CONSTANT) { use Term::ANSIColor; }

use File::Basename;
use Getopt::Long;

use Carp;
$SIG{ __DIE__ } = sub { Carp::confess( @_ ) }; # always print the stack trace when "die" activates

#no warnings 'numeric';

use strict;  use warnings;  use diagnostics;

sub main();

#"/work/Apps/Bio/bedtools/bin/intersectBed";
my $INTERSECT_BED_EXE_DEFAULT = "intersectBed";

my $SORT_FOR_INTERSECT_BED = " sort -t '\t' -k1,1 -k2,2n ";   # sorts a BED file by chromosome and then by start position in ascending order. Apparently end position is not important for intersectBed, based on the docs at: /work/Apps/Bio/bedtools/bin/intersectBed

my $DISCARD_FIRST_COLUMN = " cut -f 2- ";
my $HANDLE_CHR_PREFIX_AND_ALSO_SORT = " bed_and_gtf_handle_chr_prefix.pl "; # Alex's script! It's located in ~/TimeForScience/Lab_Tools/Perl/Tools/bed_and_gtf_handle_chr_prefix.pl

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
    # Requires: use Term::ANSIColor at the top of your file. Note: you have to
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

sub agwSafeColor($$) {
    my ($msg, $col) = @_;
    return (-t STDOUT && USE_COLORS_CONSTANT) ? Term::ANSIColor::colored("$msg", "$col") : $msg;
}

sub info($) {
    my ($msg) = @_;
    chomp($msg); # remove newline
    print STDOUT agwSafeColor("[INFO]: $msg\n", "cyan on_black");
}

sub bad($) {
    my ($msg) = @_;
    chomp($msg); # remove newline
    print STDOUT agwSafeColor("[ERROR]: $msg\n", "yellow on_red");
}

sub progressReport($) {
    my ($msg) = @_;
    chomp($msg); # remove newline
    print STDOUT agwSafeColor("[REPORT]: $msg\n", "green on_black");
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
    ## 3. $rawDownloadFileName -- the name to save the RAW file into. If it is UNDEFINED, then this will just be auto-set as the last part of the URL.
    ## 4. postprocessingCmd: something to do with the file after downloading it. Set it to "undef" if it isn't defined
    ## 5. finalName: the final filename after processing
    ## 6: allowedToUpdateFiles: whether we are allowed to actually modify the files on disk (to update them) or not
    ## 7: genomeBuildToUse: what genome build the user has specified for their input files. We will only add files that ALSO match this genome build.
    my ($shortName          , $url      , $theBuildForThisFile , $rawDownloadFileName
	, $postprocessingCmd, $finalName, $allowedToUpdateFiles, $genomeBuildToUse) = @_;

    (defined($shortName) && length($shortName) > 0) or die "Short name must be specified!";
    (defined($finalName) && length($finalName) > 0) or die "Final name must be defined!";

    (defined($allowedToUpdateFiles)) or die "allowed_to_update must be defined!";

    my $origWorkingDir = cwd();

    if (defined($genomeBuildToUse) && defined($theBuildForThisFile) && ($theBuildForThisFile ne $genomeBuildToUse)) {
	info("Not updating or adding the <$theBuildForThisFile> annotation file <$shortName> to the global annotation list, as it is for a different species.");
	return;
    }
    if (defined($genomeBuildToUse)) { info("Adding the annotation file <$shortName> to the global annotation list, as it is a valid annotation file for <$genomeBuildToUse>."); }

    my $localDir    = "${globalAnnotDir}/${theBuildForThisFile}";

    if (!$allowedToUpdateFiles) {
	info("[Skipping] re-updating of the file of <$finalName>: the --update flag was not passed into this script.");
    } else {
	if (!defined($rawDownloadFileName)) {
	    $rawDownloadFileName = File::Basename::basename($url);  # just get the LAST thing from the URL
	    info("Assuming that the unspecified download filename will be the last part of the url (setting it to <$rawDownloadFileName>)");
	}
	if (defined($url) && ($rawDownloadFileName ne File::Basename::basename($url))) {
	    bad("This is unusual; the local name of the file ($rawDownloadFileName) does not match the filename from the URL (the end of this URL: $url)!");
	}
	system("mkdir -p $localDir");
	chdir($localDir);
	if (defined($url) && $url) {
	    (not($url =~ /^(wget|curl)/)) or die "Your download url starts with 'wget' or 'curl'. It should start with ftp/http/something like that, instead. Don't put the command name (wget/curl) in the url!";
	    my $curlCmd = 'curl --remote-name ' . '"' . $url . '"' ; # downloads to whatever the remote name of the file was
	    #my $unzippedVersion = $rawDownloadFileName; $unzippedVersion    =~ s/[.](bz2|gz|zip)$//i;
	    if (agwFileHasContents("./$rawDownloadFileName")) { # changed directory ALREADY!
		progressReport("[Skipping] re-download of <$localDir/$rawDownloadFileName> from <$url> -- the downloaded file <$localDir/$rawDownloadFileName> already exists.\n");
	    } else {
		progressReport("Since we didn't find a file at <$localDir/$rawDownloadFileName>, we will start the download here...\n");
		progressReport("Downloading $url --> <$localDir/$rawDownloadFileName>...\n");
		agwSystemDieOnNonzero($curlCmd);
		progressReport("[Done] Successfully downloaded a file to this location <$localDir/$rawDownloadFileName>\n");
	    }
	}
    
	## Now the file should be downloaded, if necessary. We may also have to do a post-processing step, below.
	if (defined($postprocessingCmd) && $postprocessingCmd && length($postprocessingCmd) > 0) {
	    if (agwFileHasContents($finalName)) {
		progressReport("[Skipping] re-generation of $localDir/$finalName, as the output file already exists.\n");
	    } else {
		progressReport("Operating in local directory <$localDir>. Here is a list of all the files in this directory:\n" . join("   * ", `ls`) . "\n");
		progressReport("Now running this post-processing command:\n");
		progressReport("$postprocessingCmd\n");
		agwSystemDieOnNonzero($postprocessingCmd);
		if (not -e $finalName) { die "Uh oh, we failed to generate the file <$finalName>! Maybe the shell command had a mis-typed final filename, or perhaps something else went wrong?\n"; }
	    }
	} else {
	    # No postprocessing to do, so we just want to use the regular file as an annotation file
	}
    }

    chdir($origWorkingDir);
    addPathToGlobalAnnotationHash($shortName, "$localDir/$finalName"); # the actual full path to this file
}


sub handleUserSpecifiedDatabases($$) {
    my ($dataString, $verbose) = @_;
    # Input arg: delimited string (example: "this,is,a,string")
    # Return value: a hash, with all the elements from the string split up.
    # and: $referenceToHash: a reference (\%hash would be passed in) to a hash where we store the results.
    my %newHash = ();
    if (defined($dataString)) {
	($verbose) && info(qq{The user specified that only the following annotation databases should be used: $dataString\n});
	my @splitUpDb = split(/$DATABASE_INPUT_STRING_DELIM/, $dataString);
	for my $theItem (@splitUpDb) {
	    $newHash{$theItem} = 1; ## Let's remember that we SHOULD be using this database that the user specified on the command line!
	    $newHash{lc($theItem)} = 1; ## Let's also remember the LOWER CASE version of this name.
	}
	(scalar(keys(%newHash)) > 0) or die "The user specified some specific databases to use with the --databases (or -d) flag, but we couldn't interpret those as any actual valid databases! We ended up with ZERO databases. Note that this needs to be a comma-delimited string with NO SPACES in it.\n";
    } else {
	($verbose) && info(qq{Since no particular databases were specified, we will use ALL valid annotation databases.\n});
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
    my $intersectBedExeManuallySpecifiedLocation = undef;
    my $canAssumeSorted = 0;
    $Getopt::Long::passthrough = 1; # ignore arguments we don't recognize in GetOptions, and put them in @ARGV

    GetOptions("help|?|man" => sub { printUsageAndQuit(); }
	       , "delim|d=s" => \$delim
	       , "species|s=s" => \$genomeBuild
	       , "update!" => \$shouldUpdate
	       , "longnames|longname!" => \$longOutputNames
	       , "annotdir|a=s" => \$globalAnnotDir
	       , "databases|database|db=s" => \$databaseStr # comma-delimited string
	       , "assume-sorted!" => \$canAssumeSorted # assume sorted -by -k1,1 -k2,2n for intersectBed
	       , "ib=s" => \$intersectBedExeManuallySpecifiedLocation
	) or printUsageAndQuit();

    progressReport("-"x80 ."\n" . "Starting...!\n" . "-"x80 . "\n");

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
    if (defined($intersectBedExeManuallySpecifiedLocation) && (not -e $intersectBedExeManuallySpecifiedLocation)) { die "Error 95a: <intersectBed> was not found!\nWe  REQUIRE the program <intersectBed> to be installed. It is part of BedTools. We expected it to be located in <${intersectBedExeManuallySpecifiedLocation}>, but there was no file there at all! Check to make sure this path exists.\n"; }
    if (defined($intersectBedExeManuallySpecifiedLocation) && (not -x $intersectBedExeManuallySpecifiedLocation)) { die "Error 99x: <intersectBed> was not executable by this user!\nWe REQUIRE <intersectBed>, located in <${intersectBedExeManuallySpecifiedLocation}>, to be executable. However, the file at that location was not executable by this user! You may have to use chmod to fix this.\n"; }

    my $intersectBedExe = undef;
    if (defined($intersectBedExeManuallySpecifiedLocation)) {
	$intersectBedExe = $intersectBedExeManuallySpecifiedLocation;
    } else {
	$intersectBedExe = $INTERSECT_BED_EXE_DEFAULT;
    }
    
    my %databaseHash = handleUserSpecifiedDatabases($databaseStr, "verbose_messages");
    my $shouldUseAllDatabases = (!defined($databaseStr));

    progressReport("* Searching for annotation files in: <$globalAnnotDir>\n");
    if (defined($genomeBuild)) { progressReport("* Using the specified genome build: <$genomeBuild>\n"); }
    if (defined($genomeBuild)) { progressReport("* Searching for BED annotation files in the following file path: <$globalAnnotDir/$genomeBuild/>\n"); }

    if (scalar(@inputFiles) == 0 && !$shouldUpdate) {
	## Need to have some files to annotate!
	quitWithUsageError(">>> ERROR >>> We expected AT LEAST ONE valid file to annotate to be specified on the command line. We didn't get any, however.\n");
    }

    if ($shouldUpdate && scalar(@inputFiles) != 0) {
	quitWithUsageError(">>> ERROR >>> If you specify --update, you should NOT also specify any files to annotate! Try running again with JUST the --update flag and no filenames, or (if you already updated the files), remove the --update flag.\n");
    }

    my $DOLLAR_SIGN = '$';
    agnomenGetAnnot("Mouse_Ensembl", "ftp://ftp.ensembl.org/pub/release-68/gtf/mus_musculus/Mus_musculus.GRCm38.68.gtf.gz", "mm9", undef
		    , "zcat Mus_musculus.GRCm38.68.gtf.gz | $HANDLE_CHR_PREFIX_AND_ALSO_SORT > processed_Mus_musculus.GRCm38.68.gtf"
		    , "processed_Mus_musculus.GRCm38.68.gtf", $shouldUpdate, $genomeBuild);

    agnomenGetAnnot("Zebrafish_Ensembl", "ftp://ftp.ensembl.org/pub/release-70/gtf/danio_rerio/Danio_rerio.Zv9.70.gtf.gz", "danRer7", undef
		    , "zcat Danio_rerio.Zv9.70.gtf.gz | $HANDLE_CHR_PREFIX_AND_ALSO_SORT > processed_Danio_rerio.Zv9.70.gtf"
		    , "processed_Danio_rerio.Zv9.70.gtf", $shouldUpdate, $genomeBuild);

    agnomenGetAnnot("Chicken_Ensembl", "ftp://ftp.ensembl.org/pub/release-70/gtf/gallus_gallus/Gallus_gallus.WASHUC2.70.gtf.gz", "galGal3", undef
		    , "zcat Gallus_gallus.WASHUC2.70.gtf.gz | $HANDLE_CHR_PREFIX_AND_ALSO_SORT > processed_Gallus_gallus.WASHUC2.70.gtf"
		    , "processed_Gallus_gallus.WASHUC2.70.gtf", $shouldUpdate, $genomeBuild);

    # BedTools documentation is available at: http://bioinf.comav.upv.es/courses/sequence_analysis/bedtools.html
    agnomenGetAnnot("Human_Ensembl", undef, , "hg19", "human_ensembl_via_perl_api.tab"
		    , ( qq{ if [ ! -f human_ensembl_1_raw.tab.gz ]; then echo "\n>>>[NOTE] -- The ENSEMBL GRABBER takes about 18 hours to download data via Perl API! If you are reading this message, be prepared to wait QUITE A WHILE for the data to get downloaded!\n"; $AGNOMEN_ENSEMBL_GRABBER_SCRIPT | gzip > human_ensembl_1_raw.tab.gz ; fi; }
			. qq{ zcat human_ensembl_1_raw.tab.gz | grep -v '^[#]' | grep -v '^$DOLLAR_SIGN' | $HANDLE_CHR_PREFIX_AND_ALSO_SORT | uniq > processed_human_ensembl_via_perl_api.tab ; } # sorts a BED file by chromosome and then by start position in ascending order. Apparently end position is not important for intersectBed, based on the docs at: /work/Apps/Bio/bedtools/bin/intersectBed
		    ), "processed_human_ensembl_via_perl_api.tab", $shouldUpdate, $genomeBuild);

    agnomenGetAnnot("Human_Ensembl_Mini_Test_First_10K_Rows", undef,
		    , "hg19", "human_ensembl_via_perl_api.tab"
		    , ( qq{ head -n 10000 human_ensembl_via_perl_api.tab | $HANDLE_CHR_PREFIX_AND_ALSO_SORT > minitest.tab }
		    ), "minitest.tab", $shouldUpdate, $genomeBuild);


    # hg19: Phast Conservation 46-way placental
    # columns 2 and onward are the only ones that are interesting ones from the genome browser -- discard the FIRST column!
    agnomenGetAnnot("hg19_phast_46_placental_conservation", "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/phastCons46wayPlacental.txt.gz", "hg19", undef
		    , qq{zcat phastCons46wayPlacental.txt.gz | cut -f 2-8,10- | $HANDLE_CHR_PREFIX_AND_ALSO_SORT  > browser_phastCons46wayPlacental.bed}
		    , qq{browser_phastCons46wayPlacental.bed}, $shouldUpdate, $genomeBuild);

    # hg19: conservation track: phyloP 46-way all-species conservation
    # columns 2 and onward are the only ones that are interesting ones from the genome browser -- discard the FIRST column!
    agnomenGetAnnot("hg19_phylo_46_all_species_conservation", "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/phyloP46wayAll.txt.gz", "hg19", undef
		    , qq{zcat phyloP46wayAll.txt.gz | cut -f 2-8,10- | $HANDLE_CHR_PREFIX_AND_ALSO_SORT  > browser_phyloP46wayAll.bed}
		    , "browser_phyloP46wayAll.bed", $shouldUpdate, $genomeBuild);
    #"hg19 conservation track: PhastCons 46-way all-species conservation", "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/phastCons46way.txt.gz"
    
    #"GNF expression atlas", "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/gnfAtlas2.txt.gz"

    # "Human (hg19): Broad Institute chromatin HMM for Gm12878 cell line"
    agnomenGetAnnot("hg19_broad_chromatin_gm12878", "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/wgEncodeBroadHmmGm12878HMM.txt.gz", "hg19", undef
		    , qq{zcat wgEncodeBroadHmmGm12878HMM.txt.gz | ${DISCARD_FIRST_COLUMN} | $HANDLE_CHR_PREFIX_AND_ALSO_SORT  > browser_wgEncodeBroadHmmGm12878HMM.bed}, "browser_wgEncodeBroadHmmGm12878HMM.bed", $shouldUpdate, $genomeBuild);

    #"Human (hg19): Broad Institute chromatin HMM for H1 HESC"
    agnomenGetAnnot("hg19_broad_chromatin_h1hesc", "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/wgEncodeBroadHmmH1hescHMM.txt.gz", "hg19", undef
		    , qq{zcat wgEncodeBroadHmmH1hescHMM.txt.gz | ${DISCARD_FIRST_COLUMN} | $HANDLE_CHR_PREFIX_AND_ALSO_SORT  > browser_wgEncodeBroadHmmH1hescHMM.bed}, "browser_wgEncodeBroadHmmH1hescHMM.bed", $shouldUpdate, $genomeBuild);
    
    # "Human (hg19): Broad Institute chromatin HMM for Hepg2 cell line"
    agnomenGetAnnot("hg19_broad_chromatin_hepg2", "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/wgEncodeBroadHmmHepg2HMM.txt.gz", "hg19", undef
		    , qq{zcat wgEncodeBroadHmmHepg2HMM.txt.gz | ${DISCARD_FIRST_COLUMN} | $HANDLE_CHR_PREFIX_AND_ALSO_SORT  > browser_wgEncodeBroadHmmHepg2HMM.bed}, "browser_wgEncodeBroadHmmHepg2HMM.bed", $shouldUpdate, $genomeBuild);

    # "Human (hg19): Broad Institute chromatin HMM for Hmec cell line"
    agnomenGetAnnot("hg19_broad_chromatin_hmec", "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/wgEncodeBroadHmmHmecHMM.txt.gz", "hg19", undef
		    , qq{zcat wgEncodeBroadHmmHmecHMM.txt.gz | ${DISCARD_FIRST_COLUMN} | $HANDLE_CHR_PREFIX_AND_ALSO_SORT  > browser_wgEncodeBroadHmmHmecHMM.bed}, "browser_wgEncodeBroadHmmHmecHMM.bed", $shouldUpdate, $genomeBuild);

    # "Human (hg19): Broad Institute chromatin HMM for Hsmm cell line"
    agnomenGetAnnot("hg19_broad_chromatin_hsmm", "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/wgEncodeBroadHmmHsmmHMM.txt.gz", "hg19", undef
		    , qq{zcat wgEncodeBroadHmmHsmmHMM.txt.gz | ${DISCARD_FIRST_COLUMN} | $HANDLE_CHR_PREFIX_AND_ALSO_SORT  > browser_wgEncodeBroadHmmHsmmHMM.bed}, "browser_wgEncodeBroadHmmHsmmHMM.bed", $shouldUpdate, $genomeBuild);

    # "Human (hg19): Broad Institute chromatin HMM for Huvec cell line"
    agnomenGetAnnot("hg19_broad_chromatin_huvec", "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/wgEncodeBroadHmmHuvecHMM.txt.gz", "hg19", undef
		    , qq{zcat wgEncodeBroadHmmHuvecHMM.txt.gz | ${DISCARD_FIRST_COLUMN} | $HANDLE_CHR_PREFIX_AND_ALSO_SORT  > browser_wgEncodeBroadHmmHuvecHMM.bed}, "browser_wgEncodeBroadHmmHuvecHMM.bed", $shouldUpdate, $genomeBuild);

    # "Human (hg19): Broad Institute chromatin HMM for K562 cell line"
    agnomenGetAnnot("hg19_broad_chromatin_k562", "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/wgEncodeBroadHmmK562HMM.txt.gz", "hg19", undef
		    , qq{zcat wgEncodeBroadHmmK562HMM.txt.gz  | ${DISCARD_FIRST_COLUMN} | $HANDLE_CHR_PREFIX_AND_ALSO_SORT  > browser_wgEncodeBroadHmmK562HMM.bed}, "browser_wgEncodeBroadHmmK562HMM.bed", $shouldUpdate, $genomeBuild);

    # "Human (hg19): Broad Institute chromatin HMM for Nhek cell line"
    agnomenGetAnnot("hg19_broad_chromatin_nhek", "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/wgEncodeBroadHmmNhekHMM.txt.gz", "hg19", undef
		    , qq{zcat wgEncodeBroadHmmNhekHMM.txt.gz | ${DISCARD_FIRST_COLUMN} | $HANDLE_CHR_PREFIX_AND_ALSO_SORT  > browser_wgEncodeBroadHmmNhekHMM.bed}, "browser_wgEncodeBroadHmmNhekHMM.bed", $shouldUpdate, $genomeBuild);

    # "Human (hg19): Broad Institute chromatin HMM for Nhlf cell line"
    agnomenGetAnnot("hg19_broad_chromatin_nhlf", "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/wgEncodeBroadHmmNhlfHMM.txt.gz", "hg19", undef
		    , qq{zcat wgEncodeBroadHmmNhlfHMM.txt.gz | ${DISCARD_FIRST_COLUMN} | $HANDLE_CHR_PREFIX_AND_ALSO_SORT  > browser_wgEncodeBroadHmmNhlfHMM.bed}, "browser_wgEncodeBroadHmmNhlfHMM.bed", $shouldUpdate, $genomeBuild);
	

    #my @col1 = @{readFileColumn($filename1, 0, $delim)};
    
    #my $annot1 = "AnnotatedFeatures.gff"; #"hb.bed"; #"human_ens_manual.bed"; #"b.gtf"; #"AnnotatedFeatures.gff"; #"annot.bed";
    #progressReport("Warning: using a manually-selected annotation file here.\n");

    if ($shouldUpdate) {
	    progressReport("-------------------------------\n");
	    progressReport("[Finished updating!] Because we are UPDATING, we ignored any additional items specified on the command line.\n");
	    progressReport("-------------------------------\n");
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
		bad("Warning: the input file $input was EXPECTED to be a .bed file, but it appears not to have the .bed filename extension. Double-check to make sure this is the right file! Continuing on anyway, assuming that it is a bed file...\n");
	    }

	    while (my($annotShortName, $annotFile) = each(%GLOBAL_ANNOT_PATHS)) {
		if ($shouldUseAllDatabases || exists($databaseHash{lc($annotShortName)})) {
		    # Looks like we should include this annotation file!
		    my $tempFile = "1.agnomen.agw.in.progress." . int(rand(99999999)) . ".temp.tmp"; ## Temp file that is unlikely to already be in use!
		    my $theOutputNameBase = ($longOutputNames) ? $input : basename($input); # see if we should use the FULL output name, with directory paths, or just the filename without the path. This is an option (--longnames)
		    my $agnomenOutputFile = qq{agnomen.$theOutputNameBase.$annotShortName.out.txt};
		    $agnomenOutputFile =~ s/[;:,\/\\]/_/g; ## Remove potentially "unsafe" characters from the output filename

		    my $SORT_OPTION = $canAssumeSorted ? " -sorted " : " "; # bedIntersect has a '-sorted' option that can be used if BOTH input files have been sorted -by -k1,1 -k2,2n
		    my $cmd = qq{ ${intersectBedExe} $SORT_OPTION -wao -a ${input} -b ${annotFile} > $tempFile } . qq{ && } . qq{ /bin/mv -f $tempFile $agnomenOutputFile}; # <-- we WRITE to a temp file so that when the file actually moves to the FINAL location, it is guaranteed to be 100% finished
		    progressReport(qq{[ANNOTATING] using the data in $annotShortName\n});
		    progressReport(qq{[ANNOTATING]...\n    Running this command: $cmd\n});
		    my $intersectStatus = systemBash($cmd);
		    if ($intersectStatus != 0) {
			bad("Uh oh, non-zero exit status!");
		    } else {
			$numAnnotationsDone++;
		    }
		} else {
		    info(qq{OMITTING the annotation data in <$annotShortName>, as that database was not specified in the database string ($databaseStr).\n});
		}
	    }
	}

	if ($numAnnotationsDone == 0) {
	    bad(qq{Failure: It appears that we did not actually generate any annotation files. Check to see if perhaps the databases you have specified are invalid, or the species is not correct. Or maybe you need to run --update to generate the local databases.\n});
	} else {
	    progressReport("-------------------------------\n");
	    progressReport("[Done!]\n");
	    progressReport("-------------------------------\n");
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

  --ib=PATH/TO/intersectBed (Default: looks on user search $PATH only)
       Path to the 'intersectBed' tool. If you do not specify this option, then we try to run the command "intersectBed" with
       no directory prefix. Example: --ib=/common/share/bin/64bit/intersectBed

KNOWN BUGS:

  None known.

TO DO:

  Add ???.

--------------

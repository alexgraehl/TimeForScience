#!/usr/bin/perl

# ===================================================================================================================================
# Limited How-to tutorials are at: http://uswest.ensembl.org/info/docs/api/core/core_tutorial.html?redirect=mirror;source=www.ensembl.org
## FULL Ensembl Documentation: http://www.ensembl.org/info/docs/Doxygen/core-api/classBio_1_1EnsEMBL_1_1Transcript.html

# ===================================================================================================================================

# How to install:
# Read instructions in: http://uswest.ensembl.org/info/docs/api/api_installation.html?redirect=mirror;source=www.ensembl.org
# Download Ensembl perl API from:  http://uswest.ensembl.org/info/docs/api/api_cvs.html
# Download BioPerl from:           git clone git://github.com/bioperl/bioperl-live.git

### You will need the Perl "DBD::mysql" module to be installed for this to work!
### Here's how you can do it:  1) perl -MCPAN -we 'shell'
###                            2) install DBD::mysql

## Takes about 18 hours to download a mammalian annotation.

## Make sure to output this to STDOUT or something!

#
# Finally, filter the resulting file: grep -v '^[#]' HUMAN_CORE_ENSEMBL_HUGE | grep -v '^$'   | sort -k 2,2 -k 3,3g | uniq > HUMAN_FILTERED_1
# cut -f 2- HUMAN_FILTERED_1 > a.tmp; paste.pl -do '' 'chr' a.tmp > HUMAN_FILTERED_2; cat a.tmp >> HUMAN_FILTERED_2 ; sort HUMAN_FILTERED_2 > HUMAN_FILTERED_3




### =============================================================
### =============================================================
## You need to include some links to all the ENSEMBL Perl libraries, PLUS the bioperl libraries. Note their location on the filesystem below:
use lib "/work/Apps/perl/ensembl_perl_api/bioperl-live";
use lib "/work/Apps/perl/ensembl_perl_api/ensembl/modules";
use lib "/work/Apps/perl/ensembl_perl_api/ensembl-compara/modules";
use lib "/work/Apps/perl/ensembl_perl_api/ensembl-variation/modules";
use lib "/work/Apps/perl/ensembl_perl_api/ensembl-functgenomics/modules";
### =============================================================
### =============================================================

use strict;  use warnings;  use diagnostics;
use Bio::EnsEMBL::Registry;
use List::Util qw(min max);
use Getopt::Long;

#use Scalar::Util 'looks_like_number';  ## <-- Perl is the worst language sometimes. What an elegant function!!!
$| = 1; # flush output immediately to stdout. Otherwise you have to wait for.... something! It's annoying, that's for sure.
my $registry = 'Bio::EnsEMBL::Registry';
my $GLOBAL_BIOTYPE_PREFIX = "ts_biotype:";
my $DEFAULT_SCORE_FOR_THINGS_WITHOUT_A_SCORE = 0;

my $ALEX_ANNOTATION_DELIMITER = "|";
my $DEBUG_QUIT_AFTER_A_FEW_GENES = 0; # default is false
my $adaptorSpecies = "Human"; # default is human. Probably shouldn't be, though!

$Getopt::Long::passthrough = 1; # ignore arguments we don't recognize in GetOptions, and put them in @ARGV
GetOptions("help|?|man" => sub { die "Sorry no help is available! Note that the only option to this script right now is --species=SOMETHING."; }
	   , "debug!" => \$DEBUG_QUIT_AFTER_A_FEW_GENES
	   , "species|s=s" => \$adaptorSpecies
    ) or die "no options / wrong options given to this script (see the GetOptions part of the code";


if ($DEBUG_QUIT_AFTER_A_FEW_GENES) { print STDERR "\n\n### WARNING: QUITTING AFTER ONLY A FEW GENES! This is for DEBUGGING PURPOSES!\n\n"; }









print STDERR "### [Step 1] Loading registry (this takes 15-60 seconds). " . `date` . "\n";
## Registry info can be found at: http://uswest.ensembl.org/info/data/mysql.html?redirect=mirror;source=www.ensembl.org
$registry->load_registry_from_db( '-host'   => 'useastdb.ensembl.org' #'ensembldb.ensembl.org', #'ensembldb.ensembl.org',
				  , '-port' => '5306' # 5306 on useastdb means: Current ensembl version with Previous ensembl version ONLY
                                  , '-user' => 'anonymous');
#                                  , '-db_version' => '68' );

# http://plants.ensembl.org/info/docs/api/variation/variation_tutorial.html
# You will need to have a registry configuration file set up. By default, it takes the file defined by the ENSEMBL_REGISTRY environment variable or the file named .ensembl_init in your home directory if the former is not found. Additionally, it is possible to use a specific file (see perldoc Bio::EnsEMBL::Registry or later in this document for some examples on how to use a different file). An example of such file can be found in ensembl/modules/Bio/EnsEMBL/Utils/ensembl_init.example, and below you have a slightly modified copy of it.



foreach my $dba (@{$registry{'_DBA'}}){
    print STDERR "ZOG: " . $dba->species() . "\n"; 
    #$species_hash{lc($dba->species())} = 1;
}








my @adaps = @{Bio::EnsEMBL::Registry->get_all_adaptors()};
print STDERR ">>>>>>>>>> NUM ADAPTORS: " . scalar(@adaps) . "\n";
my $aref = Bio::EnsEMBL::Registry->get_all_aliases('Homo sapiens');
print STDERR ">>>>>>>>> GETTING ALL ALIASES\n";
for my $a (@$aref) {
    print STDERR ">>>>>>>>>>>> ALIAS IS: $a \n";
}
print STDERR ">>>>>>>>> DONE.\n";

my $alias_exists = Bio::EnsEMBL::Registry->alias_exists($adaptorSpecies);
print STDERR ">>>>>>>>>>>>>>>>>>>>>>>>\n\n" . $adaptorSpecies . "\n>>>>>>>>>>>>>>>>>>>>>>>\n";
print STDERR ">>>>>>>>>>>>>>>>>>>>>>>>\n\n" . $alias_exists . "\n>>>>>>>>>>>>>>>>>>>>>>>\n";



print STDERR "### [Step 2] Getting adaptor: " . `date` . "\n";

#my $transcript_adaptor = $registry->get_adaptor( $adaptorSpecies, 'Core', 'Transcript' );
my $slice_adaptor = $registry->get_adaptor($adaptorSpecies, 'Core', 'Slice');
#my $tr_adaptor    = $registry->get_adaptor($adaptorSpecies, 'Core', 'Transcript' );
my $gene_adaptor  = $registry->get_adaptor($adaptorSpecies, 'Core', 'Gene');

my @allSlices = @{ $slice_adaptor->fetch_all('chromosome') };
#my $slice = $slice_adaptor->fetch_by_region( 'chromosome', '20', 10e6, 10.3e6);

print STDERR "### [Step 3] Fetching... " . `date` . "\n";

#ENSG00000139597
#my $someTestGene = $gene_adaptor->fetch_by_stable_id('ENSG00000132639'); # just pick this one out of nowhere


############ ============= PROTOTYPES ============ ##############
sub featureToBedString($$$;$);
sub manualToBedString($$$$$$$$$$);
sub featureToBedStringWithManualPositions($$$$$;$);
############ ============= PROTOTYPES ============ ##############


sub intersectRegions($$$$) {
    ## Tells you the intersectinon of two regions
    ## AAAAAAAAAAAAAAA <-- region A
    ##         BBBBBBBBBBBBBBBB <-- region B
    ##         CCCCCCC <-- result

    ## Returns the LEFT and RIGHT coordinates as a pair ($left, $right);
    ## Returns (undef, undef) if there is NO INTERSECTION

    my ($aleft, $aright, $bleft, $bright) = @_;
    
    my $newLeft = max($aleft, $bleft);
    my $newRight = min($aright, $bright);

    if ($newLeft > $newRight) {
	return(undef, undef); # if there is NO OVERLAP, then just return undef, undef
    } else {
	return($newLeft, $newRight);
    }

}

sub truncateGeneRegion($$$$) {
    # Takes an input region, defined as masterLeft and masterRight,
    # and a "culling" region (cullLeft to cullRight)
    # and subtracts out the cullLeft---cullRight region from the master region.
    # Example:
    ## -------MMMMMMMMMMMMMMM     <-- master region
    ## ----CCCCCCCCCCCC------     <-- Culling region
    ##        ---------MMMMMM     <- result!
    ## FAILS if the culling region is a strict subset of the master region (but they can be equal).

    # This function is used to remove UTRs from exons, leaving only the translated region behind (hopefully).

    my ($masterLeft, $masterRight, $cullLeft, $cullRight) = @_;
    # The "master" input is the locations we want to trim, based on the "cull". So we "clip" the region
    # between cullLeft and cullRight from masterLeft and masterRight
    # Returns the new clipped region.

    if ($DEBUG_QUIT_AFTER_A_FEW_GENES) { print STDERR "Truncating $masterLeft,$masterRight by $cullLeft,$cullRight...\n"; }
    ## MMMMMMMMMMMMM
    ##        CCCCCCCCCCCC
    ## RRRRRRR           <-- and this should be the culled result

    ##       MMMMMMMMMMMMM
    ## CCCCCCCCCCCC
    ##             RRRRRRR   <-- and this should be the culled result

    ##    MMMMMMM
    ## CCCCCCCCCCCCCCC
    ##  (nothing)           <-- and this should be the culled result

    ## MMMMMMMMMMMMMMMMM
    ##     CCCCCCCCCC
    ## gives:
    ##     XXXXXXXXXX  <-- this is sort of a weird and invalid case, as it splits up the master region in two! Disallowed!

    if ($cullRight < $masterLeft || $cullLeft > $masterRight) {
	# No culling occurrs
	##          MMMMMMMMMMM
	## CCCCC
	## or:
	## MMMMMMMMMMM
	##                CCCCC
    } else {
	if ($cullRight == $masterRight) {
	    if ($cullLeft == $masterLeft || $cullLeft < $masterLeft) {
		## MMMMMMMMMMM
		## CCCCCCCCCCC
		$masterRight = undef;
		$masterLeft = undef; ## Ok, we culled the ENTIRE region
		
	    } elsif ($cullLeft > $masterLeft) {
		##      MMMMMMMMMMM
		##           CCCCCC
		##      NNNNN
		$masterRight = ($cullLeft-1);

	    } else { die "not possible to get here A"; }
	} elsif ($cullRight < $masterRight) {
	    if ($cullLeft == $masterLeft || $cullLeft < $masterLeft) {
		## MMMMMMMMMM
		## CCCCCC
		##       NNNN
		$masterLeft = ($cullRight+1);
	    } elsif ($cullLeft > $masterLeft) {
		## MMMMMMMMMM
		##   CCCCCC
		## MM------MM <-- can't handle complicated situations like this! NO CULLING will occur, instead.
		print STDERR "### CULLING TYPE B3: This is NOT a culling situation we expect to arise: [$masterLeft,$masterRight] culled by [$cullLeft,$cullRight]. NO CULLING PERFORMED.\n";
		
	    } else { die "not possible to get here B"; }
	} elsif ($cullRight > $masterRight) {
	    if ($cullLeft == $masterLeft || $cullLeft < $masterLeft) {
		## MMMMMMMMMM
		## CCCCCCCCCCCCCCC
		## (nothing left)
		$masterRight = undef;
		$masterLeft = undef; ## Ok, we culled the ENTIRE region
	    } elsif ($cullLeft > $masterLeft) {
		##  MMMMMMMMMM
		##        CCCCCCCCCCCCCCC
		##  NNNNNN
		$masterRight = ($cullLeft-1)
	    } else { die "not possible to get here C"; }
	} else {
	    die "Not possible to get here";
	}
    }

    #if ($DEBUG_QUIT_AFTER_A_FEW_GENES) { print STDERR "### result is: $masterLeft,$masterRight\n"; }
    return($masterLeft, $masterRight);
}




sub utrStringCalculate($$$$$$$) {
    my ($utr, $ex, $allExonsArrayPtr, $gene, $transcript, $fiveOrThree, $theTranscriptText) = @_;
    ## Returns a string of all the awesome UTR stuff to print, given some input items.
    ## Note that you pass in an ARRAY for $allExonsArrayPtr -- the perl prototype (\@) handles the automatic de-referencing. Perl is a nightmare language sent to punish us for all our misdeeds!!! If you don't believe me, check: http://mail.pm.org/pipermail/wellington-pm/2004-November/000211.html
    my $utrStr = '';
    if (!defined($utr)) { 
	$utrStr = "### No ${fiveOrThree}' UTR for transcript " . $transcript->stable_id() . ".\n";
    } else {
	$utrStr = (featureToBedString($utr, $gene, $transcript, "UTR${fiveOrThree}_SPANNING_WITHOUT_BREAKS") 
		   . "\t" . "ATYPE:UTR${fiveOrThree}_SPANNING_WITHOUT_BREAKS${ALEX_ANNOTATION_DELIMITER}${theTranscriptText}\n" );
	# Now we're going to go through each exon and see if it intersects the UTR---if it DOES, then we're going to output that as a UTR5_FRAGMENTED
	my $numFragmentsSoFar = 0;
	for (my $e = 0; $e < scalar(@{$allExonsArrayPtr}); $e++) {
	    my $exCheck = ${$allExonsArrayPtr}[$e];
	    my ($utrOverlapLeft, $utrOverlapRight) = intersectRegions($exCheck->start(), $exCheck->end(), $utr->start(), $utr->end());
	    if (defined($utrOverlapLeft) && defined($utrOverlapRight)) {
		# Guess we'd better output this "overlap" region as being the UTR!
		$numFragmentsSoFar++;
		my $newUTRLocation = (featureToBedStringWithManualPositions($utr, $gene, $transcript, $utrOverlapLeft, $utrOverlapRight, "UTR${fiveOrThree}_FRAGMENT_${numFragmentsSoFar}") 
				      . "\t" . "ATYPE:UTR${fiveOrThree}_FRAGMENT_${numFragmentsSoFar}${ALEX_ANNOTATION_DELIMITER}${theTranscriptText}\n" );
		$utrStr .= $newUTRLocation; ## append it to the $utr string
	    }
	}
    }
    return $utrStr;
}



my $numTscriptsTOTAL = 0;
my $numGenesTOTAL = 0;
foreach my $slice (@allSlices) {
    my $genesPtr = $gene_adaptor->fetch_all_by_Slice($slice);
    my @allGenes = @{$genesPtr}; # <-- de-reference a pointer to an array
    #my @allGenes = ($someTestGene);
    foreach my $gene (@allGenes) {
	$numGenesTOTAL++;
	my @allTranscripts = @{ $gene->get_all_Transcripts() };	
	my $nTransSoFarThisGeneOnly = 0;
	foreach my $transcript (@allTranscripts) {
	    $nTransSoFarThisGeneOnly++;
	    $numTscriptsTOTAL++;
	    my $tscriptNumOfTotal = $nTransSoFarThisGeneOnly . "/" . scalar(@allTranscripts);
	    my $inTsText = "TSCRIPT_" . $tscriptNumOfTotal; ## the "IN_TSCRIPT_4/15" text
	    print STDOUT "###\n########################################## TRANSCRIPT #############################################\n";
	    print featureToBedString($gene, $gene, "NA", undef) . "\tATYPE:GENE\n";
	    print featureToBedString($transcript, $gene, $transcript, undef) . "\tATYPE:TSCRIPT_" . $tscriptNumOfTotal . "\n";
	    my @introns = @{ $transcript->get_all_Introns() };
	    my @exons   = @{ $transcript->get_all_Exons() };
	    my $FIRST_EXON_INDEX = 0;
	    my $LAST_EXON_INDEX = (scalar(@exons)-1);
	    
	    my $utr5seq = defined($transcript->five_prime_utr()) ? $transcript->five_prime_utr()->seq() : undef;
	    my $utr3seq = defined($transcript->three_prime_utr()) ? $transcript->three_prime_utr()->seq() : undef;

	    ## Q: Why are there the awkward "try / catch" blocks below?
	    ## A:
	    ## We were getting this error, like 90% of the way through Ensembl:
	    # MSG: Start (2654896) must be less than or equal to end+1 (5510)
	    # STACK Bio::EnsEMBL::Feature::new /work/Apps/perl/ensembl_perl_api/ensembl/modules/Bio/EnsEMBL/Feature.pm:141
	    # STACK Bio::EnsEMBL::Transcript::three_prime_utr_Feature /work/Apps/perl/ensembl_perl_api/ensembl/modules/Bio/EnsEMBL/Transcript.pm:1772
	    # STACK toplevel this_script.pl:266
	    ## What this means is that SOMETHING is wrong with getting five_prime_utr_Feature() from a maybe-invalid feature... or something.
	    ## Basically, I think the annotation is wrong for something, and the code really should never actually result in an error from inquiring about a feature.
	    ## So this is a workaround! If the UTR feature is screwed up, we just totally ignore it.

	    my $utr5 = undef;
	    eval { # "try"
		$utr5 = defined($utr5seq) ? $transcript->five_prime_utr_Feature() : undef;
		1;
	    } or do { # "catch exception"
		$utr5seq = undef;
		$utr5 = undef;
	    };
	    
	    my $utr3 = undef;
	    eval { # "try"
		$utr3 = defined($utr3seq) ? $transcript->three_prime_utr_Feature() : undef;
		1;
	    } or do { # "catch exception" -- keep this! it's important (read above for description)
		$utr3seq = undef;
		$utr3 = undef;
	    };

	    my $utr5seqlen = (defined($utr5seq)) ? length($utr5seq) : 0;
	    my $utr3seqlen = (defined($utr3seq)) ? length($utr3seq) : 0;

	    my $previousExonStart = undef; # for debuggin ONLY

	    for (my $iii = 0; $iii < scalar(@exons); $iii++) {
		my $inExonText = "EXON_" . ($iii+1) . "/" . scalar(@exons);
		my $IS_FORWARD_STRAND = ($transcript->strand() > 0);
		my $ex = $exons[$iii];
		my $isFirstExon = (0 == $iii);
		my $isLastExon  = (scalar(@exons)-1 == $iii);

		if (defined($previousExonStart) && $IS_FORWARD_STRAND && $previousExonStart > $ex->start()) {
		    die "Bug--positive strand exons are expected to ALWAYS fetched from Ensembl in ***increasing*** order. This one had exon start of: ". $ex->start() . " versus a previous start of $previousExonStart.\n";
		}
		if (defined($previousExonStart) && !$IS_FORWARD_STRAND && $previousExonStart < $ex->start()) {
		    die "Bug--negative strand exons are expected to ALWAYS fetched from Ensembl in ***decreasing*** order. This one had exon start of: ". $ex->start() . " versus a previous start of $previousExonStart.\n";
		}
		$previousExonStart = $ex->start();
		
		my $geneCommonNameIfApplicable = (($gene->can('external_name')) ? $gene->external_name() : 'NA');

		my $EMPTY = '';
		my $utr5Str = ($isFirstExon) ? utrStringCalculate($utr5, $ex, \@exons, $gene, $transcript, '5', $inTsText) : $EMPTY; ## <-- only calculate (and print) this for the FIRST exon!
		my $utr3Str = ($isLastExon)  ? utrStringCalculate($utr3, $ex, \@exons, $gene, $transcript, '3', $inTsText) : $EMPTY; ## <-- only calculate (and print) this for the LAST exon!

		my ($exonNoUtrLeft, $exonNoUtrRight) = ($ex->start(), $ex->end());
		if (defined($utr3)) {
		    ($exonNoUtrLeft, $exonNoUtrRight) = truncateGeneRegion($exonNoUtrLeft, $exonNoUtrRight, $utr3->start(), $utr3->end());
		}
		if (defined($utr5) && defined($exonNoUtrLeft) && defined($exonNoUtrRight)) {
		    ($exonNoUtrLeft, $exonNoUtrRight) = truncateGeneRegion($exonNoUtrLeft, $exonNoUtrRight, $utr5->start(), $utr5->end());
		}
		
		my $trimmingOccurred = (!defined($exonNoUtrLeft) || !defined($exonNoUtrRight) || $ex->start() != $exonNoUtrLeft || $ex->end() != $exonNoUtrRight); ## some start/end location is now different! Thus, trimming actually occurred for this exon.
		if ($DEBUG_QUIT_AFTER_A_FEW_GENES) { print STDERR "### Trimming status: " . $trimmingOccurred . "\n"; }

		my $exonMinusUtrStr = $EMPTY;
		if ($trimmingOccurred) {
		    my $trimmedExonName = $ex->stable_id() . $ALEX_ANNOTATION_DELIMITER . ((defined($exonNoUtrLeft) && defined($exonNoUtrRight)) ? "WITHOUT_UTR" : "TRIMMED_TO_NOTHING");
		    my $maybeCommentThisLine = ((defined($exonNoUtrLeft) && defined($exonNoUtrRight)) ? '' : '### '); ## Maybe comment out this exon (with '###' at the front) if NOTHING is left after clipping.
		    $exonMinusUtrStr = $maybeCommentThisLine .  featureToBedStringWithManualPositions($ex, $gene, $transcript, $exonNoUtrLeft, $exonNoUtrRight, $trimmedExonName) . "\tATYPE:UTRLESS_EXON_" . $inExonText . ${ALEX_ANNOTATION_DELIMITER} . $inTsText . $ALEX_ANNOTATION_DELIMITER . "WITHOUT_UTR" . "\n";
		}

		my $utrOverlapNoteWithDelim = (!$trimmingOccurred) ? $EMPTY : ($ALEX_ANNOTATION_DELIMITER . "INCLUDES_UTR");
		my $exonStr = featureToBedString($ex, $gene, $transcript, undef) . "\tATYPE:" . $inExonText . ${ALEX_ANNOTATION_DELIMITER} . $inTsText . $utrOverlapNoteWithDelim . "\n";

		my $thisIntron = shift(@introns); #($IS_FORWARD_STRAND) ? shift(@introns) : pop(@introns); ## Grab off the START (shift) if we're looking at a forward strand item, or grab off the END (pop) if we're looking at a reverse-strand item
		my $intronStr = (!defined($thisIntron)) ? $EMPTY : featureToBedString($thisIntron, $gene, "INTRON", "INTRON") . "\tATYPE:INTRON_AFTER_EXON_" . ($iii+1) . ${ALEX_ANNOTATION_DELIMITER} . $inTsText . "\n";

		if ($IS_FORWARD_STRAND) {
		    ## write 5' utr, then exon, then possible intron, then 3' utr
		    print $utr5Str;
		    print $exonStr;
		    print $exonMinusUtrStr;
		    print $intronStr;
		    print $utr3Str;
		} else {
		    ## write 3' utr, then (optional) intron, then exon, then 5' utr
		    print $utr3Str;
		    print $intronStr;
		    print $exonStr;
		    print $exonMinusUtrStr;
		    print $utr5Str;
		}

	    }
	}
	if ($DEBUG_QUIT_AFTER_A_FEW_GENES && $numGenesTOTAL >= 7) { print "### DEBUGGING: QUITTING AFTER A FEW GENES.\n"; exit(0); }
    }
}

# The protein sequence is obtained from the translate() method. If the transcript is non-coding, undef is returned.
#my $protein = $transcript->translate();
#print "Translation: ", ( defined $protein ? $protein->seq() : "None" ), "\n";


sub featureToBedString($$$;$) {
    my ($feature, $gene, $tscript, $nameOverride) = @_;
    return featureToBedStringWithManualPositions($feature, $gene, $tscript, $feature->start(), $feature->end(), $nameOverride);
}

sub featureToBedStringWithManualPositions($$$$$;$) {
    ## Generates a BED-compatible string from an Ensembl-compatible Perl object.
    ## See Ensembl docs at: http://uswest.ensembl.org/info/docs/api/core/core_tutorial.html?redirect=mirror;source=www.ensembl.org#conventions
    # $feature should be a gene or transcript or exon, usually.
    # $gene is the gene associated with this $feature. Note that $feature can ALSO be a gene!
    #    - note: $gene ALSO be just literal text, in which case it is printed out literally!
    # $tscript is the transcript associated with this $feature.
    #    - note: $tscript can ALSO be just literal text, in which case it is printed out literally!
    #    - example of this: a GENE doesn't have a specific transcript to report, so we just set $tscript = "NA";
    my ($feature, $gene, $tscript, $left, $right, $nameOverride) = @_;

    my $name;
    if (defined($nameOverride)) {
	$name = $nameOverride; ## If a name OVERRIDE was specified, then use that.
    } else {
	## Otherewise, let's just get the name in a convenient fashion!
	$name  = ($feature->can('stable_id')) ? $feature->stable_id() : "NA"; ## If this feature HAS a stable_id method, then call it! Otherwise (for introns, mainly) don't call it!
    }
    my $seq_region = $feature->slice->seq_region_name(); ## chromosome
    my $start      = defined($left) ? $left : "undef";
    my $end        = defined($right) ? $right : "undef";
    my $strand     = $feature->strand();
    my $score      = $DEFAULT_SCORE_FOR_THINGS_WITHOUT_A_SCORE;
    my $geneEnsemblID    = ($gene->can('stable_id')) ? $gene->stable_id() : $gene; ## if $gene ISN'T an object with a stable_id, we just hope it's a string and print it directly.
    my $tsID = ($tscript->can('stable_id')) ? $tscript->stable_id() : $tscript; ## if $tscript ISN'T an object with a stable_id, we just hope it's a string and print it directly.
    my $biotype;
    if ($feature->can('biotype')) { $biotype = $GLOBAL_BIOTYPE_PREFIX . $feature->biotype(); }
    elsif ($tscript->can('biotype')) { $biotype = $GLOBAL_BIOTYPE_PREFIX . $tscript->biotype(); }
    elsif ($gene->can('biotype')) { $biotype = $GLOBAL_BIOTYPE_PREFIX . $gene->biotype(); }
    else { $biotype = $GLOBAL_BIOTYPE_PREFIX . "NO_BIOTYPE"; }

    my $geneCommonName = ($gene->can('external_name')) ? $gene->external_name() : 'NA';
#    print $gene->display_xref->display_id() . " < GS\n";

    return manualToBedString($name, $seq_region, $start, $end, $strand, $score, $geneCommonName, $geneEnsemblID, $tsID, $biotype);
}


sub manualToBedString($$$$$$$$$$) {
    # Generate a BED-comptaible string from a list of parameters. Requires specifying the "Score" parameter, which is usually the dummy value '0'.
    my ($name, $chr, $start, $end, $strand, $score, $geneCommonName, $geneID, $tscript, $biotype) = @_;
    if ($strand eq '1' || $strand eq '+1') { $strand = '+'; }  ## '1' should be reported as '+'
    if ($strand eq '-1') { $strand = '-'; } ## '-1' should be reported as '-'
    return ($chr . "\t" . $start . "\t" . $end . "\t" . $name . "\t" . $score . "\t" . $strand . "\t" . $geneCommonName . "\t" . $geneID . "\t" . $tscript . "\t" . $biotype);
}

# ====================================================

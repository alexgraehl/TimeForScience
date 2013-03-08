#!/usr/bin/perl

### =============================================================
### =============================================================
## You need to include some links to all the ENSEMBL Perl libraries, PLUS the bioperl libraries. Note their location on the filesystem below:
use lib "/work/Apps/perl/ensembl_perl_api/bioperl-live";
use lib "/work/Apps/perl/ensembl_perl_api/ensembl/modules";
use lib "/work/Apps/perl/ensembl_perl_api/ensembl-compara/modules";
use lib "/work/Apps/perl/ensembl_perl_api/ensembl-variation/modules";
use lib "/work/Apps/perl/ensembl_perl_api/ensembl-functgenomics/modules";

use lib "/work/Apps/perl/ensembl_perl_api/biomart-perl/lib";
### =============================================================
### =============================================================

use strict;  use warnings;  use diagnostics;

# An example script demonstrating the use of BioMart API.
# This perl API representation is only available for configuration versions >=  0.5 
use BioMart::Initializer;
use BioMart::Query;
use BioMart::QueryRunner;

print STDERR "Alex: note, you will have to use CPAN to install 'Log::Log4perl'  before this will work!\n";
print STDERR "Alex: note, you will have to use CPAN to install 'Readonly'  before this will work!\n";
print STDERR "Alex: as root, try: perl -MCPAN -e shell, and then install those scripts from within the interactive CPAN shell. Make sure users (not just root) can READ the resulting installed perl library files, of course, too!\n";
print STDERR "Alex: I guess you have to download a file from: http://www.biomart.org/biomart/martservice?type=registry and save it to a local path. Then set it below! Not sure though... it definitely does not work to set it to a URL, though!\n";

my $confFile = "./alex-biomart-registry.txt"; #"PATH TO YOUR REGISTRY FILE UNDER biomart-perl/conf/. For Biomart Central Registry navigate to http://www.biomart.org/biomart/martservice?type=registry";

if (!-f $confFile || (-s $confFile == 0)) {
    print STDERR "Alex: Whoops, there was no conf file. Let's try to download one...";
    my $result = system(qq{curl "http://uswest.ensembl.org/biomart/martservice?type=registry" > $confFile}); 
    #http://www.biomart.org/biomart/martservice?type=registry" > $confFile});
    ($result == 0) or die "Curl command apparently failed!";
    print STDERR "I hope that worked!";
}

#
# NB: change action to 'clean' if you wish to start a fresh configuration  
# and to 'cached' if you want to skip configuration step on subsequent runs from the same registry
#

my $globalRegistry = undef;

sub globalInitRegistryOnlyOnceEver() {
    if (!defined($globalRegistry)) {
	print STDERR "Initializing the global registry...\n";
	my $action='cached';
	my $initializer = BioMart::Initializer->new('registryFile'=>$confFile, 'action'=>$action);
	$globalRegistry = $initializer->getRegistry;
    } else {
	print STDERR "The global registry has already been initialized.\n";
    }
}

#letsDownloadStuff({ filename      => 'mouse_ens_master'
#			, species => 'mmusculus'
#		  });

#letsDownloadStuff({ filename      => 'chicken_ens_master'
#			, species => 'ggallus'
#		  });

my @latinSpecies = ('hsapiens', 'drerio', 'ggallus', 'mmusculus');
foreach my $source (@latinSpecies) {
    letsDownloadStuff({ filename => "DL_BIOMART_PRIMARY_${source}"   , species => "$source" });
    # Example: letsDownloadStuff({ filename      => 'human_ens_master', species => 'hsapiens'});
    letsDownloadStuff({ filename => "DL_BIOMART_TRANSCRIPT_${source}", species => "$source", get_transcripts_only => 1}); # yes, we ONLY want to get the transcripts + gene mapping, not any of the other data!

    foreach my $target (@latinSpecies) {
	next if ($source eq $target); # no point in mapping a species to ITSELF!
	letsDownloadStuff({ filename               => "DL_BIOMART_${source}_map_to_${target}"
				, species          => $source
				, ortholog_species => $target
			  });
    }
}



sub letsDownloadStuff {
    my ($aa) = @_; ## aa is a HASH REFERENCE. You access it like so:  $aa->{'hashkey'}
    # hash should have these elements (potentially; some may be missing!):
    # filename: REQUIRED (string)
    # species: REQUIRED (string)
    # ortholog_species: OPTIONAL (string)
    # get_transcripts_only: OPTIONAL (1 / 0 for true/false)

    if (-f $aa->{'filename'} && (-s $aa->{'filename'} > 0)) {
	print STDERR "[OMITTING] the re-download of file <" . $aa->{'filename'} . ">, as that file was already present and had non-zero size. Delete it from the command line if you REALLY want to download it again.\n";
	return; # exit early!
    }
    print STDERR "[DOWNLOADING] the file <" . $aa->{'filename'} . ">, for species \"" . $aa->{'species'} . "\"... note that this is SLOW because it's a network operation.\n";
    globalInitRegistryOnlyOnceEver();
    my $query = BioMart::Query->new('registry'=>$globalRegistry,'virtualSchemaName'=>'default');
    $query->setDataset( $aa->{'species'} . "_gene_ensembl" );  ## <-- change datasets!
    $query->addAttribute("ensembl_gene_id");

    if (exists($aa->{'get_transcripts_only'}) && defined($aa->{'get_transcripts_only'}) && $aa->{'get_transcripts_only'}) {
	print STDERR "Running in special \"GET ONLY THE GENE IDs AND TRANSCRIPTS\" mode.\n";
	$query->addAttribute("ensembl_transcript_id");
    } else {

	$query->addAttribute("external_gene_id");
	$query->addAttribute("external_gene_db");

	if (exists($aa->{'ortholog_species'}) && $aa->{'ortholog_species'}) {
	    my $orthSpecies = $aa->{'ortholog_species'};
	    print STDERR "Getting ORTHOLOG MAPPING for the SINGLE SPECIES ONLY named $orthSpecies. Note that it is intentional to only be able to get one at a time, due to the way the results are generated. Multiple results seem to generate the 'cross product' of all possible mappings.\n";

	    $query->addAttribute("${orthSpecies}_homolog_ensembl_gene"); # orthology-specific attributes
	    $query->addAttribute("${orthSpecies}_homolog_orthology_type"); # note that this must be done BEFORE
	    $query->addAttribute("${orthSpecies}_homolog_subtype");        # we change datasets (below)!
	    $query->addAttribute("${orthSpecies}_homolog_perc_id");
	    $query->addAttribute("${orthSpecies}_homolog_perc_id_r1");
	    # Now SWITCH to that other database to get the full data...
	    #$query->setDataset("${orthSpecies}_gene_ensembl"); ## <-- change datasets!
	    #addStandardAttributes($query); # for the ORTHOLOG species
	} else {
	    print STDERR "Adding the standard attributes that are incompatible with the orthology mapping.\n";
	    $query->addAttribute("chromosome_name");
	    $query->addAttribute("start_position");
	    $query->addAttribute("end_position");
	    $query->addAttribute("strand");
	    $query->addAttribute("percentage_gc_content");
	    $query->addAttribute("transcript_count");
	    $query->addAttribute("status");
	    $query->addAttribute("refseq_mrna");
	    $query->addAttribute("refseq_ncrna");
	    $query->addAttribute("unigene");
	    $query->addAttribute("description");
	}
    }

    $query->formatter("TSV");
    my $query_runner = BioMart::QueryRunner->new();
############################## GET COUNT ############################
# $query->count(1);
# $query_runner->execute($query);
# print $query_runner->getCount();
#####################################################################

############################## GET RESULTS ##########################
# to obtain unique rows only
# $query_runner->uniqueRowsOnly(1);

    my $QUERYOUT;
    open $QUERYOUT, '>', $aa->{'filename'} or die "Can't open file for writing: $!";

    $query_runner->execute($query);
    $query_runner->printHeader($QUERYOUT);
    $query_runner->printResults($QUERYOUT);
    $query_runner->printFooter($QUERYOUT);
    
    close($QUERYOUT);
#####################################################################
}

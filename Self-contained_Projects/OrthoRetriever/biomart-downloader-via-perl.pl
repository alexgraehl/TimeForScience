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
print STDERR "Alex: I guess you have to download a file from: http://www.biomart.org/biomart/martservice?type=registry and save it to a local path. Then set it below! Not sure though... it definitely does not work to set it to a URL, though!";

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
    if (!defined($globalInitRegistryOnlyOnceEver)) {
	print STDERR "Initializing the global registry...\n";
	my $action='cached';
	my $initializer = BioMart::Initializer->new('registryFile'=>$confFile, 'action'=>$action);
	$globalRegistry = $initializer->getRegistry;
    } else {
	print STDERR "The global registry has already been initialized.\n";
    }
}


letsDownloadStuff({ filename      => 'human_ens_master'
			, species => 'hsapiens'
		  });

letsDownloadStuff({ filename      => 'mouse_ens_master'
			, species => 'mmusculus'
		  });

letsDownloadStuff({ filename      => 'chicken_ens_master'
			, species => 'ggallus'
		  });

letsDownloadStuff({ filename      => 'zebrafish_ens_master'
			, species => 'drerio'
		  });

letsDownloadStuff({ filename      => 'BIOMART_hsapiens_map_to_drerio'
			, species => 'hsapiens'
			, ortholog_species => 'drerio'
		  });



sub letsDownloadStuff {
    my ($aa) = @_; ## aa is a HASH REFERENCE. You access it like so:  $aa->{'hashkey'}

    if (-f $aa->{'filename'} && (-s $aa->{'filename'} > 0)) {
	print STDERR "[OMITTING] the re-download of file <" . $aa->{'filename'} . ">, as that file was already present and had non-zero size. Delete it from the command line if you REALLY want to download it again.\n";
	return; # exit early!
    }


    print STDERR "[DOWNLOADING] the file <" . $aa->{'filename'} . ">, for species \"" . $aa->{'species'} . "\"... note that this is SLOW because it's a network operation.\n";

    globalInitRegistryOnlyOnceEver();

    my $query = BioMart::Query->new('registry'=>$globalRegistry,'virtualSchemaName'=>'default');

    $query->setDataset( $aa->{'species'} . "_gene_ensembl" );  ## <-- change datasets!

    sub addStandardAttributes($) {
	my ($theQ) = @_;
	$theQ->addAttribute("ensembl_gene_id");
	$theQ->addAttribute("chromosome_name");
	$theQ->addAttribute("start_position");
	$theQ->addAttribute("end_position");
	$theQ->addAttribute("strand");
	$theQ->addAttribute("external_gene_id");
	$theQ->addAttribute("external_gene_db");
	$theQ->addAttribute("percentage_gc_content");
	$theQ->addAttribute("transcript_count");
	$theQ->addAttribute("status");
	$theQ->addAttribute("refseq_mrna");
	$theQ->addAttribute("refseq_ncrna");
	$theQ->addAttribute("unigene");
	$theQ->addAttribute("description");
    }
    
    addStandardAttributes($theQ); # for the DEFAULT species

    print STDERR "Getting ORTHOLOG MAPPING for the SINGLE SPECIES ONLY named $ortholog_species...\n";

    if (exists($aa->{'ortholog_species'}) && $aa->{'ortholog_species'}) {
	my $orthSpecies = $aa->{'ortholog_species'};

	$query->addAttribute("${orthSpecies}_homolog_ensembl_gene"); # orthology-specific attributes
	$query->addAttribute("${orthSpecies}_homolog_orthology_type"); # note that this must be done BEFORE
	$query->addAttribute("${orthSpecies}_homolog_subtype");        # we change datasets (below)!
	$query->addAttribute("${orthSpecies}_homolog_perc_id");
	$query->addAttribute("${orthSpecies}_homolog_perc_id_r1");

	# Now SWITCH to that other database to get the full data...
	$query->setDataset("${orthSpecies}_gene_ensembl"); ## <-- change datasets!
	addStandardAttributes($query); # for the ORTHOLOG species
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

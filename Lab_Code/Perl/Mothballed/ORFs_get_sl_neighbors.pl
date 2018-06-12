#!/usr/bin/env perl

use strict;
use warnings;

# This is sort of a weird script.
# It takes in THREE things, but only two of them are in argv (the other must be passed in by STDIN).

# In a Makefile, you run it like so:

#  sl_neighbors_of_the_thing:   list_of_sls_for_each_orf     drug_sensitivity_list
#       $(RM) $@ ;
#       cat drug_sensitivity_list
#           | perl ORFs_get_sl_neighbors.pl OUTPUT_FILE  list_of_sls_for_each_orf

# Where list_of_sls_for_each_orf is list of synthetic lethal neighbors like this:
#    SOME_ORF_1     SL_NEIGHBOR_ORF_1    SL_NEIGHBOR_ORF_3    SL_NEIGHBOR_ORF_4
#    SOME_ORF_2     SL_NEIGHBOR_ORF_17   SL_NEIGHBOR_ORF_12   SL_NEIGHBOR_ORF_41
#      ...   etc

# And drug_sensitivity_list (which goes in through STDIN) is like this:
#    drug_1     SENSITIVE_ORF_1    SENSITIVE_ORF_3    SENSITIVE_ORF_4
#    drug_2     SENSITIVE_ORF_2    SENSITIVE_ORF_3    SENSITIVE_ORF_6
#      ...   etc

my $outName  = $ARGV[0];
my $slNeighborhoodsFilename  = $ARGV[1];
my $filtering = $ARGV[2];

# filtering can be "regular_sl_neighbors"
# or "neighbors_plus_members"
# or "neighbors_minus_members"

#print "Got these arguments: $outName\n$slNeighborhoodsFilename\n";

use constant NEIGHBORS_PLAIN => 1;
use constant NEIGHBORS_PLUS_MEMBERS => 2;
use constant NEIGHBORS_MINUS_MEMBERS => 3;

my $tab = "\t";
my $newline = "\n";

my $dollarSign = '$';

my $modifyLethalList = '';
	if ($filtering eq "regular_sl_neighbors") {
		$modifyLethalList = NEIGHBORS_PLAIN;
	} elsif ($filtering eq "neighbors_plus_members") {
		$modifyLethalList = NEIGHBORS_PLUS_MEMBERS;
	} elsif ($filtering eq "neighbors_minus_members") {
		$modifyLethalList = NEIGHBORS_MINUS_MEMBERS;
	} else {
		die "ERROR: An improper filtering parameter ($filtering) was passed in.";
	}


while (my $drugSensLine = <STDIN>) {

	# Each $drugSensLine should look something like this:
	# SOME_DRUG_HERE   <tab>   ORF_1   <tab>   ORF_2   <tab>   ORF_3

	#print $drugSensLine . "\n";

	$drugSensLine =~ m/^([^\t]+)/;
	my $drugName = $1; #print $1 . "\n";
	chomp($drugName);

	my $command =
		qq{echo "$drugSensLine"} .
		qq{| cut -f 2- } .
		qq{| transpose.pl -q } .
		(($modifyLethalList == NEIGHBORS_MINUS_MEMBERS) ? qq{| tee drugSens.tmp }   : '' ) .  # save the who-is-in-this-drug-sensitive-list if we want to remove them later
		qq{| join.pl -ob - $slNeighborhoodsFilename } .
		(($modifyLethalList != NEIGHBORS_PLUS_MEMBERS)  ? qq{| cut -f 2- }          : '' ) .  # if we want to ADD the members of this pathway, then don't cut out the first column (which gives all the drug-sensitive ORFs)
		qq{| sed 's/\\t/\\n/g' } .
		(($modifyLethalList == NEIGHBORS_MINUS_MEMBERS) ? qq{| join.pl -neg - drugSens.tmp}   : '') .
		qq{| sort -u } .
		qq{| sed '/^\$/d' } .
		qq{| transpose.pl -q } .
		qq{| sed 's/^/$drugName\t/' } .
		qq{ >> $outName };

#	print "COMMAND\n\n" . $command . "\n";

#		qq{| sed "s/^$dollarSign/d" } .

	system($command);
}




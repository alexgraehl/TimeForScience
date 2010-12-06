include $(HOME)/Map/Papers/Cozen07/Makefile.common

type           = $(notdir $(shell pwd))
mountain_file  = ../../Clusters/Mountains/Large/$(type)/data.tab
neighbor_file  = ../../Network/$(type)/degree.tab
orthos_file    = ../Remote/Pae_all.tab

# Supply a Makefile.common file to overwrite variables.
ifeq (Makefile.common,$(shell exists.pl Makefile.common))
include Makefile.common
endif

# Overwrite these parameters here or set them in the the Makefile.common in the
# current directory.

# File containing a list of other files that should be used for annotations.
annotation_list_file ?= annotations.tab

# If the annotations file exists, then get the
# lists of files out of it so we can set up the
# dependencies correctly.
ifneq (,$(shell exists.pl $(annotation_list_file)))
annotation_files = $(shell cat $(annotation_list_file) | cut -f 1 | grep -v '^\#')
else
annotation_files =
endif

# A file describing a network of interacting genes must be supplied. This is file
# should be stored in the variable $(link_file). This file should be a tab-
# delimited file where two columns correspond to gene identifiers describing
# pairs of interacting genes and an optional third column describes the strength
# of each interaction.
#
# The name of the file containing the network of interacting genes; default is links.tab
link_file ?= $(root)/Network/$(THISDIR)/links.tab
# The column gene 1 is in in the links file; default is 1.
col_link_gene ?= 1
# The column the gene's neighbor is in within the links file.
col_link_nbr ?= 2
# The column the link strength is in. If this is 0 then it is assumed
# that no strengths are available in which case unitary weights for each
# link are used instead. If this value is -1 it is assumed that a similarity
# file is already available in which case a link to the $(link_file) is
# created instead of making a new file.
col_link_weight ?= 3

# The similarity file is a tab-delimited file with GENE1, GENE2, and a similarity
# strength on each line. This file is created from the links file. Commented
# out the ability to rename the variable to reduce confusion.
# sim_file ?= data.sim

# The coordinates files is created by running VxOrd on the similarity file that was
# created from the links file. By default this file is named data.coord. The data.coord
# file will be tab-delimited with the first column equal to a gene identifier and the
# second and third columns equal to the X and Y coordinate of the gene in the topomap
# respectively. Commented out the ability to rename the variable to reduce confusion.
# coord_file ?= data.coord

# The description file should be tab-delimited and have, on each line:
# GENE \t NAME \t DESCRIPTION [\t ATTRIB1 \t ATTRIB2 ...] where GENE is the
# identifier of the gene and should match the identifiers used in the links file;
# NAME is the common name or symbol used for the gene; and DESCRIPTION is a short
# description of the function of the gene if known. The data can have optional
# additional columns describing attributes of the gene. The first line of the
# file is assumed to be a header line that gives the names of each column.
# If the flag $(desc_has_header) is set to 0, then the file is forced to have
# three columns of data (everything after the second field is treated as the
# description column) and a header is automatically generated.
desc_file       ?= desc.tab
desc_has_header ?= 0
# The column holding the gene identifiers.
col_desc_gene    ?= 1
# The column containing the common gene name or symbols.
col_desc_name    ?= 2
# The column containing a description of the genes functions.
col_desc_desc    ?= 3
# The columns containing additional attributes about the genes. Default is blank.
# A range or list can be specified such as 1,2-5,8
col_desc_attribs ?=

# The Microsoft Access file contains tables with gene annotations.
# access_file = data.mdb

# The connectivity file tells VxInsight how to find information in the
# database.
# connect_file = data.db

# The config file tells VxInsight how to display the information in the GUI.
# config_file = data.config

# The neighboring links to have in the nbrs.tab file that will be loaded. Typical
# values are 10, 20, or 50.
num_nbrs = 10

targets  = data.tab data.sim desc.tab \
	   data.mdb data.db data.config nbrs.tab

include $(MAPDIR)/Templates/Make/quick.mak

data.config: data.tab
	echo 'puts "Reading config file..."' \
	> $@; \
	head -n 1 $< \
	| transpose.pl -q \
	| grep -v '^X$$' \
	| grep -v '^Y$$' \
	| paste.pl 'Add_Small_Info_Field' - '1 10 green' \
	>> $@; \
	echo 'Add_Connection_Field	Links' \
	>> $@; \

data.db: data.tab
	head -n 1 $< \
	| cut -f 1 \
	| sed 's/^/Id\tautoconnect*data.mdb::Data::/' \
	| sed 's/$$/\tcache/' \
	> $@; \
	head -n 1 $< \
	| transpose.pl -q \
	| paste.pl - 'autoconnect*data.mdb::Data::' \
	| cut.pl -f 1,2,1 \
	| delim2delim.pl -f 2 '	' '' \
	| paste.pl - 'cache' \
	>> $@; \
	echo 'X	autoconnect*data.mdb::Data::K_X	cache' \
	>> $@; \
	echo 'Y	autoconnect*data.mdb::Data::L_Y	cache' \
	>> $@; \
	echo 'Links	autoconnect*data.mdb::Nbrs::B_Neighbor	connection' \
	>> $@; \

# Currently, we have to load the data.tab and the data.sim files into a Microsoft
# Access Database for VxInsight to run and be able to pull up annotations.
data.mdb: data.tab nbrs.tab
	echo 'Load data.tab and nbrs.tab into Microsoft Access and save as $@'; \
	echo 'Name the table 'Data' for the data.tab file.'; \
	echo 'Name the table 'Nbrs' for the nbrs.tab file.'; \

data.tab: $(annotation_list_file) data.coord
	cp $(mountain_file) ../mountains.tmp; \
	cat $(neighbor_file) \
	| cut -f 1,3- \
	| sed 's/	/ /g' \
	| sed 's/ /	/' \
	> ../neighbors.tmp; \
	cat $(orthos_file) \
	| cut -f 1,13,14,15,16 \
	| sed 's/	/ /g' \
	| sed 's/ /	/' \
	> ../orthos.tmp; \
	join_batch.pl -ob -b NaN -dos $(annotation_list_file) \
	> 1.tmp; \
	head -n 1 1.tmp \
	> $@; \
	cat data.coord \
	| dos2unix \
	| join.pl - 1.tmp \
	>> $@; \

nbrs.tab: data.sim data.tab
	head -n 1 data.tab \
	| cut -f 1 \
	| paste.pl - 'B_Neighbor' \
	> $@; \
	cut.pl -f 2,1,3 $< \
	| cat - $< \
	| sort -k 3,3 -t '	' -g -r \
	| cut -f 1-2 \
	| expand.pl \
	> 1.tmp; \
	cut -f 1 1.tmp \
	> 2.tmp; \
	cut -f 2- 1.tmp \
	| cut -f 1-$(num_nbrs) \
	| paste 2.tmp - \
	| flatten.pl \
	>> $@; \
	rm -f 1.tmp 2.tmp; \

# Currently, we have to run VxOrd on a Windows computer. Hopefully someone
# will implement the VxOrd algorithm so we can call it in the pipeline!
data.coord: data.sim
	echo 'Run VxOrd on $< and save the output as $@'

# Create the similarity file from the provided links file.
# If the similarity file exists with GENE1 \t GENE2 \t WEIGHT
# on each line, just link the data.sim file to this pre-existing
# file.
ifeq (-1,$(col_link_weight))

# Make sure the $(link_file) does not already equal the data.sim file!
ifneq (data.sim,$(link_file))
data.sim: $(link_file)
	rm -f $@; \
	ln -s $< $@;
endif

else

# If no weight column exists, paste in a 1 after each gene pair on each
# line.
ifeq (0,$(col_link_weight))
data.sim: $(link_file)
	cat.pl $< \
	| dos2unix \
	| grep -v '^#' \
	| cut.pl -f $(col_link_gene),$(col_link_nbr) \
	| paste.pl - '1' \
	> $@;

# A weight column exists, so just cut it out of the links file.
else
data.sim: $(link_file)
	cat.pl $< \
	| dos2unix \
	| grep -v '^#' \
	| cut.pl -f $(col_link_gene),$(col_link_nbr),$(col_link_weight) \
	> $@;
endif

endif

ifneq (desc.tab,$(desc_file))

ifeq (1,$(desc_has_header))
desc.tab: $(desc_file)
	cat.pl $< \
	| dos2unix \
	| cut.pl -f $(col_desc_gene),$(col_desc_name),$(col_desc_desc),$(col_desc_attribs) \
	> $@;
else
desc.tab: $(desc_file)
	cat.pl $< \
	| dos2unix \
	| cut.pl -f $(col_desc_gene),$(col_desc_name),$(col_desc_desc),$(col_desc_attribs) \
	| delim2delim.pl -f 3- '	' ' ' \
	| cap.pl Id,Gene,Description \
	> $@;
endif

endif

$(annotation_list_file):
	echo 'Make a file that lists each data file you want to supply for annotations.'; \
	echo 'Should conform to the format needed by join_batch.pl:'; \
	join_batch.pl --help; \



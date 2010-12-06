include $(HOME)/Map/Makefile.common
include $(HOME)/Map/Data/Makefile.common
include $(HOME)/Map/Papers/ExpressionNetwork/Makefile.common

ifneq (,$(shell exists.pl Makefile.common))
   include Makefile.common
endif

# For this makefile to work, you must make a
# Makefile.common file in the same directory with
# the variable definition:
#
# pub_id = PUBID
#
# where PUBID is equal to an SMD publication identifier.

ifeq (,$(gene_field))
   gene_field = $(SMD_$(ORGANISM_UPPER)_GENE_FIELD)
endif

ifeq (,$(gene_field))
   gene_field = NAME
endif

ifeq (,$(alias))
   alias = $(SMD_$(ORGANISM_UPPER)_ALIAS)
endif

ifeq (,$(alias))
   alias_file  = $(MAP_DATA)/Aliases/Gene/$(ORGANISM)/$(alias)/data.tab
   map_aliases = join.pl -1 2 $(alias_file) - \
                 | cut -f 2-
else
   map_aliases = cat
endif

pub_meta    = Remote/publication_$(pub_id).meta
smd_dir     = Remote/share
targets     = data.tab

all: $(targets)

clean:
	$(RM) $(wildcard *.tmp) $(targets)

make:

remote:
	download_smd_pub_meta.pl $(pub_id); \

data.tab: $(pub_meta)
	download_smd_pub_data.pl $(pub_meta); \
	join_smd.pl -g '$(gene_field)' Remote \
	| $(map_aliases) \
	> data.tab; \
	$(RMDIR) $(smd_dir); \
	$(RM) $(wildcard Remote/*.xls); \


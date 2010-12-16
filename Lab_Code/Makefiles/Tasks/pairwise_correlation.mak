include $(HOME)/Map/Makefile.common
include $(HOME)/Map/Data/Makefile.common

include $(HOME)/Map/Templates/Make/run.mak

DATASET           = CompendiumSlim
MEG_FILES         = $(foreach org, $(ORGANISMS), $(org).tab)
EXAMPLE_MEG_FILE  = $(word 1, $(MEG_FILES))

MEG_LISTS         = $(foreach org, $(ORGANISMS), $(org).lst)
EXAMPLE_MEG_LIST  = $(word 1, $(MEG_LISTS))

all: data.tab

clean:
	$(RM) data.tab

data.tab: data.map
	$(MAP_ERAN) data.map; \
	cat data.tmp \
	| grep -v '\-1.79769e+308' \
	| grep -v 'NaN' \
	| select.pl -gte '$(MIN_DIMENSIONS)' -k 4 \
	> $(TEMP_1); \
	mv $(TEMP_1) data.tmp; \
	cut -f 1-3 data.tmp \
	> $(TEMP_1); \
	cut -f 1 $(TEMP_1) \
	| join.pl -o 'NO_MEG' - $(MEG_DIR)/$(ORGANISM)/data.meg \
	| cut -f 2 \
	> $(TEMP_2); \
	cut -f 2 $(TEMP_1) \
	| join.pl -o 'NO_MEG' - $(MEG_DIR)/$(ORGANISM)/data.meg \
	| cut -f 2 \
	> $(TEMP_3); \
	cut -f 3 $(TEMP_1) \
	| paste $(TEMP_2) $(TEMP_3) - \
	| grep -v 'NO_MEG' \
	| order_keys.pl \
	| uniq.pl -f 1,2 -s 3,max \
	> data.tab; \

# data.map: template.map
# 	bind.pl template.map org=$(ORGANISM) dataset=$(DATASET) \
# 	> data.map

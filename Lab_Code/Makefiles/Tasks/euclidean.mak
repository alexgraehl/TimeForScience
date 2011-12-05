include $(HOME)/Map/Makefile.common
include $(HOME)/Map/Data/Makefile.common

STDEV            = 1.96
MIN_DIMENSIONS   = 50

COMPENDIUMS       = $(foreach org, $(ORGANISMS), \
		    $(MAP_DATA)/Expression/$(org)/Compendium/data.tab )

all: $(org1).tab

make:
	$(RM) data.map; \
	$(LINK) $(MAP_TEMPLATES)/Runs/euclidean.map data.map; \

include $(HOME)/Map/Templates/Make/run.mak

$(org1).tab: $(org1)_genes.tab
	date > $(org1).tab.time; \
	$(foreach org, $(ORGANISMS), \
	   cut -f 1 $(org)_genes.tab \
	   | join.pl -o 'NO_MEG' - $(MEG_DIR)/$(org)/data.meg \
	   | cut -f 2 \
	   > $(org)_1.tmp; \
	   cut -f 2 $(org)_genes.tab \
	   | join.pl -o 'NO_MEG' - $(MEG_DIR)/$(org)/data.meg \
	   | cut -f 2 \
	   > $(org)_2.tmp; \
	   cut -f 3,4 $(org)_genes.tab \
	   | paste $(org)_1.tmp $(org)_2.tmp - \
	   | grep -v 'NO_MEG' \
	   | order_keys.pl \
	   | sort -t '	' -k 1,2 \
	   | uniq.pl -o -f 1,2 -s 3,max \
	   | select.exe -k 3 -gte -0.3 \
	   > $(org).tab; \
	   $(RM) $(org)_1.tmp $(org)_2.tmp; \
	) \
	date >> $(org1).tab.time; \

$(org1)_genes.tab: $(COMPENDIUMS) data.map
	date > $(org1)_genes.tab.time; \
	bind.pl d.map organisms=$(ORGANISMS_COMMA) -exe $(MAP_EXE); \
	$(foreach org, $(ORGANISMS), \
	   cat $(org)_genes.tab \
	   | grep -v '\-1.79769e+308' \
	   | grep -v 'NaN' \
	   | select.exe -k 4 -gte '$(MIN_DIMENSIONS)' \
	   > $(org)_1.tmp; \
	   mv $(org)_1.tmp $(org)_genes.tab; )
	date >> $(org1)_genes.tab.time; \

data.tab:
	cut -f 1,2 $(foreach org, $(orgs), $(org).tab ) \
	| sort -t '	' -k 1,2 -u \
	$(foreach org, $(orgs), \
	   | join.pl -1 1,2 -2 1,2 -o 'NaN	NaN' - $(org).tab \
	   | cut.pl -f 1--2 ) \
	> data.tab


split_combo = $(shell echo '$(combo)' | sed 's/_/ /g')

$(org1)_pval.tab: $(org1)_genes_pval.tab
	$(foreach org, $(ORGANISMS), \
	   cut -f 1 $(org)_genes_pval.tab \
	   | join.pl - $(MEG_DIR)/$(org)/data.meg \
	   > left.tmp; \
	   cut -f 2 $(org)_genes_pval.tab \
	   | join.pl - $(MEG_DIR)/$(org)/data.meg \
	   > right.tmp; \
	   cut -f 3- $(org)_genes_pval.tab \
	   | paste left.tmp right.tmp - \
	   | cut -f 2,4,7 \
	   | order_keys.pl \
	   | grep -v NaN \
	   > $(org)_pval.tab; ) \
	   $(RM) *.tmp;

$(org1)_genes_pval.tab: $(org1)_genes.tab
	$(foreach org, $(ORGANISMS), \
	   pearson2pvalue.pl -k 3 -dim 4 $(org)_genes.tab \
	   | paste $(org)_genes.tab - \
	   > $(org)_genes_pval.tab; )

cutoffs.stat: $(STAT_FILES)
	cat *.stat \
	| awk '{print $$2 + $(STDEV)*$$3;}' \
	> tmp_cutoffs; \
	echo '$(ORGANISMS)' \
	| space2tab.pl \
	| transpose.pl \
	| paste - tmp_cutoffs \
	> cutoffs.stat; \
	rm -f tmp_cutoffs;

## Used


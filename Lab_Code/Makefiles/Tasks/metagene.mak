include $(HOME)/Map/Makefile.common
include $(HOME)/Map/Data/Makefile.common

ORGANISMS             = $(subst _,$(space),$(notdir $(PWD)))
blast_dir             = $(MAP_DATA)/Blast/Protein/Result
meg1                  = $(word 1, $(ORGANISMS))/data.meg
info1                 = $(word 1, $(ORGANISMS))/meg_info.tab
exe                   = $(MAP_RELEASE)

all:
	make blast.tab; \
	make $(meg1); \
	make $(info1); \
	make megs.tab; \
	$(RM) blast.tab; \
	make data.tab; \

clean:
	$(RMDIR) $(ORGANISMS) data.tab megs.tab; \

ex: $(org1)_ex.lst

multiplicity.tab:
	cut -f 1,2 data.tab \
	| grep -v 'No Fly gene' \
	| space2tab.pl \
	| row_stats.pl -count - \
	| cut -f 2 \
	| body.pl 2 -1 \
	| sort \
	| uniq -c \
	| sed 's/^[ ][ ]*//' \
	| space2tab.pl \
	| cut.pl -f 2,1 \
	> f; \
	cut -f 1,3 data.tab \
	| grep -v 'No Human gene' \
	| space2tab.pl \
	| row_stats.pl -count - \
	| cut -f 2 \
	| body.pl 2 -1 \
	| sort \
	| uniq -c \
	| sed 's/^[ ][ ]*//' \
	| space2tab.pl \
	| cut.pl -f 2,1 \
	> h; \
	cut -f 1,4 data.tab \
	| grep -v 'No Worm gene' \
	| space2tab.pl \
	| row_stats.pl -count - \
	| cut -f 2 \
	| body.pl 2 -1 \
	| sort \
	| uniq -c \
	| sed 's/^[ ][ ]*//' \
	| space2tab.pl \
	| cut.pl -f 2,1 \
	> w; \
	cut -f 1,5 data.tab \
	| grep -v 'No Yeast gene' \
	| space2tab.pl \
	| row_stats.pl -count - \
	| cut -f 2 \
	| body.pl 2 -1 \
	| sort \
	| uniq -c \
	| sed 's/^[ ][ ]*//' \
	| space2tab.pl \
	| cut.pl -f 2,1 \
	> y; \
	cut -f 1 f h w y \
	| sort \
	| uniq \
	| join.pl -o NaN - f \
	| join.pl -o NaN - h \
	| join.pl -o NaN - w \
	| join.pl -o NaN - y \
	> multiplicity.tab; \

data.tab: $(foreach org, $(ORGS), $(org)/data.meg)
	cut -f 2 $(foreach org, $(ORGS), $(org)/data.meg) \
	| sort -u \
	> $(tmp); \
	$(foreach org, $(ORGANISMS), \
	   join.pl -1 2 $(org)/data.meg $(tmp) \
	   | expand.pl \
	   | delim2delim.pl -f 2- '	' ' ' \
	   > $(org)_$(tmp); \
	) \
	cat $(tmp) \
	$(foreach org, $(ORGANISMS), \
	   | join.pl -o 'No $(org) gene' - $(org)_$(tmp) \
	) \
	| cap.pl MetaGene,$(ORGANISMS_COMMA) \
	> data.tab; \
	$(RM) $(foreach org, $(ORGS), $(org)_$(tmp)) $(tmp); \


$(info1):
	$(foreach org, $(ORGANISMS), \
	   cd $(org); \
	   cut -f 1 data.meg \
	   | join.pl -o 'No Name	No Description' - $(GENE_DIR)/Description/$(org)/data.tab \
	   | paste data.meg - \
	   | cut.pl -f 2,4- \
	   | sort -t '	' -k 1,1 -u \
	   | cap.pl MetaGene,Name,Description \
	   > meg_info.tab; \
	   cd ..; \
	) \

$(meg1): $(blast_dir)/$(ORG1)/$(ORG2)/data.tab
	make blast.tab; \
	bind.pl -exe $(exe) $(HOME)/Map/Templates/MetaGene/reciprocal.map \
	   blast_table=blast.tab; \
	\
	$(foreach org, $(ORGANISMS), \
	   $(MKDIR) $(org); \
	   cat $(org).meg \
	   | sort -t '	' -k 1,1 \
	   > $(org)/data.meg; \
	   $(RM) $(org).meg; \
	   cut -f 1 $(org)/data.meg \
	   > $(org)/data.lst; \
	) \

blast.tab: $(blast_dir)/$(ORG1)/$(ORG2)/data.tab
	cat $(shell echo '$(foreach org1, $(ORGANISMS), \
	                     $(foreach org2, $(ORGANISMS), \
	                        $(blast_dir)/$(org1)/$(org2)/data.tab \
	                        $(blast_dir)/$(org2)/$(org1)/data.tab \
	                     ) \
	                  )' \
	            | space2tab.pl \
		    | transpose.pl -q \
		    | $(FILTER_SAME_SUBDIRS) \
	   ) \
	> blast.tab; \

genes.tab: $(foreach org, $(ORGANISMS), $(org)/data.meg)
	echo 'Organism	Gene	MetaGene' \
	> genes.tab; \
	$(foreach org, $(ORGANISMS), \
	   concat.pl -p '$(org)	' -sf $(org)/data.meg \
	   >> genes.tab; \
	) \

megs.tab: $(foreach org, $(ORGANISMS), $(org)/data.meg)
	cut -f 2 $(foreach org, $(ORGANISMS), $(org)/data.meg) \
	| sort -u \
	| lin.pl -a -i 0 \
	| cap.pl MetaGene,Dummy \
	> megs.tab; \

$(org1)_ex.lst: megs.tab $(foreach org, $(ORGANISMS), $(org)/data.meg)
	$(foreach org, $(ORGS), \
	   join.pl -2 2 -neg megs.tab $(org)/data.meg \
	   | cut -f 1 \
	   | body.pl 2 -1 \
	   | sort -u \
	   > $(org)_ex.lst; \
	) \



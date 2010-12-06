include $(HOME)/Map/Makefile.common

DEGREE                 = 10

# The number of standard deviations that gives 1 in a million chance
# of getting a significant interaction (from a two-tailed gaussian).
STDEV                  = 4.8916

RANDOM_DIR             = ../Permuted
REAL_DIR               = ../Real

COMBINATIONS           = $(shell all_combinations.pl -d _ Fly Human Worm Yeast)

MEG_FILE         = $(MAP_DATA)/MetaGene/Reciprocal/$(ORGANISMS_MERGED)/$(org)/data.meg

STAT_FILES             = $(foreach org, $(ORGANISMS), $(RANDOM_DIR)/$(org).stat)
EXAMPLE_STAT_FILE      = $(word 1, $(STAT_FILES))

CORRELATION_FILES      = $(foreach org, $(ORGANISMS), $(REAL_DIR)/$(org).tab)
EXAMPLE_CORR_FILE      = $(word 1, $(CORRELATION_FILES))

SIG_GENE_FILES              = $(foreach org, $(ORGANISMS), $(org)_genes.tab)
EXAMPLE_SIG_GENE_FILE       = $(word 1, $(SIG_GENE_FILES))

SIG_MEG_FILES         = $(foreach org, $(ORGANISMS), $(org).tab)
EXAMPLE_SIG_MEG_FILE  = $(word 1, $(SIG_MEG_FILES))

SIG_MEG_LISTS         = $(foreach combo, $(COMBINATIONS), $(combo).lst)
EXAMPLE_SIG_MEG_LIST  = $(word 1, $(SIG_MEG_LISTS))

MAKE_SIG_FILES    = $(foreach org, $(ORGANISMS), \
		       cat $(RANDOM_DIR)/$(org).stat \
		       | awk '{ print $$2 + $(STDEV)*$$3; }' \
		       | cat - $(REAL_DIR)/$(org).tab \
		       | threshold.pl - -k 3 \
		       > $(org)_genes.tab; )

MAKE_SIG_MEG_FILES = $(foreach org, $(ORGANISMS), \
		         join.pl $(org)_genes.tab $(MEG_FILE) \
		         | join.pl - $(MEG_FILE) -1 2 \
		         | cut.pl -f 4,5,3 \
		         | order_keys.pl \
			 > tmp1; \
			 cat tmp1 \
			 | cut -f 1,2 \
			 | tab2space.pl \
			 > tmp2; \
			 cat tmp1 \
			 | cut -f 3 \
			 | paste tmp2 - \
		         > $(org).tab; \
			 rm -f tmp1 tmp2; )

all:
	make sig

make:
	$(RM) data.map; \
	$(LINK) $(MAP_TEMPLATES)/Run/pearson.map data.map; \

clean:
	rm -f *.tab; \

$(word 1, $(ORGANISMS))_triangles.tab: $(word 1, $(ORGANISMS)).tab
	$(foreach org, $(ORGANISMS), \
	   cat $(org).tab \
	   | space2tab.pl \
	   | triangles.pl \
	   > $(org)_triangles.tab; )

fhw_triangles.ps: fhw_triangles.dot
	dot -Tps fhw_triangles.dot \
	> fhw_triangles.ps

fhw_triangles.dot: fhw_triangles.tab
	cat fhw_triangles.tab \
	| tab2dotty.pl -desc $(MEG_DIR)/Info/Description/data.tab \
	> fhw_triangles.dot

fhw_triangles.tab: Fly_Human_Worm.tab
	cat Fly_Human_Worm.tab \
	| space2tab.pl \
	| cut -f 1,2 \
	| triangles.pl \
	> fhw_triangles.tab

fwy_triangles.ps: fwy_triangles.dot
	dot -Tps fwy_triangles.dot \
	> fwy_triangles.ps

fwy_triangles.dot: fwy_triangles.tab
	cat fwy_triangles.tab \
	| tab2dotty.pl -desc $(MEG_DIR)/Info/Description/data.tab \
	> fwy_triangles.dot

fwy_triangles.tab: Fly_Worm_Yeast.tab
	cat Fly_Worm_Yeast.tab \
	| space2tab.pl \
	| cut -f 1,2 \
	| triangles.pl \
	> fwy_triangles.tab

sig: $(EXAMPLE_SIG_GENE_FILE)
	make sig_meg

sig_meg: $(EXAMPLE_SIG_MEG_FILE) $(EXAMPLE_SIG_MEG_LIST)
	make overlaps.tab

split_combo = $(shell echo '$(combo)' | sed 's/_/ /g')

meg_overlaps.tab:
	rm -f meg_overlaps.tab; \
	\
	$(foreach combo, $(COMBINATIONS), \
	   cat $(combo).lst \
	   | $(NUM_LINES) \
	   | cap.pl '$(combo)' \
	   | transpose.pl \
	   >> meg_overlaps.tab; )

###############################################################################

$(EXAMPLE_SIG_GENE_FILE): $(STAT_FILES)
	$(MAKE_SIG_FILES)

$(EXAMPLE_SIG_MEG_FILE): $(SIG_GENE_FILES)
	$(MAKE_SIG_MEG_FILES)

$(EXAMPLE_SIG_MEG_LIST):
	$(foreach combo, $(COMBINATIONS), \
	   cat $(combo).tab \
	   | space2tab.pl \
	   | cut -f 1 \
	   > meg1.tmp; \
	   cat $(combo).tab \
	   | space2tab.pl \
	   | cut -f 2 \
	   | cat meg1.tmp - \
	   | sort -k 1,1 -u \
	   > $(combo).lst; \
	   rm -f meg1.tmp; )

###############################################################################

sig_overlaps: sig_combinations.lst
	join_combinations.pl sig_combinations.lst

overlaps: orig_combinations.lst
	join_combinations.pl orig_combinations.lst

sig_combinations.lst: $(SIG_MEG_FILES)
	/bin/ls $(SIG_MEG_FILES) \
	| all_combinations.pl - \
	> tmp_all_combinations; \
	cat tmp_all_combinations \
	| sed 's/.tab//g' \
	| sed 's/ /_/g' \
	| sed 's/$$/.tab/' \
	| paste - tmp_all_combinations \
	> sig_combinations.lst; \
	rm -f tmp_all_combinations; \


overlaps.tab: $(SIG_MEG_FILES)
	overlap_combinations.pl -k 1,2 -di '\s+' $(SIG_MEG_FILES) \
	| cut -f 1,2 \
	| sed 's/.tab//g' \
	| sed 's/,/_/g' \
	| body.pl 2 -1 \
	> overlaps.tab

multi_cellular.tab: Fly_Human_Worm.tab Fly_Human_Worm_Yeast.tab
	join.pl -neg Fly_Human_Worm.tab Fly_Human_Worm_Yeast.tab \
	> multi_cellular.tab

multi_cellular.lst: multi_cellular.tab
	cat multi_cellular.tab \
	| space2tab.pl \
	| cut -f 1 \
	> meg1.tmp; \
	cat multi_cellular.tab \
	| space2tab.pl \
	| cut -f 2 \
	| cat meg1.tmp - \
	| sort -u \
	> multi_cellular.lst; \
	rm -f meg1.tmp; \

multi_cellular.dot: multi_cellular.tab
	cat multi_cellular.tab \
	| space2tab.pl \
	| tab2dotty.pl -desc $(MEG_DESC) \
	> mutli_cellular.dot \

multi_cellular_desc.tab: multi_cellular.lst
	cat multi_cellular.lst \
	| join.pl -o 'NO DESCRIPTION' - $(MEG_DESC) \
	> multi_cellular_desc.tab; \

fwy$(DEGREE).pdf: fwy$(DEGREE).ps
	ps2pdf fwy$(DEGREE).ps fwy$(DEGREE).pdf \

fwy$(DEGREE).ps: fwy$(DEGREE).dot
	dot -Tps fwy$(DEGREE).dot \
	> fwy$(DEGREE).ps

fwy$(DEGREE).dot: fwy$(DEGREE).tab
	cat fwy$(DEGREE).tab \
	| tab2dotty.pl -desc $(MEG_DIR)/Info/Description/data.tab \
	> fwy$(DEGREE).dot

fwy$(DEGREE).tab:
	cat Fly_Worm_Yeast.tab \
	| space2tab.pl \
	| cut -f 1 \
	> 1.tmp; \
	\
	cat Fly_Worm_Yeast.tab \
	| space2tab.pl \
	| cut -f 2 \
	> 2.tmp; \
	\
	cat 1.tmp 2.tmp \
	| sort \
	| uniq -c \
	| sed 's/^[ ][ ]*//' \
	| select.pl -op '>=' $(DEGREE) -k 1 \
	| cut -f 2 \
	> fwy10.lst; \
	\
	cat Fly_Worm_Yeast.tab \
	| space2tab.pl \
	| cut.pl -f 1 \
	| lin.pl -a \
	| join.pl fwy10.lst - \
	| cut -f 2 \
	> 1.tmp; \
	\
	cat Fly_Worm_Yeast.tab \
	| space2tab.pl \
	| cut.pl -f 2 \
	| lin.pl -a \
	| join.pl fwy10.lst - \
	| cut -f 2 \
	> 2.tmp; \
	\
	cat 1.tmp 2.tmp \
	| sort -u \
	> 12.tmp; \
	\
	cat Fly_Worm_Yeast.tab \
	| space2tab.pl \
	| lin.pl \
	| join.pl 12.tmp - \
	| cut.pl -f 2,3 \
	| sort -k 1,1 -k 2,2 -t '	' \
	> fwy10.tab; \
	\






include $(HOME)/Map/Makefile.common
include $(HOME)/Map/Data/Makefile.common

include $(HOME)/Map/Templates/Make/run.mak

# The number of standard deviations that gives 1 in a million chance
# of getting a significant interaction (from a two-tailed gaussian).
STDEV = 4.8916

COMPENDIUMS       = $(foreach org, $(ORGANISMS), \
		    $(MAP_DATA)/Expression/$(org)/Compendium/permuted.tab )

METAGENES_FILE    = $(MAP_DATA)/MetaGene/Reciprocal/$(ORGANISMS_MERGED)/$(org)/data.meg

STAT_FILES        = $(foreach org, $(ORGANISMS), $(org).stat)
EXAMPLE_STAT_FILE = $(word 1, $(STAT_FILES))
CORRELATION_FILES = $(foreach org, $(ORGANISMS), $(org).tab)
EXAMPLE_CORR_FILE = $(word 1, $(CORRELATION_FILES))
COMPUTE_STATS     = $(foreach org, $(ORGANISMS), \
		    cat $(org).tab | \
		    stats.pl -k 3 > \
		    $(org).stat; )

SIG_FILES         = $(foreach org, $(ORGANISMS), $(org)_sig.tab)
EXAMPLE_SIG_FILE  = $(word 1, $(SIG_FILES))

SIG_META_FILES         = $(foreach org, $(ORGANISMS), $(org)_sig_meta.tab)
EXAMPLE_SIG_META_FILE  = $(word 1, $(SIG_META_FILES))

MAKE_SIG_FILES    = $(foreach org, $(ORGANISMS), \
		       cat $(org).stat \
		       | awk '{ print $$2 + $(STDEV)*$$3; }' \
		       | cat - $(org).tab \
		       | threshold.pl - -k 3 \
		       > $(org)_sig.tab; )

MAKE_SIG_META_FILES = $(foreach org, $(ORGANISMS), \
		         join.pl $(org)_sig.tab $(METAGENES_FILE) \
		         | join.pl - $(METAGENES_FILE) -1 2 \
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
		         > $(org)_sig_meta.tab; \
			 rm -f tmp1 tmp2; )

all: data.time
	make stat

make:

stat: $(EXAMPLE_STAT_FILE)
	make cutoffs.stat; \
	make sig

sig: $(EXAMPLE_SIG_FILE)
	make sig_meta

sig_meta: $(EXAMPLE_SIG_META_FILE)
	rm -rf overlaps; \
	make overlaps

overlaps: $(SIG_META_FILES) overlaps.mak
	mkdir -p overlaps; \
	rm -f overlaps/Makefile; \
	cd overlaps; \
	ln -s ../overlaps.mak Makefile; \
	make

$(EXAMPLE_SIG_META_FILE): $(SIG_FILES)
	$(MAKE_SIG_META_FILES)

$(EXAMPLE_STAT_FILE): $(EXAMPLE_CORR_FILE)
	$(COMPUTE_STATS) \

$(EXAMPLE_CORR_FILE): $(COMPENDIUMS)
	$(MAP_EXE) data.map

$(EXAMPLE_SIG_FILE): $(STAT_FILES)
	$(MAKE_SIG_FILES)

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



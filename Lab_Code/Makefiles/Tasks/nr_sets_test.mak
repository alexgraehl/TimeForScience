
thresholds     = $(shell nums.pl 10 | transpose.pl -q | lin.pl -0 -i 10 | cut -f 1 | sed 's/^\([0-9]\)$$/0\1/')
CHILDREN       = $(foreach t,$(thresholds),NR$(t))
CHILD_MAKEFILE = $(MAPDIR)/Templates/Make/nr_sets.mak
lacks          = $(foreach c,$(CHILDREN),$(c)/lack_of_specificity.tab)
gene2sizes     = $(foreach c,$(CHILDREN),$(c)/gene2sizes.tab)
targets        = data.tab lists_t.tab

include $(MAPDIR)/Templates/Make/quick.mak

data.tab: $(lacks) $(gene2sizes)
	rm -f cover.tmp; \
	rm -f lack.tmp; \
	$(foreach c,$(CHILDREN),\
	   wc -l $(c)/gene2sizes.tab \
	   | cut -f 1 -d ' ' \
	   | cap.pl $(c) \
	   | transpose.pl -q \
	   >> cover.tmp; \
	   projection.pl -f Median $(c)/lack_of_specificity.tab \
	   | tail +2 \
	   | cap.pl $(c) \
	   | transpose.pl -q \
	   >> lack.tmp; \
	) \
	paste cover.tmp lack.tmp \
	| cut -f 1,2,4 \
	| sed 's/NR//' \
	| cap.pl Threshold,Coverage,LackOfSpecificity \
	> $@; \

lists_t.tab: ../../lists_t.tab
	cat $< \
	| grep -v '^#' \
	| rand_lines.pl \
	> $@; \

Overlaps: ../../Overlaps
	rm -f $@; \
	ln -s $< $@; \

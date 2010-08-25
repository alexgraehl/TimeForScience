root           = ../..
org            = $(notdir $(shell pwd))
org_dir        = $(root)/$(org)
CHILDREN       = $(shell grep '^ORGANISMS\>' $(MAPDIR)/Makefile.common | cut -f 2 -d = | sed 's/$(org)//g')
CHILD_MAKEFILE = $(HOME)/Map/Templates/Make/genesets_ortho_grandchild.mak
ortho_lists    = $(addsuffix /lists.tab,$(CHILDREN))
lists          = $(wildcard $(org_dir)/lists.tab)
targets        = lists_t.tab lists.tab



lists.tab: lists_t.tab
	sets.pl $< \
	> $@; \

lists_t.tab: $(ortho_lists) $(lists)
	cat $(ortho_lists) $(lists) \
	| sets.pl \
	> $@; \


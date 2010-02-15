root           = ../..
org            = $(notdir $(shell pwd))
org_dir        = $(root)/$(org)
CHILDREN       = $(foreach o,$(shell grep '^ORGANISMS\>' $(MAPDIR)/Makefile.common | cut -f 2 -d = | sed 's/$(org)//g'),$(wildcard $(root)/$(o)/data*.tab))
CHILDREN       = Fly Hpylori Human Mouse Yeast
CHILD_MAKEFILE = $(HOME)/Map/Templates/Make/interactions_ortho_grandchild.mak
ortho_links    = $(addsuffix /data.tab,$(CHILDREN))
org_links      = $(wildcard $(org_dir)/data.tab)
targets        = data.tab

include $(MAPDIR)/Templates/Make/quick.mak

a:
	echo '$(CHILDREN)'

data.tab: $(ortho_links) $(org_links)
	cat $(ortho_links) $(org_links) \
	| order_keys.pl \
	| sort -k 1,2 -t '	' -u \
	> $@; \


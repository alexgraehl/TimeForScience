input    = ../lists_t.tab
targets  = data.fa

include $(MAPDIR)/Templates/Make/quick.mak

data.fa: $(input)
	sets_overlap.pl -do '	'  -p .99  -q $< $< \
	> $@; \
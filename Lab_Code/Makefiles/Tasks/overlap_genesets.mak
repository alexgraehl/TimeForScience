
# Significance cutoff
pval_cut = 0.05

in       = ../lists_t.tab
geneset  = $(GENESETS)/$(THISDIR)/$(ORGANISM)/lists_t.tab
targets  = data.fa data.tab

include $(MAPDIR)/Templates/Make/quick.mak

data.tab: data.fa
	cat $< \
	| fasta2tab.pl \
	| sed 's/\\,/;/g' \
	| sort -k 2,2 -t ',' -g \
	> $@; \

data.fa: $(in) $(geneset)
	sets_overlap.pl $(geneset) $< \
	  -p $(pval_cut) \
	> $@; \


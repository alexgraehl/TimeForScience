# file: modes.mak
# last edit: chrisw on 25APR08
# This template is a Makefile for running MODES on a network to get modules.  The network to use is indicated by the variable $(in).  It should be a tab-delimited file with at least 2 colummns, one for each node in the edge.  The file can also have comment lines indicated by a leading "#" character.  A comment line may be used for a header line, for example.  All comment lines will be ignored in processing.

# Also note that the ordering of links in the file seems to have an effect on the modules output by MODES. Although it's not confirmed with the authors of MODES, this is probably due to the way MODES finds seeds to build modules from.

SPECIES       = $(GRANDDIR)
in            = ../data.tab.gz
modes_exe     = /projects/sysbio/apps/x86_64/bin/MODES/bin/modes
num_genes     = $(shell zcat $(in) | grep -v '^\#' | cut -f 2 | gzip | zcat - $(in) | grep -v '^\#' | cut -f 1 | sort -u | wc -l)
num_edges     = $(shell zcat $(in) | grep -v '^\#' | wc -l)
d_modes_param = 0.25
s_modes_param = 100
c_modes_param = 0.90
g_modes_param = 4

targets  = data.tab stats.tab lists_t.tab modules_t.tab data50.tab lists_t_names.tab

include $(MAPDIR)/Templates/Make/quick.mak

stats.tab: data.tab
	cat $< \
	| cap.pl dummy \
	| row_stats.pl -count \
	| transpose.pl -q \
	| row_stats.pl -count -mean -median -std \
	| transpose.pl -q \
	| sed -e '1d' \
	| cap.pl Stat,Value \
	> $@; \

# Link for compatability
lists_t.tab: data.tab
	rm -f $@; \
	ln -s $< $@; \

# Link for compatability
modules_t.tab: data.tab
	rm -f $@; \
	ln -s $< $@; \

# Create an adjacency matrix representation of the input network, then feed it to MODES, then map the MODES output back to genes.
data.tab: $(in)
	zcat $< \
	| grep -v '^\#' \
	> data.tmp;
	\
	cut.pl -f 2,1 data.tmp \
	| cat - data.tmp \
	| cut -f 1,2 \
	| expand.pl \
	| sets.pl -o m -e 0 -m 1 \
	> 1.tmp;
	\
	head -n 1 1.tmp \
	| transpose.pl -q \
	| join.pl - 1.tmp \
	| tee adj_matrix.tmp \
	| cut -f 2- \
	| sed -e '1d' \
	> adj_matrix_no_headers.tmp;
	\
	$(modes_exe) \
	   -m 0 \
	   -i adj_matrix_no_headers.tmp \
	   -n $(num_genes) \
	   -o mods.tmp \
	   -g $(g_modes_param) \
	   -d $(d_modes_param) \
	   -s $(s_modes_param) \
	   -c $(c_modes_param); \
	cut -f 1 adj_matrix.tmp \
	| sed -e '1d' \
	| lin.pl -0 \
	> index2id.tmp;
	\
	cat mods.tmp \
	| delim2delim.pl -f 1,2 '	' _ \
	| flatten.pl \
	| join.pl -1 2 -2 1 - index2id.tmp \
	| cut -f 2- \
	| expand.pl \
	| sed 's/^/Module/' \
	> $@;
	\
	rm data.tmp 1.tmp adj_matrix.tmp adj_matrix_no_headers.tmp index2id.tmp mods.tmp ; \

d.tab: $(in)
	zcat $< \
	| grep -v '^\#' \
	| cut -f 1,2 \
	> data.tmp; \
	cut -f 2 data.tmp \
	| cat - data.tmp \
	| cut -f 1 \
	| sort -u \
	| lin.pl -0 -a \
	> index2id.tmp; \
	cut -f 1 data.tmp \
	| join.pl - index2id.tmp \
	| cut -f 2 \
	> left.tmp; \
	cut -f 2 data.tmp \
	| join.pl - index2id.tmp \
	| cut -f 2 \
	> right.tmp; \
	paste left.tmp right.tmp \
	> in.tmp; \
	$(modes_exe) \
	   -m 0 \
	   -i in.tmp \
	   -y $(num_edges) \
	   -o mods.tmp \
	   -g $(g_modes_param) \
	   -d $(d_modes_param) \
	   -s $(s_modes_param) \
	   -c $(c_modes_param) \
	; \
	cat mods.tmp \
	| delim2delim.pl -f 1,2 '	' _ \
	| flatten.pl \
	| join.pl -1 2 -2 1 - index2id.tmp \
	| cut -f 2- \
	| expand.pl \
	| sed 's/^/Module/' \
	> $@; \
	rm *.tmp; \

data50.tab: $(in)
	zcat $< \
	| grep -v '^\#' \
	> data.tmp; \
	cut.pl -f 2,1 data.tmp \
	| cat - data.tmp \
	| cut -f 1,2 \
	| expand.pl \
	| sets.pl -o m -e 0 -m 1 \
	> 1.tmp; \
	head -n 1 1.tmp \
	| transpose.pl -q \
	| join.pl - 1.tmp \
	| tee adj_matrix.tmp \
	| cut -f 2- \
	| sed -e '1d' \
	> adj_matrix_no_headers.tmp; \
	$(modes_exe) \
	   -m 0 \
	   -i adj_matrix_no_headers.tmp \
	   -n $(num_genes) \
	   -o mods.tmp \
	   -g $(g_modes_param) \
	   -d 0.50 \
	   -s $(s_modes_param) \
	   -c $(c_modes_param); \
	cut -f 1 adj_matrix.tmp \
	| sed -e '1d' \
	| lin.pl -0 \
	> index2id.tmp; \
	cat mods.tmp \
	| delim2delim.pl -f 1,2 '	' _ \
	| flatten.pl \
	| join.pl -1 2 -2 1 - index2id.tmp \
	| cut -f 2- \
	| expand.pl \
	| sed 's/^/Module/' \
	> $@; \
	rm -f *.tmp; \

lists_t_names.tab: lists_t.tab $(HOME)/Map/Data/Gene/Description/$(SPECIES)/data.tab
	cat $< \
	| sets.pl \
	| cut -f 1 \
	| join.pl - $(HOME)/Map/Data/Gene/Description/$(SPECIES)/data.tab -o NONE \
	| paste.pl - $< \
	| cut.pl -f 2,1 \
	| sed 's/^NONE	//' \
	| sed 's/^	//' \
	| cut -f 1 \
	> 1.tmp; \
	\
	cat $< \
	| sets.pl \
	| paste.pl 1.tmp - \
	| cut -f 1,3 \
	| sets.pl \
	| sort -u \
	> $@; \


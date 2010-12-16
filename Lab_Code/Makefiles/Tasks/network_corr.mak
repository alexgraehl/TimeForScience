in             = ../../data.tab
num_expts      = $(shell head -n 1 $(in) | transpose.pl -q | sed -e '1d' | wc -l | sed 's/^[ ][ ]*//' | cut -f 1 -d ' ')
# The minimum number of experiments for Spearman calculation.
min_expts      = 30
max_top        = 500
tops           = 50
max_top_tab_gz = top$(max_top).tab.gz
top_tab_gzs    = $(foreach t,$(tops),top$(t).tab.gz)
# num_rand     = 10000
num_rand       = 10001
stdevs         = 4
norm_ops       = -trim 80 -normcol -normrow -flat 25
#norm_ops       = -trim 80 -normrow -flat 25
data_ops       = -topk $(max_top) $(norm_ops)
rand_ops       = -permcol -arb $(num_rand) $(norm_ops)
pearson_exe    = /cse/faculty/jstuart/src/stuartlab/genie_release/Programs/$(MACHTYPE)/pearson
pearson        = nice $(pearson_exe) $(if $(findstring -,$(shell echo '$(min_expts)-$(num_expts)' | bc)),-spearman,)
targets        = top500.tab.gz nodes.tab probs.tab 

include /projects/sysbio/map/Templates/Make/quick.mak

a:
	echo $(max_top_tab_gz); \

$(max_top_tab_gz): $(in)
	rm -f run.txt; \
	tail -n +2 $< \
	| wc -l \
	| sed 's/^[ ]*//' \
	| cap.pl 'No. Genes =' \
	| transpose.pl -q \
	| sed 's/\t/ /g' \
	>> run.txt; \
	echo 'No. Experiments = $(num_expts)' \
	>> run.txt; \
	echo 'Pearson command = $(pearson)' \
	>> run.txt; \
	echo 'Standard Deviations = $(stdevs)' \
	>> run.txt; \
	echo "$(pearson) $< $(data_ops)" \
	>> run.txt; \
	date >> run.txt; \
	\
	$(pearson) $< $(data_ops) \
	| order_keys.pl \
	| gzip \
	> $@; \
	date >> run.txt; \

nodes.tab: top500.tab.gz
	zcat $< \
	| tail -n +2 \
	> data.tmp; \
	cat data.tmp \
	| cut -f 2 \
	| cat data.tmp - \
	| cut -f 1 \
	| sort -u \
	> $@; \
	rm -rf data.tmp; \

probs.tab: nodes.tab top500.tab.gz
	wc -l $<  \
	| cut -f 1 -d ' ' \
	> stats.tmp; \
	zcat top500.tab.gz \
	| sed -e 1d  \
	| cut -f 1,2 \
	| sort -u \
	| wc -l \
	| cut -f 1 -d ' ' \
	>> stats.tmp; \
	cat stats.tmp \
	| transpose.pl \
	| perl -ne 'chomp; @tabs = split (/\t/); $$den = ($$tabs[0]*($$tabs[0]-1))/2; $$con = $$tabs[1]/$$den; print "$$tabs[0]\t$$tabs[1]\t$$con\n";' \
	| cap.pl Nodes,Links,Connectivity \
	> $@; \

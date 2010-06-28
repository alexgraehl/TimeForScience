in          = ../../data.tab
num_expts   = $(shell head -n 1 $(in) | transpose.pl -q | tail +2 | wc -l | sed 's/^[ ][ ]*//' | cut -f 1 -d ' ')
# The minimum number of experiments for Spearman calculation.
min_expts   = 30
topk        = 500
# num_rand   = 10000
num_rand    = 10001
stdevs      = 4
# norm_ops   = -trim 80 -normcol -normrow -flat 25
norm_ops    = -trim 80 -normrow -flat 25
data_ops    = -topk $(topk) $(norm_ops)
rand_ops    = -permrow -arb $(num_rand) $(norm_ops)
pearson_exe = /cse/faculty/jstuart/src/stuartlab/genie_release/Programs/$(MACHTYPE)/pearson -abs
pearson     = nice $(pearson_exe) $(if $(findstring -,$(shell echo '$(min_expts)-$(num_expts)' | bc)),-spearman,)
targets     = data.tab.gz top.tab.gz stats.tab

include $(HOME)/Map/Templates/Make/quick.mak

data.tab.gz: top.tab.gz stats.tab
	zcat $< \
	| head -n 1 \
	> header.tmp; \
	zcat $< \
	| tail +2 \
	| select.exe -abs -k 5 -gte $(shell tail +2 stats.tab | cut -f 1) \
	| cut.pl -f 1-2,1- \
	| order_keys.pl \
	| uniq.pl -k 1,2 \
	| cut -f 3- \
	| cat header.tmp - \
	| gzip \
	> $@; \
	rm *.tmp; \

top.tab.gz: $(in) stats.tab
	rm -f run.txt; \
	tail +2 $< \
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
	cat stats.tab \
	>> run.txt; \
	echo "$(pearson) $< $(data_ops)" \
	>> run.txt; \
	date >> run.txt; \
	\
	$(pearson) $< $(data_ops) \
	| gzip \
	> $@; \
	date >> run.txt; \

stats.tab: $(in)
	$(pearson) $< $(rand_ops) \
	| tail +2 \
	| cut -f 5 \
	| cap.pl dummy \
	| transpose.pl -q \
	| cap.pl dummy \
	| row_stats.pl -mean -std \
	| tail +2 \
	| cut -f 2- \
	| awk '{print $$1+$(stdevs)*$$2"\t"$$1"\t"$$2;}' \
	| cap.pl Cutoff,Mean,StdDev \
	> $@; \



in       = ../../data.tab
genes    = ../../Filter/genes.lst
expts    = ../../Filter/expts.lst
# R_exe  = /cluster/home/sugnet/local/kolossus/bin/R
R_exe    = /cluster/home/kpollard/R/R-2.0.1/bin/R
r_templ  = $(HOME)/Map/Templates/R/runHopachM.R
targets  = data.tab hopach.out run.R

all: $(targets)

clean:
	rm -f $(targets) date.txt

data.tab: hopach.out $(in)
	cut -f 1 $(in) \
	| paste - hopach.out \
	| tail +2 \
	| sort -k 2,2 -t '	' -n \
	> $@; \
	sort -u hopach.out \
	| lin.pl -a \
	| join.pl -1 2 $@ - \
	| cut -f 2- \
	> $@.tmp; \
	mv $@.tmp $@; \

hopach.out: run.R
	projection.pl -f $(expts) $(in) \
	> 1.tmp; \
	cut -f 1 $(in) \
	| paste - 1.tmp \
	> 2.tmp; \
	head -n 1 2.tmp \
	> in.tmp; \
	join.pl $(genes) 2.tmp \
	>> in.tmp; \
	\
	date > date.txt; \
	cat run.R \
	| $(R_exe) --no-save; \
	date >> date.txt; \
	rm -f *.tmp; \

run.R: $(in) $(r_templ)
	bind.pl $(r_templ) input_file=in.tmp \
	> $@; \



root       = ../../..
org1       = $(notdir $(shell cd ..; pwd))
org2       = $(notdir $(shell pwd))
blast_dir  = $(HOME)/Map/Data/Blast/Protein/Result
blast      = $(blast_dir)/$(org1)/$(org2)/best.tab
org2_dir   = $(root)/$(org2)
lists      = $(org2_dir)/lists.tab
targets    = lists.tab lists_t.tab

all: $(targets)

clean:
	rm -f $(targets) $(wildcard *.tmp)

make:

include $(HOME)/Map/Templates/Make/quick.mak

lists_t.tab: lists.tab
	sets.pl -h 0 $< \
	> $@; \

lists.tab: $(lists) $(blast)
	cut.pl -f 2,1 $(blast) \
	| join.pl - $< \
	| cut -f 2- \
	> $@; \



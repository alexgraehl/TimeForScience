root       = ../../..
org1       = $(notdir $(shell cd ..; pwd))
org2       = $(notdir $(shell pwd))
blast_dir  = $(HOME)/Map/Data/Blast/Protein/Result
blast      = $(blast_dir)/$(org1)/$(org2)/best.tab
org2_dir   = $(root)/$(org2)
links      = $(org2_dir)/data.tab
targets    = data.tab

all: $(targets)

clean:
	rm -f $(targets) $(wildcard *.tmp)

make:

include $(HOME)/Map/Templates/Make/quick.mak

data.tab: $(links) $(blast)
	cut.pl -f 2,1 $(blast) \
	> blast.tmp; \
	cut -f 1 $< \
	| join.pl -o NONE - blast.tmp \
	| cut -f 2 \
	> left.tmp; \
	cut -f 2 $< \
	| join.pl -o NONE - blast.tmp \
	| cut -f 2 \
	> right.tmp; \
	paste left.tmp right.tmp \
	| grep -v NONE \
	| order_keys.pl \
	| sort -k 1,2 -t '	' -u \
	> $@; \
	rm blast.tmp left.tmp right.tmp; \



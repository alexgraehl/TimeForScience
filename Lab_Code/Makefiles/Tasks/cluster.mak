empty      =
space      = $(empty) $(empty)
root       = ../..
genes_file = $(root)/Filter/genes.lst
expts_file = $(root)/Filter/expts.lst
input      = $(root)/data.tab
method     = $(word 1, $(subst _,$(space),$(notdir $(shell pwd))))
metric     = $(word 2, $(subst _,$(space),$(notdir $(shell pwd))))
clusters   = $(word 3, $(subst _,$(space),$(notdir $(shell pwd))))
em         = false
genexpress = true
template   = $(HOME)/Map/Templates/Runs/cluster.map

ifneq (,$(shell exists.pl $(genes_file)))
   genes = $(genes_file)
else
   genes = '#undef'
endif

ifneq (,$(shell exists.pl $(expts_file)))
   expts = $(expts_file)
else
   expts = '#undef'
endif

targets    = data.gxp data.tab

all: $(targets)

clean:
	rm -f $(targets) map.log date.txt map.xml $(wildcard *.tmp)

data.tab: data.gxp
	grep '<Gene Id' $< \
	| sed 's/^..*ORF=//' \
	| sed 's/"/\t/g' \
	| cut -f 2 \
	> 1.tmp; \
	grep '<Attributes AttributesGroupId' $< \
	| sed 's/^..*Value=//' \
	| sed 's/"/\t/g' \
	| cut -f 2 \
	| paste 1.tmp - \
	> $@; \
	rm *.tmp; \

data.gxp: $(input)
	make map.xml; \
	date > date.txt; \
	map_learn map.xml; \
	date >> date.txt; \

map.xml: $(input)
	bind.pl $(template) \
	   output_file=data.gxp \
	   expression_file=$< \
	   num_clusters=$(clusters) \
	   cluster_method=$(method) \
	   cluster_metric=$(metric) \
	   print_genexpress=$(genexpress) \
	   gene_list=$(genes) \
	   experiment_list=$(expts) \
	   em=$(em) \
	   em=$(em) \
	> $@; \



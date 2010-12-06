include $(HOME)/Map/Makefile.common
include $(HOME)/Map/Data/Makefile.common

ifneq (,$(shell exists.pl ./Makefile.common))
include ./Makefile.common
endif

ifndef ($(platform_regex))
   platform_regex = '.'
endif

ifndef ($(root))
   root = ..
endif

platform_pubs  = $(shell grep $(platform_regex) $(root)/*/desc.tab \
                         | grep '^experiment_type' \
                         | cut -f 1 \
                         | delim2delim.pl '/' '	' \
                         | cut -f 3 \
                         | sort -u \
                  )


ifeq (,$(shell exists.pl publications.lst))
  sel_pubs  = $(notdir $(shell find $(root) -type d -maxdepth 1 \
                                | grep -v 'Compendium' \
                                | grep -v '\.\.$$' \
                                | grep -v '\.\.\/$$' \
                                | grep '\/' \
                         ) \
                )
else
  sel_pubs  = $(shell cat publications.lst)
endif

pubs   = $(shell echo $(platform_pubs)          \
                 | space2tab.pl            \
                 | transpose.pl -q         \
                 > $(tmp);                 \
                 echo $(sel_pubs)          \
                 | space2tab.pl            \
                 | transpose.pl -q         \
                 | cat $(tmp) -            \
                 | sort -u                 \
                 > generated_pubs.lst;     \
                 $(RM) *$(tmp);            \
                 cat generated_pubs.lst;   \
           )

data_files   = $(foreach pub, $(pubs), $(root)/$(pub)/data.tab)
targets      = data.tab

all: $(targets)

clean:
	$(RM) $(targets)

data.tab: $(data_files) header.lst
	cat header.lst \
	| transpose.pl -q \
	> data.tab; \
	join_multi.pl $(data_files) \
	| body.pl 2 -1 \
	>> data.tab; \

header.lst: $(data_files)
	echo 'Gene' \
	> header.lst; \
	$(foreach pub, $(pubs), \
	   head -n 1 $(root)/$(pub)/data.tab \
	   | transpose.pl -q \
	   | body.pl 2 -1 \
	   | concat.pl -p '$(ORGANISM) $(pub) ' -sf - \
	   >> header.lst; \
	) \



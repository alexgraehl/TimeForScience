
ifeq (,$(doc))
doc = doc
endif

styles       = sty cls tex
style_files  = $(wildcard $(addprefix $(style_dir)/,$(addprefix *.,$(styles))))
style_links  = $(notdir $(style_files))
side_effects = $(doc).aux $(doc).log
targets      = $(doc).pdf $(style_links) $(fig_links)

include $(MAPDIR)/Templates/Make/quick.mak

$(doc).pdf: $(doc).tex $(style_links) $(fig_links)
	pdflatex $(basename $<); \

$(style_links): %: $(style_dir)/%
	rm -f $@; \
	ln -s $< $@; \

$(doc).tex:
	echo 'You need to write this; see $(wildcard $(style_dir)/*.tex) for an example.';

ifneq (,$(fig_links))

fig_nums = $(shell nums.pl $(words $(fig_links)))

$(word 1, $(fig_links)): $(fig_files)
	$(foreach i,$(fig_nums), \
	   rm -f $(word $(i),$(fig_links)); \
	   ln -s $(word $(i),$(fig_files)) $(word $(i),$(fig_links)); \
	) \

endif


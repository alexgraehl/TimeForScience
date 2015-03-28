
in           = ../desc.tab
in2          = ../data.tab
fields       = title description organism sample_organism
gsms         = $(shell grep '^GSM' $(in) | cut -f 1)
syn          = $(SYSBIODIR)/users/dsam/src/perl/Map/Expression/syn.tab
# mmtx_exe   = MMTx --KSYear=2006 -anu --best_mappings_only -I --show_treecodes -x -Q
mmtx_exe     = MMTx --KSYear=2006 -anu --show_treecodes --show_cuis -x -Q
clean_exe    = perl $(SYSBIODIR)/users/dsam/src/perl/Map/Expression/clean_geo_desc_for_mmtx.pl -syn $(syn)
txtToTab_exe = perl $(SYSBIODIR)/users/dsam/src/perl/Map/Expression/geo_umls_txt_to_tab.pl
cluster_exe  = cluster-eisen -g 1 -e 1 -m a
threshold    = 600
targets      = data.tab.gz data.txt.gz data_cleaned.cdt data_cleaned.list data_cleaned_semgrp.list

include $(MAPDIR)/Templates/Make/quick.mak


data.tab.gz: $(in)
	rm -f $@; \
	touch $(basename $@); \
	$(foreach x, $(fields) $(gsms), \
	   grep '^$(x)	' $< \
	   | cut -f 2- \
	   | $(mmtx_exe) -f \
	   | paste.pl $(x) - \
	   >> $(basename $@); \
	) \
	gzip $(basename $@); \

data.txt.gz: $(in)
	rm -f $@; \
	touch $(basename $@); \
	$(foreach x, $(fields) $(gsms), \
	   grep '^$(x)	' $< \
	   | cut -f 2- \
	   | $(mmtx_exe) \
	   | paste.pl $(x) - \
	   >> $(basename $@); \
	) \
	gzip $(basename $@); \

# Generate a clean version of MMtx UMLS mapping
# - run perl script to clean up description file
#   e.g. add spaces to non-word characters
# - retrieve both phrase and single word mappings
data_cleaned.txt.gz: $(in) $(in2)
	rm -f $@ $(basename $@); \
	touch $(basename $@); \
	head -1 $(in2) \
	| flatten.pl \
	| cut.pl -f 2 \
	| delim2delim.pl ";" "	" \
	| flatten.pl \
	| cat $(in) - \
	> in.tmp ; \
	$(clean_exe) in.tmp \
	| uniq.pl -f 1- \
	| sed 's/\t\s*/\t/g' \
	> cleaned_data.in; \
	$(foreach x, $(fields) $(gsms), \
	   rm -f in1.tmp in2.tmp ; \
	   grep '^$(x)	' cleaned_data.in \
	   | cut -f 2- \
	   > in1.tmp ; \
	   grep '^$(x)	' cleaned_data.in \
	   | delim2delim.pl " " " .	" \
	   | flatten.pl \
	   | cut -f 2- \
	   > in2.tmp ; \
	   echo " . " \
	   | cat - in2.tmp \
	   | cat in1.tmp - \
	   | $(mmtx_exe) \
	   | paste.pl $(x) - \
	   >> $(basename $@); \
	) \
	gzip $(basename $@); \
	rm -f *.tmp ; \

# Convert txt output of MMTx to tab-delimited form
data_cleaned.tab.gz: data_cleaned.txt.gz
	$(txtToTab_exe) $< \
	| gzip \
	> $@; \

# Foreach term, get GSMs assigned with that term above threshold
# add organism as a term
data_cleaned.list: data_cleaned.tab.gz $(in)
	zcat $< \
	| sort -t '	' -k 2,2 -g -r \
	| select.pl -k 2 -gte $(threshold) \
	| uniq.pl -f 1,4 \
	> in.tmp ; \
	grep "^sample_organism" in.tmp \
	| cut.pl -f 4 \
	> org.tmp ; \
	add_column.pl org.tmp -s "$(shell grep "^GSM" $(in) | cut.pl -f 1 | transpose.pl)" \
	> org2.tmp ; \
	grep "^GSM" in.tmp \
	| cut.pl -f 4,1 \
	| uniq.pl -f 1,2 \
	| expand.pl \
	| cat org2.tmp - \
	| uniq.pl -f 1 \
	> $@; \
	rm -f *.tmp ; \

# Generate clustered heatmap of terms assigned to GSMs
data_cleaned.cdt: data_cleaned.list
	flatten.pl data_cleaned.list \
	> in.tmp ; \
	add_column.pl in.tmp -s '1' \
	| sets.pl -i p -o m -m 3 -e 0 -h 0 \
	> data_cleaned.tmp ; \
	$(cluster_exe) -f data_cleaned.tmp ; \
	rm -f *.tmp ; \

# Retrieve UMLS term with semantic network grouping
# - for GSM and organism only
# - terms about threshold
data_cleaned_semgrp.list: data_cleaned.tab.gz
	zcat $< \
	| tail -n +2 \
	> in.tmp ; \
	grep "^GSM" in.tmp \
	> gsm.tmp ; \
	grep "^sample_organism" in.tmp \
	| cat - gsm.tmp \
	| select.pl -k 2 -gte $(threshold) \
	| cut.pl -f 4,6 \
	| uniq.pl -f 1,2 \
	| expand.pl \
	> 1.tmp ; \
	cut.pl -f 2- 1.tmp \
	| perl -i -p -e 's/\t/\|/g' \
	| paste - 1.tmp \
	| cut.pl -f 2,1 \
	> $@ ; \


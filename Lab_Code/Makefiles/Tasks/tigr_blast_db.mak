include $(HOME)/Map/Makefile.common
include $(HOME)/Map/Data/Makefile.common

URL           = $(TIGR_FTP_GENE_SEQ_$(ORGANISM_UPPER))
targets       = data.nsq

all: $(targets)

clean:
	$(RM) $(wildcard *.tmp) $(targets)

include $(HOME)/Map/Templates/Make/data.mak

tigr_abbrev   = $(TIGR_GENE_INDEX_$(ORGANISM_UPPER)_ABBREV)
ego_abbrev    = $(EGO_SEQ_$(ORGANISM_UPPER)_ABBREV)
unzipped      = $(shell unzip -l $(REMOTE_TARGET) \
		  | grep '$(tigr_abbrev)\.[0-9]' \
		  | sed 's/^[ ][ ]*//' \
		  | space2tab.pl \
		  | cut -f 4 \
		)

data.nsq: $(REMOTE_TARGET)
	unzip $(REMOTE_TARGET) $(unzipped); \
	cat $(unzipped) \
	| fasta2stab.pl \
	| cut.pl -f 2,1 \
	| space2tab.pl \
	| cut.pl -f 2,1 \
	| stab2fasta.pl \
	> data; \
	formatdb -i data -p F; \
	$(RM) data $(unzipped); \


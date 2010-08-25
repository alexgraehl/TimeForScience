include $(MAPDIR)/Makefile.common
include $(MAPDIR)/Data/Makefile.common

targets = data.tab.gz best.tab

all: $(targets)

clean:
	rm -f $(targets)

make:

remote:
	cd $(FASTA_QUERY_DIR); make remote; make; \
	cd $(BLAST_DATABASE_DIR); make remote; make; \

remote_clean:

back:
	mkdir -p Backup; \
	cp $(targets) Backup; \

back_clean:
	rm -rf Backup; \

best.tab: data.tab.gz
	zcat $< \
	| uniq.pl \
	> best.tab; \


data.tab.gz: $(FASTA_QUERY) $(BLAST_DATABASE_FILES)
	$(BLAST) \
	| $(BLAST_SORT) \
	| gzip \
	> $@; \


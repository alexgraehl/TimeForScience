include $(MAPDIR)/Makefile.common
include $(MAPDIR)/Data/Makefile.common

# File containing multiple FASTA sequences
FASTA_DIR  = $(MAP_DATA)/Protein/$(ORGANISM)
FASTA_FILE = $(FASTA_DIR)/$(FASTA_PROTEIN_FILE)

all: data.phr

clean:
	$(CLEAN_BLAST_DATABASE)

remote:
	cd $(FASTA_DIR); \
	make remote; \
	make

remote_clean:

data.phr: $(FASTA_FILE)
	$(FORMATDB); \
	rename.pl 's/$(FASTA_PROTEIN_FILE)/data/' $(FASTA_PROTEIN_FILE).*


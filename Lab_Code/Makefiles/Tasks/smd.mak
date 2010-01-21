include $(HOME)/Map/Makefile.common
include $(HOME)/Map/Data/Makefile.common
include $(HOME)/Map/Papers/ExpressionNetwork/Makefile.common

smd_archives = $(wildcard $(REMOTE_DIR)/*.tar*)
smd_archive1 = $(word 1, $(smd_archives))
smd_files    = $(wildcard $(REMOTE_DIR)/*_$(SMD_ARRAY_FILE_SUFFIX))
smd_file1    = $(word 1, $(smd_files))

ifneq (,$(shell exists.pl Makefile.common))
   include Makefile.common
endif

ifndef ($(DATA_FIELD))
   DATA_FIELD = LOG_RAT2N_MEAN
endif

ifndef ($(ID_FIELD))
   ID_FIELD = 'Cluster ID'
endif

ifndef ($(GENE_FIELD))
   GENE_FIELD = 'Gene Name'
endif

ifndef ($(PUB))
   PUB = 0
endif

URL      = ftp://genome-ftp.stanford.edu/pub/smd/publications/$(PUB)
targets  = data.tab experiments.lst


all: $(targets)


include $(HOME)/Map/Templates/Make/data.mak

remote:
	$(MKDIR) $(REMOTE_DIR); \
	$(foreach u, $(URL), \
	   $(WGET_ALL) -P $(REMOTE_DIR) '$(u)'; \
	) \
	make unpack; \

slim:
	$(RMDIR) $(REMOTE_DIR)

clean:
	$(RM) $(wildcard *.tmp) $(targets)

ifneq (,$(smd_file1))

unpack: $(smd_file1)

$(smd_file1): $(smd_archive1)
	remote_smd.pl $(REMOTE_DIR) $(SMD_ARRAY_FILE_SUFFIX); \

else

unpack:
	remote_smd.pl $(REMOTE_DIR) $(SMD_ARRAY_FILE_SUFFIX); \

endif

data.tab: $(REMOTE_DIR)/publication_$(PUB).meta
	smdxls2tab.pl $(ID_FIELD) $(GENE_FIELD) $(DATA_FIELD) $(REMOTE_DIR) \
	> data.tab; \

experiments.lst: data.tab
	head -n 1 data.tab \
	| transpose.pl -q \
	| body.pl 3 -1 \
	> experiments.lst \



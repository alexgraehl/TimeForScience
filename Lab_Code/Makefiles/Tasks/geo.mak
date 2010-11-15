include $(MAPDIR)/Makefile.common
include $(MAPDIR)/Data/Makefile.common

TIMESTAMPING    = OFF

REMOTE_DIR      = $(MAP_DATA)/Expression/Any/Geo/Remote
# platform_url    = ftp://ftp.ncbi.nih.gov/pub/geo/data/geo/soft
# dataset_url     = ftp://ftp.ncbi.nih.gov/pub/geo/data/gds/soft_gz
# series_url     = ftp://ftp.ncbi.nih.gov/pub/geo/DATA/SOFT/by_series
platform_url    = ftp://ftp.ncbi.nih.gov/pub/geo/DATA/SOFT/by_platform
series_url     = ftp://ftp.ncbi.nih.gov/pub/geo/DATA/SOFT/GDS
datasets        = $(shell cat datasets.fa | fasta2tab.pl | cut -f 1,3- | flatten.pl | cut -f 2 | sort -u)
dataset1        = $(word 1, $(datasets))
dataset_suffix  = .soft.gz
platforms       = $(shell cat datasets.fa | fasta2tab.pl | cut -f 1 | sort -u)
platform1       = $(word 1, $(platforms))
platform_suffix = _family.soft.gz
targets         = data.tab datasets.fa

include $(MAPDIR)/Templates/Make/quick.mak

clean:
	rm -f data.tab

remote:
	mkdir -p $(REMOTE_DIR); \
	$(foreach p, $(platforms), \
	   $(REMOTE_FTP) -P $(REMOTE_DIR) $(platform_url)/$(p)/$(p)$(platform_suffix); \
	) \
	$(foreach d, $(datasets), \
	   $(REMOTE_FTP) -P $(REMOTE_DIR) $(series_url)/$(d)$(dataset_suffix); \
	) \

data.tab: $(REMOTE_DIR)/$(platform1)$(platform_suffix) datasets.fa
	parse_geo_gds.pl -o $(ORGANISM) \
	   -d desc.tab \
	   $(REMOTE_DIR) datasets.fa \
	| sed 's/NaN//g' \
	> data.tmp; \
	cat data.tmp \
	| cut -f 1 \
	> orig_genes.txt; \
	cat orig_genes.txt \
	| perl -ne 's/; .*//g; print;' \
	| paste.pl - data.tmp \
	| cut -f 1,3- \
	> $@; \
	rm -rf data.tmp; \

datasets.fa:
	echo 'You need to make this file.'


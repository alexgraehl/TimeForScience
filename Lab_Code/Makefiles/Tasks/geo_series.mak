TIMESTAMPING    = OFF

REMOTE_DIR      = $(MAP_DATA)/Expression/Any/Geo/Remote
series_url     = ftp://ftp.ncbi.nih.gov/pub/geo/DATA/SOFT/by_series
platform_url    = ftp://ftp.ncbi.nih.gov/pub/geo/DATA/SOFT/by_platform
series        = $(shell cat series.fa | fasta2tab.pl | cut -f 1,3- | flatten.pl | cut -f 2 | sort -u)
series1        = $(word 1, $(series))
series_suffix  = _family.soft.gz
platforms       = $(shell cat series.fa | fasta2tab.pl | cut -f 1 | sort -u)
platform1       = $(word 1, $(platforms))
platform_suffix = _family.soft.gz
targets         = data.tab series.fa

include $(MAPDIR)/Templates/Make/quick.mak

clean:
	rm -f data.tab

remote:
	mkdir -p $(REMOTE_DIR); \
	$(foreach p, $(platforms), \
	   $(REMOTE_FTP) -P $(REMOTE_DIR) $(platform_url)/$(p)/$(p)$(platform_suffix); \
	) \
	$(foreach d, $(series), \
	   $(REMOTE_FTP) -P $(REMOTE_DIR) $(series_url)/$(d)/$(d)$(series_suffix); \
	) \

data.tab: $(REMOTE_DIR)/$(platform1)$(platform_suffix) series.fa
	parse_geo_gds.pl -o $(org_type) \
	   -d desc.tab \
	   $(REMOTE_DIR) series.fa \
	| sed 's/NaN//g' \
	> data.tab; \

series.fa:
	echo 'You need to make this file.'


include $(HOME)/Map/Makefile.common
include $(HOME)/Map/Data/Makefile.common

exe     = $(MAP_RELEASE)
xml     = $(MAP_TEMPLATES)/Runs/pearson.map
in      = ../data.tab
topk    = 50
mind    = 5
targets = data.tab t.tab

all: $(targets)

clean:
	$(RM) $(wildcard *.tmp) $(targets)

make:

include $(HOME)/Map/Templates/Make/data.mak

t.tab: data.tab
	cut -f 3,4 data.tab \
	| r2t.exe \
	| paste data.tab - \
	| cut -f 1,2,5 \
	> t.tab; \


data.tab: $(in) $(xml)
	date > time.txt; \
	bind.pl $(xml) top_k=$(topk) in=$(in) out=$(tmp) -exe $(exe); \
	cat $(tmp) \
	| grep -v '\-1.79769e+308' \
	| grep -v 'NaN' \
	| select.exe -k 4 -gte '$(mind)' \
	> data.tab; \
	date >> time.txt; \
	$(RM) $(tmp); \

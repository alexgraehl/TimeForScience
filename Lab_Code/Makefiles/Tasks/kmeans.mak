
clusters    = $(shell pwd | cut.pl -f -1 -d /)
distance    = $(shell pwd | cut.pl -f -2 -d / | tr A-Z a-z)
start       = sample
maxiter     = 30
emptyaction = singleton
display     = off
params      = 'Start','$(start)','Distance','$(distance)','Maxiter',$(maxiter)','Display','$(display)','EmptyAction','$(emptyaction)'

input       = ../../../../data.tab
targets     = data.tab startup.m data_t.tab

include $(HOME)/Map/Templates/Make/quick.mak

data_t.tab: data.tab
	sets.pl $< \
	> $@; \

data.tab: $(input) startup.m
	cat $(input) \
	| sed 's/\t\t/\tNaN\t/g' \
	| sed 's/\t\t/\tNaN\t/g' \
	| sed 's/\t$$/\tNaN/' \
	| tail +2 \
	| cut -f 2- \
	> data.tmp; \
	\
	cat $(input) \
	| cut -f 1 \
	| tail +2 \
	> genes.tmp; \
	\
	matlab -nojvm; \
	\
	rm genes.tmp data.tmp; \

startup.m: $(input)
	echo 'clear all;' \
	> $@; \
	echo "X = load('data.tmp');" \
	>> $@; \
	echo "genes = load_strings('genes.tmp');" \
	>> $@; \
	echo "I = kmeans(X,$(clusters),$(params));" \
	>> $@; \
	echo "printmatrix(I,genes,[],'data.tab','w');" \
	>> $@; \
	echo "exit;" \
	>> $@; \


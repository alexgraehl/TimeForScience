
in      = ../data.tab
targets = run.m X.dat Y.dat data.tab

include $(MAPDIR)/Templates/Make/quick.mak

data.tab: Y.dat $(in)
	head -n 1 $< \
	| transpose.pl -q \
	| lin.pl \
	| cut -f 1 \
	| cap.pl Key \
	| transpose.pl -q \
	> $@; \
	cut -f 1 $(in) \
	| tail +2 \
	| paste - $< \
	| sed 's/NaN//g' \
	>> $@; \

Y.dat: X.dat run.m
	cp run.m startup.m; \
	echo 'exit;' \
	>> startup.m; \
	/usr/local/bin/matlab -nojvm -glnx86; \
	rm startup.m; \

X.dat: $(in)
	KNNimpute $< 15 0 1.tmp; \
	cut -f 2- 1.tmp \
	| tail +2 \
	| fill_nan.pl \
	> X.dat; \
	rm 1.tmp; \

run.m: $(MAPDIR)/Templates/Matlab/decorrelate.m
	rm -f $@; \
	ln -s $< $@; \


in        = ../lists_t.tab
query_tol = $(shell echo '$(THISDIR)' | grep 'Nr' | sed 's/Nr\([0-9][0-9]*\)/\1/')
targets   = lists_t.tab big2small.tab small2big.tab stats.tab

include $(MAPDIR)/Templates/Make/quick.mak

ifeq (,$(query_tol))
query_tol = 50
endif

lists_t.tab: big2small.tab
	rm -f $@; \
	ln -s $< $@; \

big2small.tab: $(in)
	cat $< \
	| cap.pl dummy \
	| row_stats.pl -count \
	| tail +2 \
	| cut -f 2 \
	| paste - $< \
	| sort -k 1,1 -t '	' -n -r \
	| cut -f 2- \
	| nr_sets_josh.pl -qtol $(query_tol) \
	> $@; \
	$(RM) *.tmp; \

small2big.tab: $(in)
	cat $< \
	| cap.pl dummy \
	| row_stats.pl -count \
	| tail +2 \
	| cut -f 2 \
	| paste - $< \
	| sort -k 1,1 -t '	' -n \
	| cut -f 2- \
	| nr_sets_josh.pl -qtol $(query_tol) \
	> $@; \
	$(RM) *.tmp; \

# Coverage: The number of original members (genes) that are included in at least one set
# of the resulting non-redudant set.
# Specificity: (N-size)/(N-1) where N is the total number of members.
# RelativeSize: The 1/(size of the average set that any member (gene) is a member of).
# Higher specificities mean genes are in smaller clusters.
stats.tab: lists_t.tab $(in)
	wc -l $< \
	| cut -f 1 -d ' ' \
	| paste.pl - `wc -l $(in) | cut -f 1 -d ' '` \
	| awk '{print ($$2-$$1)/$$2*100;}' \
	| sed 's/^/Reduction\t/' \
	> $@; \
	flatten.pl $(in) \
	| cut -f 2 \
	| sort -u \
	| wc -l \
	| cut -f 1 \
	> old.tmp; \
	flatten.pl $< \
	| cut -f 2 \
	| sort -u \
	| wc -l \
	| cut -f 1 \
	> new.tmp; \
	paste new.tmp old.tmp \
	| awk '{print $$1/$$2*100;}' \
	| paste.pl 'Coverage' - \
	>> $@; \
	\
	flatten.pl $< \
	> 1.tmp; \
	join.pl 1.tmp $< \
	| cut -f 2- \
	| flatten.pl \
	| sed 's/^\([^\t][^\t]*\)\t\1/SELF\t\1\t\1/' \
	| grep -v '^SELF' \
	| expand.pl \
	| cap.pl dummy \
	| row_stats.pl -count \
	| tail +2 \
	| cut -f 2 \
	| paste.pl - `cat old.tmp` \
	| awk '{print ($$2-$$1)/($$2-1)*100;}' \
	| cap.pl dummy \
	| transpose.pl -q \
	| cap.pl dummy \
	| row_stats.pl -mean -std \
	| tail +2 \
	| cut -f 2- \
	| cap.pl Specificity_Mean,Specificity_Std \
	| transpose.pl -q \
	>> $@; \
	\
	cap.pl $< \
	| row_stats.pl -count \
	| tail +2 \
	> sizes.tmp; \
	flatten.pl $< \
	| join.pl -rev $< sizes.tmp \
	| sed 's/\t/@@@/' \
	| sets.pl \
	| flatten.pl \
	| sed 's/@@@/\t/' \
	| cut -f 1,3 \
	| expand.pl \
	| cap.pl dummy \
	| row_stats.pl -mean \
	| tail +2 \
	> new_ave_set.tmp; \
	\
	cap.pl $(in) \
	| row_stats.pl -count \
	| tail +2 \
	> sizes.tmp; \
	flatten.pl $(in) \
	| join.pl -rev $(in) sizes.tmp \
	| sed 's/\t/@@@/' \
	| sets.pl \
	| flatten.pl \
	| sed 's/@@@/\t/' \
	| cut -f 1,3 \
	| expand.pl \
	| cap.pl dummy \
	| row_stats.pl -mean \
	| tail +2 \
	> old_ave_set.tmp; \
	join.pl new_ave_set.tmp old_ave_set.tmp \
	> join.tmp; \
	cut -f 2,3 join.tmp \
	| awk '{print $$2/$$1*100;}' \
	| cap.pl dummy \
	| transpose.pl -q \
	| cap.pl dummy \
	| row_stats.pl -mean -std \
	| tail +2 \
	| cut -f 2- \
	| cap.pl RelativeSize_Mean,RelativeSize_Std \
	| transpose.pl -q \
	>> $@; \
	rm 1.tmp old.tmp new.tmp sizes.tmp old_ave_set.tmp new_ave_set.tmp join.tmp; \


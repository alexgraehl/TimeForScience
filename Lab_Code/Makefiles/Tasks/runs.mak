
include $(HOME)/Map/Makefile.common

RUN_DIRS = $(shell bind.pl $(RUN_META) | fasta2tab.pl | cut -f 1)
CD_DIRS  = $(addprefix cd@@@, $(RUN_DIRS))

all: 
	$(subst @@@, ,$(addsuffix ; make; cd ..;, $(CD_DIRS)))

clean:
	$(subst @@@, ,$(addsuffix ; make clean; cd ..;, $(CD_DIRS)))

# Eval under Map data hierarchy
ifneq (,$(EVAL_EXISTS))
eval_dir:
	$(MAKE_MAP_EVAL_DIR); \
	cd Eval; \
	$(RM) Makefile; \
	make_propagate.pl; \
	cd $(EVAL2RUN_DIR); \
	$(subst @@@, ,$(addsuffix ; make eval_dir; cd ..;, $(CD_DIRS)));
else
eval_dir:
	$(MAKE_EVAL_DIR); \
	$(subst @@@, ,$(addsuffix ; make eval_dir; cd ..;, $(CD_DIRS)));
endif

eval_dir_clean:
	cd Eval; \
	make maps_clean;

eval:
	$(subst @@@, ,$(addsuffix ; make eval; cd ..;, $(CD_DIRS)))

eval_clean:
	$(subst @@@, ,$(addsuffix ; make eval_clean; cd ..;, $(CD_DIRS)))

maps: $(RUN_META)
	bind.pl $(RUN_META) | create_maps.pl -

maps_clean:
	$(RMDIR) $(RUN_DIRS)


child_runs:
	$(addprefix cp $(MAP_TEMPLATES)/Runs/$(RUN_META) ,$(addsuffix ;, $(RUN_DIRS)))

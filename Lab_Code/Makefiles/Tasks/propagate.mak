
include $(HOME)/Map/Makefile.common

RUN_DIRS = $(shell /bin/ls -p | grep "/" | sed 's/\///' | transpose.pl)
CD_DIRS  = $(addprefix cd@@@, $(RUN_DIRS))

all: 
	$(subst @@@, ,$(addsuffix ; make; cd ..;, $(CD_DIRS)))

clean:
	$(subst @@@, ,$(addsuffix ; make clean; cd ..;, $(CD_DIRS)))

eval_dir:
	$(subst @@@, ,$(addsuffix ; make eval_dir; cd ..;, $(CD_DIRS)))

eval_dir_clean:
	$(subst @@@, ,$(addsuffix ; make eval_dir_clean; cd ..;, $(CD_DIRS)))

eval:
	$(subst @@@, ,$(addsuffix ; make eval; cd ..;, $(CD_DIRS)))

eval_clean:
	$(subst @@@, ,$(addsuffix ; make eval_clean; cd ..;, $(CD_DIRS)))

maps:
	$(subst @@@, ,$(addsuffix ; make maps; cd ..;, $(CD_DIRS)))

maps_clean:
	$(RMDIR) $(RUN_DIRS)

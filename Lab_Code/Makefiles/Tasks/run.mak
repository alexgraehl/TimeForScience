include $(HOME)/Map/Makefile.common

ifeq (,$(shell exists.pl after.exe))

all: time.txt

else

all: time.txt
	./after.exe

endif

ifeq (,$(shell exists.pl before.exe))

time.txt:
	date > time.txt; \
	setenv | grep "^HOST=" | cut -f2 -d "=" >> time.txt; \
	$(MAP_EXE) data.map; \
	date >> time.txt;

else

time.txt:
	before.exe; \
	date > time.txt; \
	$(MAP_EXE) data.map; \
	date >> time.txt;

endif

clean:
	$(RM) $(RUN_RESULTS) time.txt map.log;

ifneq (,$(EVAL_EXISTS))
eval_dir:
	$(MAKE_MAP_EVAL_DIR); \
	cd Eval; \
	$(RM) data.run; \
	$(LINK) $(dir $(EVAL2RUN_DIR))data.eval data.run; \
	make_runs.pl; \
	make maps;

else
eval_dir:
	$(MAKE_EVAL_DIR)
endif

eval_dir_clean:
	$(RMDIR) $(EVAL_DIR); \
	$(RM) Eval

eval:
	cd Eval; \
	make

eval_clean:
	cd Eval; \
	make clean

misc:
	../misc.pl

make:


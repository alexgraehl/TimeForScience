include $(HOME)/Map/Makefile.common

ifeq (,$(shell exists.pl after.exe))

all: data.time

else

all: data.time
	./after.exe

endif

ifeq (,$(shell exists.pl before.exe))

data.time:
	date > data.time; \
	setenv | grep "^HOST=" | cut -f2 -d "=" >> data.time; \
	$(MAP_EXE) data.map; \
	date >> data.time;

else

data.time:
	before.exe; \
	date > data.time; \
	$(MAP_EXE) data.map; \
	date >> data.time;

endif

clean:
	$(RM) $(RUN_RESULTS) data.time map.log;

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


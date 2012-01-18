include $(HOME)/Map/Makefile.common
include $(HOME)/Map/Data/Makefile.common
include $(HOME)/Map/Templates/Make/parent.mak

all:
	$(subst @@@, ,$(addsuffix ; make; cd ..;, $(CD_CHILD_DIRS)))

decorrelated:
	$(foreach org, $(ORGANISMS), \
	   cd $(org); \
	   make decorrelated; )


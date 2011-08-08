
#CHILDREN = Correlation CityBlock
CHILDREN = Correlation

include $(HOME)/Map/Templates/Make/quick.mak

make: Lib/child.mak
	$(foreach c, $(CHILDREN), \
	   mkdir -p $(c); \
	   cd $(c); \
	   rm -f Makefile; \
	   ln -s $(CHILD_MAKEFILE) Makefile; \
	   make make; \
	   cd ..; \
	) \

Lib/child.mak: Makefile
	mkdir -p Lib; \
	echo 'CHILDREN       = 10 20 30' \
	> $@; \
	echo 'CHILD_MAKEFILE = $(HOME)/Map/Templates/Make/kmeans.mak' \
	>> $@; \
	echo 'include $(HOME)/Map/Templates/Make/quick.mak' \
	>> $@; \


# DIRS = cd Fly/    \
#        cd Human/  \
#        cd Mouse/  \
#        cd Worm/   \
#        cd Yeast/  \
#

DIRS = cd Human/    \
       cd Mouse/  \
       cd Worm/   \

MAKE_DIRS = $(DIRS:/=; make_command & cd ..;)

TARGETS       = $(MAKE_DIRS:make_command=make)
CLEAN         = $(MAKE_DIRS:make_command=make clean)
REMOTE        = $(MAKE_DIRS:make_command=make remote)
REMOTE_CLEAN  = $(MAKE_DIRS:make_command=make remote_clean)

all: 
	$(TARGETS)

clean:
	$(CLEAN)

remote:
	$(REMOTE)

remote_clean:
	$(REMOTE_CLEAN)


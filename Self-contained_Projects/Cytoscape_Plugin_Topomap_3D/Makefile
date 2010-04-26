
#JAVA3D=/System/Library/Java/Extensions/vecmath.jar:/System/Library/Java/Extensions/j3dcore.jar:/System/Library/Java/Extensions/j3dutils.jar

#JAVALANG=/System/Library/Frameworks/JavaVM.framework/Versions/CurrentJDK/Classes/charsets.jar:/System/Library/Frameworks/JavaVM.framework/Versions/CurrentJDK/Classes/classes.jar:/System/Library/Frameworks/JavaVM.framework/Versions/CurrentJDK/Classes/dt.jar:/System/Library/Frameworks/JavaVM.framework/Versions/CurrentJDK/Classes/jce.jar:/System/Library/Frameworks/JavaVM.framework/Versions/CurrentJDK/Classes/jsse.jar:/System/Library/Frameworks/JavaVM.framework/Versions/CurrentJDK/Classes/laf.jar:/System/Library/Frameworks/JavaVM.framework/Versions/CurrentJDK/Classes/sunrsasign.jar:/System/Library/Frameworks/JavaVM.framework/Versions/CurrentJDK/Classes/ui.jar

# Two environment variables must be set:
#	CYTOSCAPE_DIR = location of cytoscape application (the folder where the application is located)
#	CYTOSCAPE_SRC = location of cytoscape source code (base folder--NOT the src folder, but the enclosing folder)

ifeq (cvaske, $(USER))
  CYTOSCAPE_DIR=/Applications/cytoscape-v2.1
  CYTOSCAPE_SRC=$(HOME)/src/cytoscapeSource-v2.1
endif

ifeq (alexgw, $(USER))
  CYTOSCAPE_DIR=$(HOME)/Alex_Lab/Cytoscape_Things/Cytoscape_from_CVS/cytoscape
  CYTOSCAPE_SRC=$(HOME)/Alex_Lab/Cytoscape_Things/Cytoscape_from_CVS/cytoscape
endif

CYTOSCAPE_PLUGINS=$(CYTOSCAPE_DIR)/plugins
CYTOSCAPE_JAR    =$(CYTOSCAPE_SRC)/cytoscape.jar
COLTGINY         =$(CYTOSCAPE_SRC)/lib/coltginy.jar
GINY             =$(CYTOSCAPE_SRC)/lib/giny.jar
MORE_JAR		 =$(CYTOSCAPE_SRC)/lib/phoebe.jar:$(CYTOSCAPE_SRC)/lib/cytoscape-graph-legacy.jar:$(CYTOSCAPE_SRC)/lib/concurrent.jar:$(CYTOSCAPE_SRC)/lib/cytoscape-task.jar
HUD_JAR          =./jars/newdawnhud.jar

CLASSES_DIR		= classes
APPNAME			= Topomap3DPlugin
MANIFEST_NAME	= Manifest.txt
CVS_TAR_ITEMS   = CVS_Tarred_Items
XCODE_TAR       = CVS_XCode_Project
# put the resulting classes in their own directory--don't clutter up the source folder
PUT_CLASSES_IN_SEPARATE_DIR = -d $(CLASSES_DIR)
# -source 1.4 is required to make "assert" work
SOURCE_VERSION              = -source 1.4

CFLAGS          = $(SOURCE_VERSION) $(PUT_CLASSES_IN_SEPARATE_DIR)
# -deprecation

export CLASSPATH=.:$(CYTOSCAPE_JAR):$(JAVALANG):$(JAVA3D):$(GINY):$(COLTGINY):$(MORE_JAR):$(HUD_JAR)

IMPORTANT_SOURCE_FILES =  *.java $(MANIFEST_NAME) $(CVS_TAR_ITEMS).tar $(XCODE_TAR).tar README.html Makefile

# vpath: notify "make" that classes are located in a different folder
vpath %.class $(CLASSES_DIR)

all: setdirs classes

setdirs:
	mkdir -p $(CLASSES_DIR)

# :=  <--- expand the right side only once
#  =  <--- expand the right side every time we see it

FILES = $(wildcard *.java)
classes: $(FILES:.java=.class)

#CLASSES = $(wildcard *.class)

%.class: %.java
	javac $(CFLAGS) $<

show_variables:
	echo Cytoscape dir: $(CYTOSCAPE_DIR)
	echo Cytoscape src: $(CYTOSCAPE_SRC)
	echo Classpath: $(CLASSPATH)

#MYVAR = $(shell find org -type f | grep -v CVS)
jar: all
	jar cmf $(MANIFEST_NAME) $(APPNAME).jar -C $(CLASSES_DIR) . -C lib .

	
#	jar i $(APPNAME).jar
#-C lib .

# Jar:
# -C dir
# Temporarily changes directories to dir while processing the following inputfiles argument. Multiple -C dir inputfiles sets are allowed.

look: # look at the jar file...
	jar tf $(APPNAME).jar

clean: bu
	rm -f $(APPNAME).jar $(CLASSES_DIR)/*.class

cyto: jar
	cp -f $(APPNAME).jar $(CYTOSCAPE_PLUGINS)/$(APPNAME).jar

runmac: cyto
	open $(CYTOSCAPE_DIR)/Cytoscape.app

# Xmx512M sets the maximum amount of memory that Java is allowed to allocate
run: cyto
	java -Xmx512M -jar $(CYTOSCAPE_DIR)/cytoscape.jar cytoscape.CyMain --JLD $(CYTOSCAPE_DIR)/plugins -i testdata.sif

rundebug: cyto
	java -Xdebug -Xnoagent -Djava.compiler=NONE -Xrunjdwp:transport=dt_shmem,address=jdbconn,server=y,suspend=n -Xmx512M -jar $(CYTOSCAPE_DIR)/cytoscape.jar cytoscape.CyMain --JLD $(CYTOSCAPE_DIR)/plugins -i testdata.sif

# Standalone: makes a double-clickable jar file (but won't work for cytoscape)	
standalone: all
	jar cmf $(MANIFEST_NAME) Topomap3D_standalone_app.jar -C $(CLASSES_DIR) .

runstandalone: standalone
	java -jar Topomap3D_standalone_app.jar

bu: # backup
	tar -cf backup.tar $(IMPORTANT_SOURCE_FILES)
	gzip backup.tar
	mkdir -p Backup
	mv -f backup.tar.gz Backup/

tarthings:
	tar -cf $(CVS_TAR_ITEMS).tar jars lib
	tar -cf $(XCODE_TAR).tar *.xcodeproj

update:
	cvs update
	tar -xf $(CVS_TAR_ITEMS).tar
	tar -xf $(XCODE_TAR).tar

commit: tarthings bu
	cvs commit -m "no message" $(IMPORTANT_SOURCE_FILES)

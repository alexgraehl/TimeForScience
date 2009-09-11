#!/bin/bash
#### This is a common init file for people in the lab

#### If you make any changes to this file, make SURE to launch
#### a new terminal and make sure you didn't get any error messages.
#### Any problems here will affect the entire lab!

# This makes pipelines return an error if any command exits with an error code
# By default BASH only returns the exit status of the rightmost command
# which does not play well with make
set -o pipefail
export GREP_OPTIONS=--color=auto

# "cut -f1 -d," removes the part after the comma on machines
# like SUNW,something (soe.ucsc.edu is an example). The comma messes up setenv.
export MACHTYPE=`uname -i | cut -f1 -d,`

export CVSROOT=/projects/compbio3/cvsroot/stuartlab
export SYSBIODIR=/projects/sysbio
export COMPBIODIR=/projects/compbio
export MAPDIR=$SYSBIODIR/map

export SYSBIOBINJAVA=$SYSBIODIR/apps/java
export SYSBIOBINPERL=$SYSBIODIR/apps/perl
export SYSBIOBINPYTHON=$SYSBIODIR/dev/python
export SYSBIOINMATLAB=$SYSBIODIR/apps/matlab
export SYSBIOMATLABPATH=$SYSBIODIR/dev/matlab

##
## If a JAVA_HOME is already specified, then let the user keep it
##
if [ -z ${JAVA_HOME} ]; then
    export JAVA_HOME=$SYSBIODIR/apps/$MACHTYPE/jdk/current
fi
export JAVABIN=${JAVA_HOME}/bin

## **********  Groovy  ***********
## The following options are good for performance of long-running programs. 
## export JAVA_OPTS="-Xms500m -Xmx800m -XX:CompileThreshold=0 -server"
## GROOVY_HOME is what groovy uses to find it's lib and other support files.
export SYSBIOGROOVY=$SYSBIODIR/apps/groovy
export GROOVY_HOME=$SYSBIODIR/apps/groovy

export SYSBIOLIB=$SYSBIODIR/apps/$MACHTYPE/lib
export SYSBIOPKGCONFIG=$SYSBIOLIB/pkgconfig
export SYSBIOBIN=$SYSBIODIR/apps/$MACHTYPE/bin
export XERCESCROOT=$SYSBIOLIB/xerces-c-src_2_7_0
export GRAPHVIZLIB=$SYSBIOLIB/graphviz

export SYSBIODEVPERL=$SYSBIODIR/dev/perl
export SYSBIOPERLLIB=$SYSBIODIR/lib/perl5
export SYSBIORLIB=$SYSBIODIR/apps/$MACHTYPE/Rlib
export MATLABPATH=$SYSBIOINMATLAB/MatlabSetupFiles:{$SYSBIOMATLABPATH}

#***  export MYPERLDIR=$SYSBIODEVPERLadded by Junior Josh to fix "Can't locate /lib/libfile.pl in @INC" from bash bug.
export MYPERLDIR=$SYSBIODEVPERL

export SYSBIO_LD_PATH=${SYSBIOLIB}:${XERCESCROOT}/lib:${GRAPHVIZLIB}:${GENIEDIR}/Lib/$MACHTYPE

# If we are x86_64, allow fallback to i386
case $MACHTYPE in
    x86_64)
    export SYSBIOLIB=$SYSBIOLIB:$SYSBIODIR/apps/i386/lib:$SYSBIODIR/apps/x86_64/mysql/lib
    export SYSBIOBIN=$SYSBIOBIN:$SYSBIODIR/apps/i386/bin:$SYSBIODIR/apps/x86_64/mysql/bin
    export SYSBIO_LD_PATH=${SYSBIO_LD_PATH}:${JAVA_HOME}/jre/lib/amd64
    ;;
    *)
    ;;
esac

if [ -z ${PERL5LIB} ]; then
   export PERL5LIB="${SYSBIOPERLLIB}"
else
   export PERL5LIB="${PERL5LIB}:${SYSBIOPERLLIB}"
fi

if [ -z ${GENIEDIR} ]; then
    export GENIEDIR=/cse/faculty/jstuart/src/stuartlab/genie_release
fi

## External environment variables needed for other programs

if [ -z ${LD_LIBRARY_PATH} ]; then
    export LD_LIBRARY_PATH=".:${SYSBIO_LD_PATH}"
else
    export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${SYSBIO_LD_PATH}"
fi

if [ -z ${PKG_CONFIG_PATH} ]; then
    export PKG_CONFIG_PATH="${SYSBIOPKGCONFIG}"
else
    export PKG_CONFIG_PATH="${SYSBIOPKGCONFIG}:${PKG_CONFIG_PATH}"
fi

export SYSBIOPATH=$SYSBIOBINPYTHON/Tools:$JAVABIN:$SYSBIOBIN:$SYSBIOBINJAVA:$SYSBIODEVPERL/Tools:$SYSBIODIR/dev/R_shell:$SYSBIODEVPERL/Map:$SYSBIODEVPERL/web:$SYSBIODIR/dev/alone/bin/$MACHTYPE:$COMPBIODIR/bin:$GENIEDIR/Programs/$MACHTYPE:$SYSBIOGROOVY/bin

## MMTx is a program written by NIH to map free text to UMLS terms. This
## environment variable is needed for MMTx to run correctly. It points to
## the home directory of where the MMTx java classes reside. See the
## documentation at:
##
## http://mmtx.nlm.nih.gov/Download/download.shtml
export MMTX_PATH="${SYSBIOBINJAVA}/MMTx"

if [ -z $PATH ] ; then
    export PATH=${SYSBIOPATH}
else
    export PATH=${SYSBIOPATH}:${PATH}
fi

if [ -z $R_LIBS ] ; then
    export R_LIBS=$SYSBIORLIB
else
    export R_LIBS=$SYSBIORLIB:$R_LIBS
fi

##
## Repair any broken permisions in the R library repository
##
/projects/sysbio/etc/fix_sysbio_r_libs_permissions.sh

if [[ -n "$PS1" && ("jig" = "$HOSTNAME") ]]; then
    echo '**************************************'
    echo '** Beware!                          **'
    echo '**                                  **'
    echo '** You have logged into JIG.        **'
    echo '**                                  **'
    echo '** Please note that this is the     **'
    echo '** FILE SERVER.                     **'
    echo '**                                  **'
    echo '** Running jobs on this machine     **'
    echo '** may slow down the sysbio         **'
    echo '** filesystem for everyone else!    **'
    echo '**                                  **'
    echo '** If you do not specifically want  **'
    echo '** to be using JIG, you should log  **'
    echo '** out and use another machine.     **'
    echo '**                                  **'
    echo '**************************************'
fi

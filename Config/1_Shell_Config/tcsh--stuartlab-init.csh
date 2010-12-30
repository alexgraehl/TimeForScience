#!/bin/tcsh

## This is a common init file for people in the lab
## Everyone should include it from their own ~/.cshrc file.
## with this line:   source /projects/sysbio/etc/stuartlab-init.csh

## WARNING ## WARNING ## WARNING ## WARNING ## WARNING
## WARNING ## WARNING ## WARNING ## WARNING ## WARNING
##
## DIRE WARNING FOR PEOPLE CHANGING THIS FILE:
##
## Changes to this file affect EVERYONE IN THE LAB.
## If there is an error in this file, it might prevent anyone else from
## logging in / running things.
##
## THEREFORE, if you change this file, ALWAYS LAUNCH A NEW TERMINAL
## and make sure you didn't get any error messages.
## Any problems here will affect the entire lab!
##
## Also, notify note that there is a stuartlab-init.sh file for bash users.
## Anything changed in this file WILL NOT APPLY TO BASH USERS. The bash
## init file will have to be modified separately.
## WARNING ## WARNING ## WARNING ## WARNING ## WARNING
## WARNING ## WARNING ## WARNING ## WARNING ## WARNING

# ADDITIONAL SEVERE WARNING FOR PEOPLE ADDING TO THIS FILE:
# If you put any commands in this file without checking to make sure $?prompt
# is true, then you will break the pager program "less"!!!
# If you need to print or do anything at all here, you MUST CHECK TO SEE
# if this is an interactive shell.
# DO NOT PRINT ANY OUTPUT TO NON-INTERACTIVE SHELLS!!!
# The solution here is that instead of just saying <do_something>
# You would say:    if ($?prompt && your_other_conditions) then <do_something> endif



if ($OSTYPE == darwin) then
   setenv MACHTYPE         `uname`
else
   # "cut -f1 -d," removes the part after the comma on machines
   # like SUNW,something (soe.ucsc.edu is an example). The comma messes up setenv.
   setenv MACHTYPE         `uname -i | cut -f1 -d,`
endif

setenv CVSROOT          /projects/compbio3/cvsroot/stuartlab
setenv SYSBIODIR        /projects/sysbio
setenv COMPBIODIR       /projects/compbio
setenv MAPDIR           $SYSBIODIR/map

# GREP_OPTIONS: --color=auto means "colorize the matching text when outputting to the terminal, but don't write the color characters when exporting to a file." You *don't* want to set this to --color=all, or it will mess up your redirects.
setenv GREP_OPTIONS     --color=auto

setenv SYSBIODEV        $SYSBIODIR/lab_apps
setenv SYSBIOBINJAVA    $SYSBIODIR/apps/java
## If a JAVA_HOME is already specified, then let the user keep it
if (! $?JAVA_HOME) then
    setenv JAVA_HOME $SYSBIODIR/apps/$MACHTYPE/jdk/current
endif
setenv JAVABIN ${JAVA_HOME}/bin

## Groovy 
setenv GROOVY_HOME    $SYSBIODIR/apps/groovy   

setenv SYSBIOBINPERL    $SYSBIODIR/apps/perl
setenv SYSBIOBINPYTHON  $SYSBIODIR/apps/python
setenv SYSBIOINMATLAB   $SYSBIODIR/apps/matlab
setenv SYSBIOMATLABPATH $SYSBIODIR/dev/matlab

setenv SYSBIOLIB        $SYSBIODIR/apps/$MACHTYPE/lib
setenv SYSBIOBIN        $SYSBIODIR/apps/$MACHTYPE/bin
setenv XERCESCROOT      $SYSBIOLIB/xerces-c-src_2_7_0
setenv GRAPHVIZLIB      $SYSBIOLIB/graphviz

setenv SYSBIODEVPERL    $SYSBIODEV/perl
setenv SYSBIOPERLLIB    $SYSBIODIR/lib/perl5

set SYSBIO_CPAN_INSTALL_LOCATIONS=${SYSBIOPERLLIB}/lib64/perl5/site_perl/5.8.5/x86_64-linux-thread-multi/:${SYSBIOPERLLIB}/lib/perl5/site_perl/5.8.5/

setenv OTHER_PEOPLE_PERL_APPS    $SYSBIODIR/apps/perl/bin

setenv SYSBIORLIB       $SYSBIODIR/apps/$MACHTYPE/Rlib
setenv MATLABPATH       $SYSBIOINMATLAB/MatlabSetupFiles:{$SYSBIOMATLABPATH}
setenv SYSBIOMYSQL      $SYSBIODIR/apps/x86_64/mysql
setenv SYSBIOMYSQLBIN   $SYSBIOMYSQL/bin


# Aliases (commands)
alias startmysql        $SYSBIODIR/etc/init.d/mysql.server start
alias stopmysql         $SYSBIODIR/etc/init.d/mysql.server stop

alias whichshell        echo $0   # Tells you which shell you used. Similar to "echo $SHELL"

alias ltab              sheet.py

# Tattle: who is using up all the CPU on this machine??
# Prints out a little box giving up the top "ps" output, sorted by CPU and formatted with sheet.pl.
alias tattle "echo -n '\033[45m\033[33m'; ps aux | tail -n +2 | sort --reverse -k 3,3 | head -n 5 | perl -p -e 's/[ ]+/\t/g' | cut -f 1,3,4,11 | cap.pl 'USER,CPU,MEM,TASK' | sheet.pl --ht=75 --trunc=60 | tail -n +2 ; echo -n '\033[0m'"

## MMTx is a program written by NIH to map free text to UMLS terms. This
## environment variable is needed for MMTx to run correctly. It points to
## the home directory of where the MMTx java classes reside. See the
## documentation at:
##
## http://mmtx.nlm.nih.gov/Download/download.shtml
setenv MMTX_PATH $SYSBIOBINJAVA/MMTx

# Perl 5 lib is where perl looks for @INC (include files)
# If you use CPAN to install a new module, it will go into
# SYSBIO_CPAN_INSTALL_LOCATION.
# Note that CPAN is invoked with <perl -MCPAN -e shell>
# You may have to configure CPAN separately, as the default
# install location might require you to have root access.
set ALL_PERL_LIB_LOCATIONS="${SYSBIOPERLLIB}:${SYSBIO_CPAN_INSTALL_LOCATIONS}"

if(! $?PERL5LIB) then
   setenv PERL5LIB "{$ALL_PERL_LIB_LOCATIONS}"
else
   setenv PERL5LIB "${PERL5LIB}:${ALL_PERL_LIB_LOCATIONS}"
endif

if(! $?MYPERLDIR) then
   setenv MYPERLDIR $SYSBIODEVPERL
endif

if(! $?GENIEDIR) then
   setenv GENIEDIR $SYSBIODEV/genie_release
endif

## External environment variables needed for other programs

setenv SYSBIO_LD_PATH "${XERCESCROOT}/lib:${GRAPHVIZLIB}:${GENIEDIR}/Lib/$MACHTYPE"

# Set up the path. This must happen BEFORE the switch statements below (which rely on SYSBIOPATH being set initially).

setenv SYSBIOPATH "$JAVABIN $GENIEDIR/Programs/$MACHTYPE $SYSBIOBIN $SYSBIOBINJAVA $SYSBIODEVPERL/Tools $SYSBIODEV/python/Tools $SYSBIODEV/R_shell $SYSBIODEV/shell_scripts $SYSBIODEVPERL/web $SYSBIODEVPERL/Map $SYSBIODEV/alone/bin/$MACHTYPE $MMTX_PATH/nls/mmtx/bin /cse/faculty/jstuart/src/stuartlab/alone/Tools/x86_64/ $COMPBIODIR/bin $SYSBIOMYSQLBIN $OTHER_PEOPLE_PERL_APPS $GROOVY_HOME/bin"

# If we are x86_64, allow fallback to i386
switch ($MACHTYPE)
	case x86_64:
		setenv SYSBIOLIB "$SYSBIOLIB $SYSBIODIR/apps/i386/lib $SYSBIODIR/apps/x86_64/mysql/lib"
		setenv SYSBIOPATH "$SYSBIOPATH $SYSBIODIR/apps/i386/bin $SYSBIODIR/apps/x86_64/mysql/bin"
		setenv SYSBIO_LD_PATH "${SYSBIO_LD_PATH}:${JAVA_HOME}/jre/lib/amd64"
	breaksw
	default: # did NOT recognize the machine architecture! (not x86_64, anyway)
		if ($?prompt) then
			echo 'Note: you have just logged into a non-x86_64 machine. Some programs that were compiled specifically for x86_64 may not work.\n'
		endif
	breaksw
endsw

# Note that the current directory (.) is included in the library path.
if(! $?LD_LIBRARY_PATH) then
   setenv LD_LIBRARY_PATH ".:${SYSBIO_LD_PATH}"
else
   setenv LD_LIBRARY_PATH "${LD_LIBRARY_PATH}:${SYSBIO_LD_PATH}"
endif

if(! $?path) then
   set path=($SYSBIOPATH)
else
   set path=($path $SYSBIOPATH)
endif

if(! $?R_LIBS) then
    setenv R_LIBS $SYSBIORLIB
else
    setenv R_LIBS "${SYSBIORLIB}:${R_LIBS}"
endif

if ($?prompt) then
# Sometimes we only want to do things for INTERACTIVE shells. Things
# that go in here should NOT break "less" and other programs.

# Attempt to set the user's group to sysbio. This will prevent the annoying
# situations in which an undergrad edits a file and then the group gets set
# to "ugrads" and no one else can edit it.
#newgrp sysbio

# note: the above command BREAKS EVERYTHING, so do not uncomment it
# until you can figure out why

endif


# Below: checks to see if this is an INTERACTIVE shell ($?prompt)
# and then if you are on a machine that is used for some specific non-compute-job purpose.
# If you ARE, then it displays a warning. Note that it pattern-matches so that both
# "jig.cse.ucsc.edu" and "jig" will both match.
if ($?prompt && ($HOST =~ 'jig*')) then
    echo -n "\033[41m\033[33m" # Set colors: 40m-47m are for background, 30m-37m are for foreground
    echo '**************************************'
    echo '** You have logged into JIG.        **'
    echo '**                                  **'
    echo '** Running jobs on this machine     **'
    echo '** may cause filesystem woes!       **'
    echo '**                                  **'
    echo '** This is the BACKUP server, with  **'
    echo '** files in /mnt/.snapsots .        **'
    echo '** Hustle contains a mirror, but    **'
    echo '** the mirror is NOT a backup.      **'
    echo '**                                  **'
    echo '** If you do not specifically want  **'
    echo '** to be using JIG, you should log  **'
    echo '** out and use another machine.     **'
    echo '**                                  **'
    echo '**************************************'
    echo -n "\033[0m" # Clear colors: restore colors back to clear
    echo 'The above message was generated when /projects/sysbio/etc/stuartlab-init.csh was loaded by your shell.\n'
endif

if ($?prompt && ($HOST =~ 'hustle*')) then
    echo -n "\033[45m\033[30m" # Set colors: 40m-47m are for background, 30m-37m are for foreground
    echo '**************************************'
    echo '** You have logged into HUSTLE.     **'
    echo '**                                  **'
    echo '** Please note that this was the    **'
    echo '** NEW FILE SERVER as of 2008.      **'
    echo '** It has a mirror of sysbio, but   **'
    echo '** this is a RAID mirror and not a  **'
    echo '** backup.                          **'
    echo '**                                  **'
    echo '** Running jobs on this machine     **'
    echo '** will slow down the sysbio        **'
    echo '** filesystem for everyone else!    **'
    echo '**                                  **'
    echo '** If you do not specifically want  **'
    echo '** to be using this machine, you    **'
    echo '** should use a different one.      **'
    echo '**                                  **'
    echo '**************************************'
    echo -n "\033[0m" # Clear colors: restore colors back to clear
    echo 'The above message was generated when /projects/sysbio/etc/stuartlab-init.csh was loaded by your shell.\n'
endif

if ($?prompt && ($HOST =~ 'disco*')) then
    echo -n "\033[42m\033[30m" # Set colors: 40m-47m are for background, 30m-37m are for foreground
    echo '**************************************'
    echo '** You have logged into DISCO.      **'
    echo '**                                  **'
    echo '** Please note that this is the     **'
    echo '** WEB SERVER as of 2008.           **'
    echo '**                                  **'
    echo '** Running jobs on this machine     **'
    echo '** or restarting it may impact      **'
    echo '** people using the web services!   **'
    echo '**                                  **'
    echo '** Also, some web-related programs  **'
    echo '** may *not* automatically come up  **'
    echo '** again after a reboot!            **'
    echo '**                                  **'
    echo '**************************************'
    echo -n "\033[0m" # Clear colors: restore colors back to clear
    echo 'The above message was generated when /projects/sysbio/etc/stuartlab-init.csh was loaded by your shell.\n'
endif

# Inform the user that this file finished loading. Only do this if the shell is interactive
# (otherwise LESS gets messed up for some reason)
if ($?prompt) then
    set history=32767   # Set the history length
    set autologout=0 # Don't auto logout! (Like on waterdance, in particular)
    set histdup=prev   # Don't put two duplicate commands in a row into the history
    set rmstar      # show a WARNING prompt before executing rm *
    set ignoreeof  # disable ctrl-d from logging out
    set color
    set colorcat # Note, if you have problems with colors when you run shell stuff from within your emacs shells, put this in your .emacs file: (add-hook 'shell-mode-hook 'ansi-color-for-comint-mode-on)
    set fignore=(\~ .o \#)   # fignore: suffixes to ignore for file completion
    /projects/sysbio/etc/fix_sysbio_r_libs_permissions.sh  # <-- run this script to fix R's permissions
    echo 'Finished loading config file </projects/sysbio/etc/stuartlab-init.csh>.'
endif

#alias make   make --warn-undefined-variables
#alias makeq  make --quiet --warn-undefined-variables





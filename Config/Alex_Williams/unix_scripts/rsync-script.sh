#!/bin/sh
# Author: Brice Burgess - bhb@iceburg.net
# rbackup.sh -- secure backup to a remote machine using rsync.

# Directories to backup. Separate with a space. Exclude trailing slash!
SOURCES="/Users/alexgw/Alex_Lab"
# can put more items in sources: "/a/b/c /b/c/d /d/e/f"

# IP or FQDN of Remote Machine
RMACHINE=roar.cse.ucsc.edu

# Remote username
RUSER=alexgw

# Location of passphraseless ssh keyfile
RKEY=/Users/alexgw/.rsync-key

# Directory to backup to on the remote machine. This is where your backup(s) will be stored
# Exclude trailing slash!
RTARGET="/projects/sysbio/alexgw/rsync-backups/from-onion"

# Your EXCLUDE_FILE tells rsync what NOT to backup. Leave it unchanged, missing or
# empty if you want to backup all files in your SOURCES. If performing a
# FULL SYSTEM BACKUP, ie. Your SOURCES is set to "/", you will need to make
# use of EXCLUDE_FILE. The file should contain directories and filenames, one per line.
# An example of a EXCLUDE_FILE would be:
# /proc/
# /tmp/
# /mnt/
# *.SOME_KIND_OF_FILE
#EXCLUDE_FILE="/path/to/your/exclude_file.txt"

# Comment out the following line to disable verbose output
VERBOSE="-v"

#######################################
########DO_NOT_EDIT_BELOW_THIS_POINT#########
#######################################

if [ ! -f $RKEY ]; then
  echo "Couldn't find ssh keyfile!"
  echo "Exiting..."
  exit 2
fi

if ! ssh -i $RKEY $RUSER@$RMACHINE "test -x $RTARGET"; then
  echo "Target directory on remote machine doesn't exist or bad permissions."
  echo "Exiting..."
  exit 2
fi

echo "Verifying Sources..."
for source in $SOURCES; do
        echo "Checking $source..."
        if [ ! -x $source ]; then
     echo "Error with $source!"
     echo "Directory either does not exist, or you do not have proper permissions."
     exit 2
   fi
done

if [ -f $EXCLUDE_FILE ]; then
EXCLUDE="--exclude-from=$EXCLUDE_FILE"
fi

echo "Sources verified. Running rsync..."
for source in $SOURCES; do

  # Create directories in $RTARGET to mimick source directory hiearchy
  if ! ssh -i $RKEY $RUSER@$RMACHINE "test -d $RTARGET/$source"; then
    ssh -i $RKEY $RUSER@$RMACHINE "mkdir -p $RTARGET/$source"
  fi

  rsync $VERBOSE $EXCLUDE -a --delete -e "ssh -i $RKEY" $source/ $RUSER@$RMACHINE:$RTARGET/$source/

done

exit 0
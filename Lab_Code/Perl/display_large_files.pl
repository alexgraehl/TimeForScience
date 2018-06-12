#!/usr/bin/env perl

# Script by Alex Williams

use strict;

my $FORCE_RUN = 0;

my $DOLLAR = '$'; # literal dollar sign
my $OUTPUT_FILE = "${DOLLAR}HOME/list_of_my_large_files.txt"; # note that the \$ is a LITERAL dollar sign, not a variable

my $TODAY = `date +%m-%d-%G`; # runs this command and gets the result of it
chomp($TODAY);


system("touch ${OUTPUT_FILE}"); # create the file if it doesn't already exist (doesn't matter if it already does)

my $LAST_RUN_DAY = `head -n 1 --silent $OUTPUT_FILE`;
chomp($LAST_RUN_DAY);


if (!defined($LAST_RUN_DAY) || (0==length($LAST_RUN_DAY))) { $LAST_RUN_DAY = "(none found)"; }

print "Today is: ${TODAY}\n";
print "Last run was on: ${LAST_RUN_DAY}\n";


if ($TODAY eq $LAST_RUN_DAY && (not $FORCE_RUN)) {
    print "Not running the space calculation script, because we already ran it today.\n";
    exit(0);
} else {
    
    my $DIRS = "${DOLLAR}HOME/* /projects/sysbio/`whoami`/* /projects/sysbioback/`whoami`/*";
    # Displays files with size more than 10,000 KB (~10 MB)

    my $OLD_FILE = "${OUTPUT_FILE}-old.tmp";
    
    system("mv -f ${OUTPUT_FILE} ${OLD_FILE}");

    system("echo $TODAY > $OUTPUT_FILE"); # add this to the top of the output file
    system("(du -ms ${DIRS} | sort -n | egrep '^[0-9]{2}' >> $OUTPUT_FILE) >& /dev/null"); # >& /dev/null");
    
    my $TEMP = "${OUTPUT_FILE}.tmp";

    system("cat ${OUTPUT_FILE} ${OLD_FILE} >> ${TEMP}");
    system("mv -f ${TEMP} ${OUTPUT_FILE}");

    system("rm -f $OLD_FILE $TEMP");
    
    # Note that this redirects stdout to the list file, and stderr (>&) to /dev/null.    
    
    
    
    
    
    
}




exit(0);




#!/bin/bash

# By Alex, 2015.
# This script checks the current time and decides whether or not to allow disk drives to spin down or not.
# Useful if you have a lot of external drives; otherwise OS X VERY SLOWLY (10+ seconds) waits for them to spin
#    up again whenever you trigger an open/save dialog.

# Basically allows the drive to sleep only when you are not likely to be using the computer.

# Note: must be run as root (use the "setuid" bit)

# Probably should add it to cron and run it every ***1st minute of every hour.***

# HOW TO USE:
# 1) First change the ownership of this script:
# SCRIPTNAME="this_script" &&  sudo chown root $SCRIPTNAME  &&  sudo chmod 744  $SCRIPTNAME  &&  sudo chmod u+s  $SCRIPTNAME

# 2) Then add it to the ROOT CRONTAB:
# env EDITOR=nano sudo crontab -e
# 1 <TAB> * <TAB> * <TAB> * <TAB> * <TAB> /full/path/to/this/script
# Make sure to use tabs, NOT spaces! There should be four asterisks after the "1".
# e.g.  1  *  *  *  *  (command path)

# I put the script in /Users/myname/bin/disksleep_set_by_time_agw.sh

NOW=$(date "+%u-%H-%M") # DAY_OF_WEEK:HOUR24:MINUTE60  (e.g 1:18:22 or 8:18:22)
# %u: day of week (1..7); 1 represents Monday, 6 is Saturday, 7 is Sunday

NOWARR=(${NOW//-/ }) # Split the result of the 'date' command by the hyphen
WEEKDAY=${NOWARR[0]} # NUMERIC: 1 = Monday... 6 = Saturday, 7 = Sunday
HR=${NOWARR[1]}  # 00 to 23
MIN=${NOWARR[2]} # 00 to 59

HR=$(( 10#$HR ))   # Prevent "08" and "09" from being interpreted as an invalid OCTAL number due to its leading zeroes
MIN=$(( 10#$MIN )) # Prevent "08" and "09" from being interpreted as an invalid OCTAL number due to its leading zeroes


if (("$WEEKDAY" >= "6")); then
    ISWEEKEND=1
else
    ISWEEKEND=0
fi

WEEKEND_START_SLEEP=3 # 3 AM. MUST BE AFTER MIDNIGHT OR THE COMPARISONS WON'T WORK.
WEEKEND_END_SLEEP=10  # 10 AM. MUST BE GREATER THAN "START"

WORKDAY_START_SLEEP=3 # 3 AM. MUST BE AFTER MIDNIGHT OR THE COMPARISONS WON'T WORK.
WORKDAY_END_SLEEP=17  # 5 PM. MUST BE GREATER THAN "START"

#echo $ISWEEKEND
#echo $WEEKDAY
#echo $HR
#echo $MIN

# pmset -c means only set it for when we are CONNECTED TO WALL POWER (-c). Does not affect disk spin when on battery.

if (( "$ISWEEKEND" == "1" && "$HR" >= "$WEEKEND_START_SLEEP" && "$HR" < "$WEEKEND_END_SLEEP" )); then
    # Ok, sleep the disk, it's a WEEKEND and it's between the "start" (inclusive) and "stop" (inclusive) sleep hours
    pmset -c disksleep 10
    echo "$NOW $ISWEEKEND: turning ON WEEKEND disk sleep at" `date` >> "/disksleep_status.txt"
    exit 0;
fi;

if (( "$ISWEEKEND" == "0" && "$HR" >= "$WORKDAY_START_SLEEP" && "$HR" < "$WORKDAY_END_SLEEP" )); then
    # Ok, sleep the disk, it's a WORKDAY and it's between the "start" (inclusive) and "stop" (inclusive) sleep hours
    pmset -c disksleep 10
    echo "$NOW $ISWEEKEND: turning ON WORKDAY disk sleep at" `date` >> "/disksleep_status.txt"
    exit 0; # our work here is done
fi;

echo "$NOW $ISWEEKEND: turning OFF disk sleep at" `date` >> "/disksleep_status.txt"

pmset -c disksleep 0 # I guess we are in the "do not sleep" period so turn OFF disk sleeping!


#!/usr/bin/env bash

# By Alex, 2015.
#
# Basically allows spinning (non-SSD) drives to sleep only when you are UNLIKELY to be using the computer.
# This script checks the current time and decides whether or not to allow disk drives to spin down or not.
# Useful if you have non-SSD external hard drives; otherwise OS X VERY SLOWLY (10+ seconds) waits for them
# to spin up again whenever you trigger an open/save dialog.
#
# *** Note: must be run as root (use the "setuid" bit) ***
#
# FYI: In the Mac OS 26 "Energy" control panel (as of August 2026), the UI-facing name of
#      this option is "Put hard tidks to sleep when possible"

# Recommendation: Probably should add it to cron and run it every ***1st minute of every hour.***

# HOW TO USE:
# 1) First change the ownership of this script:
#    SCRIPTPATH="this_script" &&  sudo chown root $SCRIPTPATH  &&  sudo chmod 744  $SCRIPTPATH  &&  sudo chmod u+s  $SCRIPTPATH
#
#    (I put the script in "/Users/myname/bin/disksleep_set_by_time_agw.sh", so that's what I would use for the SCRIPTPATH.

# 2) Then edit the ROOT CRONTAB to run this script periodically:
#     env EDITOR=nano sudo crontab -e
#     17<TAB>*<TAB>*<TAB>*<TAB>*<TAB>/full/path/to/this/script
#     Make sure to use tabs, NOT spaces! There should be four asterisks after the "17".
#     Visually, the result will look like this:
#     17   *    *    *    *    (command path)
#     Again, you must use TABS and not spaces.

NOW=$(date "+%u-%H-%M") # DAY_OF_WEEK:HOUR24:MINUTE60  (e.g 1:18:22 or 8:18:22)
# %u: day of week (1..7); 1 represents Monday, 6 is Saturday, 7 is Sunday

NOWARR=(${NOW//-/ }) # Split the result of the 'date' command by the hyphen
WEEKDAY=${NOWARR[0]} # NUMERIC: 1 = Monday... 6 = Saturday, 7 = Sunday
HR=${NOWARR[1]}  # 00 to 23
MIN=${NOWARR[2]} # 00 to 59

HR=${HR#0}   # strip leading 0--prevents "08" and "09" from being interpreted as an invalid OCTAL number due to its leading zeroes
MIN=${MIN#0} # strip leading 0--prevents "08" and "09" from being interpreted as an invalid OCTAL number due to its leading zeroes

if (("$WEEKDAY" >= "6")); then
    ISWEEKEND=1
else
    ISWEEKEND=0
fi

WEEKEND_START_SLEEP=3 # 3 AM. MUST BE AFTER MIDNIGHT OR THE COMPARISONS WON'T WORK.
WEEKEND_END_SLEEP=10  # 10 AM. MUST BE GREATER THAN "START"

WORKDAY_START_SLEEP=3 # 3 AM. MUST BE AFTER MIDNIGHT OR THE COMPARISONS WON'T WORK.
WORKDAY_END_SLEEP=17  # 5 PM. MUST BE GREATER THAN "START"

# pmset -c means only set it for when we are CONNECTED TO WALL POWER (-c). Does not affect disk spin when on battery.

if (( "$ISWEEKEND" == "1" && "$HR" >= "$WEEKEND_START_SLEEP" && "$HR" < "$WEEKEND_END_SLEEP" )); then
    # Ok, sleep the disk, it's a WEEKEND and it's between the "start" (inclusive) and "stop" (inclusive) sleep hours
    pmset -c disksleep 10
    #echo "$NOW $ISWEEKEND: turning ON WEEKEND disk sleep at" `date` >> "/disksleep_status.txt"
    exit 0;
fi;

if (( "$ISWEEKEND" == "0" && "$HR" >= "$WORKDAY_START_SLEEP" && "$HR" < "$WORKDAY_END_SLEEP" )); then
    # Ok, sleep the disk, it's a WORKDAY and it's between the "start" (inclusive) and "stop" (inclusive) sleep hours
    pmset -c disksleep 10
    #echo "$NOW $ISWEEKEND: turning ON WORKDAY disk sleep at" `date` >> "/disksleep_status.txt"
    exit 0; # our work here is done
fi;

#echo "$NOW $ISWEEKEND: turning OFF disk sleep at" `date` >> "/disksleep_status.txt"

pmset -c disksleep 0 # I guess we are in the "do not sleep" period so turn OFF disk sleeping!


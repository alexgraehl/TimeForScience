#!/bin/bash

# This just updates the sysbio directory to the latest CVS stuff.

# AFTER YOU COMMIT CHANGES TO CVS, even in your own home directory copy of the lab CVS
#   YOU SHOULD RUN THIS SCRIPT.

# Example:
# Suppose you were editing in your local copy of the CVS repository
#   ~/myLocalCvs/
# And you modified this file:
#   ~/myLocalCvs/science.pl
# Now you say:
#   cvs commit -m 'fixed some bugs'
# Now the CVS repository is updated.
# HOWEVER!
# The sysbio/lab_apps directory is NOT the master CVS repository--it is also a checked-out copy.
# Therefore, until someone runs "cvs update" in that directory (or in their own local copy of CVS),
# the only person who can use your newly-updated "science.pl" script is YOU.
# Also, some people don't have their own CVS checked-out copies of everything, so they won't have
# any idea that "science.pl" even exists until some individual runs cvs update in lab_apps.

# Therefore!

# You should run update_sysbio_dev_perl.sh ,
# and then everyone will be able to run the latest version of your new script.

# Note that there is still one issue:
# Some people may have their own CVS-repository-checked-out-copy of "science.pl" in their home
# directories.

# Unfortunately, they will still be using the old version of science.pl, since their home directory
# copy of science.pl probably comes earlier in the $PATH than the lab_apps copy.

# Regrettably, there is no way to inform people that a newer version is available, so you just have
# to have them run "cvs update" in their home directories.

cd /projects/sysbio/lab_apps/
cvs update -Pd


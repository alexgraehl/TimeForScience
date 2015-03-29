Q: What is in this repo?
* Enhanced unix tools (usually implemented in Perl). "cut.pl," for example, is like "cut," but allows you to easily re-order the output.
* Enhanced set-operation files. Joining / selecting specfic items / performing database operations on files is much easier than trying to just use "join." Also, "sets.pl" for turning a sets into lists of pairs, and vice-versa.
* Graphical / Plotting tools: programs like plot.R and various R scripts for making nice figures.
* UNIX configuration files. This is a good place to find usable emacs / shell / etc configurations


Q: Can people who are not in the lab read these files?

A: YES. This repository may be checked out by anyone, anywhere. It is not private at all. So don't put anything here that will cause your paper to get scooped, and try to avoid hilariously inappropriate code comments!

Q: What shell (environment) variables do I need to set to make use of this stuff?

A: Here are the BASH commands for setting things up:
export TIME_FOR_SCIENCE_DIR=/path/to/your/git/checkout/of/TimeForScience
export MYPERLDIR=${TIME_FOR_SCIENCE_DIR}/Lab_Code/Perl/

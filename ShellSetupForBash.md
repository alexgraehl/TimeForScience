## Setting up the $PATH for bash ##

`bash` is the default shell on Mac OS X. You can also run it at any time from the instructional machines by typing `bash` on the command line.

If you just check files out with Mercurial, some programs (Perl ones in particular) will not work properly, because various paths must be properly set.

Add these lines to the end of your `~/.bash_profile` . Put them at the very end, unless you have some reason not to.

```
## These lines make it possible to use files from the Mercurial source
## found at http://code.google.com/p/timeforscience/
export TIME_FOR_SCIENCE_DIR=$HOME/TimeForScience
export PATH=$TIME_FOR_SCIENCE_DIR/Lab_Code/Perl/Tools:$PATH
export PATH=$TIME_FOR_SCIENCE_DIR/Lab_Code/Perl/Scientific:$PATH
export PATH=$TIME_FOR_SCIENCE_DIR/Lab_Code/Python/Tools:$PATH
export PATH=$TIME_FOR_SCIENCE_DIR/Lab_Code/Shell:$PATH
export PATH=$TIME_FOR_SCIENCE_DIR/Lab_Code/R:$PATH
```

**Note that if your repository is named something OTHER than `TimeForScience`, you will need to change the `$HOME/TimeForScience` to `$HOME/WhateverYourNewNameIs`.**

**Do not change the `TIME_FOR_SCIENCE_DIR` variable name, no matter what your directory is named**
## Setting up the $path for tcsh and csh ##

`tcsh` is the default shell on the lab and instructional machines. You can always switch to `bash` by typing `bash` at the command prompt. `bash` is better in almost every way, but most users won't actually care about the (relatively minor) differences.

If you just check files out with Mercurial, some programs (Perl ones in particular) will not work properly, because various paths must be properly set.

You want to add the following lines to the very end of your `~/.cshrc` file, which is a file that is loaded every time a `tcsh` shell begins.

```
## These lines are useful for letting you use code from the lab Mercurial
## repository located at http://code.google.com/p/timeforscience/
if(! $?path) then
        set path=""
endif

setenv TIME_FOR_SCIENCE_DIR $HOME/TimeForScience

## Make sure all the lines below end with a "\", or else this command will not work!
## Do not put any comments on the same line after the backslash!
set path=($TIME_FOR_SCIENCE_DIR/Lab_Code/Perl/Tools \
        $TIME_FOR_SCIENCE_DIR/Lab_Code/Perl/Scientific \
        $TIME_FOR_SCIENCE_DIR/Lab_Code/Python/Tools \
        $TIME_FOR_SCIENCE_DIR/Lab_Code/Shell \
        $TIME_FOR_SCIENCE_DIR/Lab_Code/R \
        $path) ## <-- this line is the only one that does not have to end with a backslash

## End of mercurial-related lines
```


**Note that if your repository is named something OTHER than `TimeForScience`, you will need to change the `$HOME/TimeForScience` to `$HOME/WhateverYourNewNameIs`.**

**Do not change the `TIME_FOR_SCIENCE_DIR` variable name, no matter what your directory is named**
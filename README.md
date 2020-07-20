Q: What is in this repo?
* Enhanced unix tools (usually implemented in Perl or Python). 
  * `cut.pl` is like "cut," but allows you to easily re-order the output.
* Enhanced set-operation files. Joining / selecting specfic items / performing database operations on files is much easier than trying to just use "join." Also, `sets.pl` for turning a sets into lists of pairs, and vice-versa.
* Graphical / Plotting tools: programs like `plot.R` and various R scripts for making nice figures.
* UNIX configuration files. This is a good place to find usable emacs / shell / etc configurations. Do not put anything secret in here, or anything even seeminly innocuous like a server name!

Q: Can EVERYONE read these files?
A: YES. This repository may be checked out by anyone, anywhere. It is not private at all. So don't put anything here that will cause your paper to get scooped, and try to avoid hilariously inappropriate code comments!

Q: What shell (environment) variables do I need to set to make use of this stuff?

A: Here are the BASH commands for setting things up:
* `export TIME_FOR_SCIENCE_DIR=/path/to/your/git/checkout/of/TimeForScience`
* `export MYPERLDIR=${TIME_FOR_SCIENCE_DIR}/Lab_Code/Perl/`

## HOW TO CHECK OUT THIS REPO with SSH.
1. Get your hopefully-already-existing SSH key: cat ~/.ssh/id_rsa.pub
2. Add that key to Github (the "gear" icon in the top right, then "SSH Keys")
3. Test the connection in the terminal with ssh -T git@github.com
4. Now check out TimeForScience into your home directory:
  * `cd ~/ ; git clone https://github.com/alexgraehl/TimeForScience.git`
5. Finally, set your .bashrc to have an important TIME_FOR_SCIENCE_DIR environment variable, which some scripts rely on, by running the following TWO commands:
  * Command 1: `echo 'export TIME_FOR_SCIENCE_DIR="$HOME/TimeForScience"' >> ~/.bashrc`
  * Command 2: `echo 'export PATH="$PATH:$TIME_FOR_SCIENCE_DIR/Lab_Code/Perl/Tools:$TIME_FOR_SCIENCE_DIR/Lab_Code/Perl/Scientific:$TIME_FOR_SCIENCE_DIR/Lab_Code/Python/Tools:$TIME_FOR_SCIENCE_DIR/Lab_Code/Shell:$TIME_FOR_SCIENCE_DIR/Lab_Code/R"' >> ~/.bashrc ; source ~/.bashrc`


## How to set up your `$PATH` to use these tools

1. Find out what your shell is, by typing:
     echo $SHELL
   on the command line.

2. If the result of step 1 is "/bin/bash" or anything else ending in "bash", then you can just add these lines to the end of your ~/.bash_profile:

```export TIME_FOR_SCIENCE_DIR=$HOME/TimeForScience
export PATH=$TIME_FOR_SCIENCE_DIR/Lab_Code/Perl:$PATH
export PATH=$TIME_FOR_SCIENCE_DIR/Lab_Code/Python:$PATH
export PATH=$TIME_FOR_SCIENCE_DIR/Lab_Code/Shell:$PATH
export PATH=$TIME_FOR_SCIENCE_DIR/Lab_Code/R:$PATH
```

3. Check to see if this worked: open a new shell and type "which sets_overlap.pl" . It should say something like:
	/Users/yourname/TimeForScience/Lab_Code/Perl/sets_overlap.pl
	(This is what we want it to be)

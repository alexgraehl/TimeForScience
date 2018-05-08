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

HOW TO CHECK OUT THIS REPO with SSH.
0. You should use SSH and not HTTPS. With HTTPS, you have to type your password every time. Annoying!
1. Get your hopefully-already-existing SSH key: cat ~/.ssh/id_rsa.pub
2. Add that key to Github (the "gear" icon in the top right, then "SSH Keys")
3. Test the connection in the terminal with ssh -T git@github.com
4. To check out TimeForScience into your home directory:
  * `cd ~/ ; git clone git@github.com:alexgraehl/TimeForScience.git`
5. Finally, set your .bashrc to have an important TIME_FOR_SCIENCE_DIR environment variable, which some scripts rely on, by running the following TWO commands:
  * Command 1: `echo 'export TIME_FOR_SCIENCE_DIR="$HOME/TimeForScience"' >> ~/.bashrc`
  * Command 2: `echo 'export PATH="$PATH:$TIME_FOR_SCIENCE_DIR/Lab_Code/Perl/Tools:$TIME_FOR_SCIENCE_DIR/Lab_Code/Perl/Scientific:$TIME_FOR_SCIENCE_DIR/Lab_Code/Python/Tools:$TIME_FOR_SCIENCE_DIR/Lab_Code/Shell:$TIME_FOR_SCIENCE_DIR/Lab_Code/R"' >> ~/.bashrc ; source ~/.bashrc`

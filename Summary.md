# Time For Science: UNIX tools and other lab tools #

### Basic Information ###

Automatically generated information about this repository, including a link to how to find your google code password, can be found at: http://code.google.com/p/timeforscience/source/checkout (google code password can be found at http://code.google.com/hosting/settings )

### Q: What files are in this repository? ###
  1. Enhanced unix tools (usually implemented in Perl). "cut.pl," for example, is like "cut," but allows you to easily re-order the output.

  1. Enhanced set-operation files. Joining / selecting specfic items / performing database operations on files is much easier than trying to just use "join." Also, "sets.pl" for turning a sets into lists of pairs, and vice-versa.

  1. Graphical / Plotting tools: programs like plot.R and various R scripts for making nice figures.

  1. UNIX configuration files. This is a good place to find usable emacs / shell / etc configurations

### Q: Can people who are not in the lab read these files? ###

A: YES. This repository may be checked out by anyone, anywhere. It is not private at all. So don't put anything here that will cause your paper to get scooped, and try to avoid hilariously inappropriate code comments!


---


# How to use this code #

## First, you need Mercurial installed on your system! ##

If you are on the Mac or PC: http://mercurial.berkwood.com/

If you are on Linux, type  ` sudo apt-get install mercurial `. You will need administrator privileges to use this command.

## For the VERY FIRST TIME you want a Mercurial repository: ##

```
hg  clone  https://timeforscience.googlecode.com/hg/  ~/TimeForScience
```
This creates a new folder named `TimeForScience` in your home directory. You can see what's in it by typing `ls ~/TimeForScience`. Don't ever delete this folder unless you really know what you're doing!

You can move or rename the `TimeForScience` directory, but note that later instructions
assume that you did in fact name your repository `TimeForScience` and put it in your
home directory. If you did not, you will have to remember what you called it and change
the later instructions accordingly!


---

## Set up Terminal / Console Shell path ##
Now you will need to set up your shell's path, so that you can use the command line to run all the programs you just download.

There are two main shells you might be using: `bash` and `tcsh` (or `csh`).

1. To find out which shell you are using, open a terminal on the machine you are going to use and type
```
echo $SHELL
```

The result will probably either end in `bash` or end in `tcsh`.

2a. If it ends in `bash`, follow the instructions on the ShellSetupForBash page, and then come back for step 3.

2b. If it ends in `tcsh` or `csh`, follow the instructions on the ShellSetupForCsh page instead, and then come back for step 3.

3. Now open a **new** terminal window (you must make a new one before your new settings will take effect), and type
{{which sets\_overlap.pl`
}}

If the result is something like:
`/Users/yourname/*TimeForScience*/Lab_Code/Perl/Tools/sets_overlap.pl`

Then you are done.

If it is something like

If it says:
`/Users/yourname/cvs_directory/*perl/tools*/sets_overlap.pl`
Then something went wrong and nothing will work until you figure out what's wrong. Try opening a new terminal window again just in case. (This means that the $PATH is not actually set to use the updated Mercurial repository,
but is instead using the old lab CVS repository.)




---

## Updating your copy of the code ##

Occasionally you should update your code to get the latest bugfixes and new files. Here is the procedure:

First command:
```
cd ~/TimeForScience
```
> You must be in the source code repository for this to work! If you named your repository something other than `TimeForScience,` then you will obviously need to use that name instead.

Second command:
```
hg pull  &&  hg update
```
> `hg pull` checks the google code server to find out what has been changed.
> `hg update` actually overwrites the local files with the updated ones. This is non-destructive--your old files can still be obtained by "rolling back" the update (which is beyond the scope of this page).

> If you get the error message
```
abort: There is no Mercurial repository here (.hg not found)!
```
> then you are probably not in the `~/TimeForScience` directory or in one of its subdirectories. Type `pwd` to make sure you are in a child directory of `~/TimeForScience`.


---

## Uploading your own changes ##

Procedure for uploading your changes to the google code repository:

First command:
```
   cd ~/TimeForScience
```
> You must be in the source code repository for this to work!

Second command:
```
hg pull  &&  hg update
```
> `hg pull` checks the google code server to find out what has been changed.
> `hg update` actually overwrites the local files with the updated ones. This is non-destructive--your old files can

Third command (only required if you have completely new files to add):
```
hg add MY_NEW_FILE
```
> This tells Mercurial that you have a totally new file that you made. Make sure to put it in a sensible location!) You can skip this step if you aren't creating a totally new file.

Fourth and final command:
```
hg commit -m "Fixed the bug in this file, so I am uploading it"   &&  hg push
```
> `hg commit` tells your local computer that you have a change to make, but it does **not** send it to the server yet! (it is **not like** CVS commit!)

> `hg push` sends your updated file's changes to the server.

> Note: you will have to type your username and password here, UNLESS you put your password in your config file. Read on for instructions on how to do that.


---

## What is my password for committing changes? ##

Answer: it is an auto-generated password that can be found at: https://code.google.com/hosting/settings


---

## Solve Password Annoyance ##
To avoid having to type in the annoying auto-generated password every time, you can save your password in a local file. Note that anyone on the filesystem can read this file, so don't use a password you care about! The password is stored as plaintext! Anyone can read it!

There are two steps:

1. Here are the two important lines that should be in the file `~/TimeForScience/.hg/hgrc` :

```
[paths]
default = https://USERNAMEHERE:PASSWORDHERE@timeforscience.googlecode.com/hg/
```
> Obviously you will need to replace the `USERNAMEHERE:PASSWORDHERE` text with your username/password combination. Note that **your password will be stored in plain text in a file that anyone can read**, so it's probably **not safe** on the shared machines, although you can mitigate the risk with the `chmod` command in the next step.


2. For a tiny bit of additional security, now run:
```
chmod -R og-rwx ~/TimeForScience/.hg*
```
> to make it "very easy to read your password" instead of "incredibly easy."
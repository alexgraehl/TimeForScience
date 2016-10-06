# -*-Sh-*- <-- tells emacs what kind of syntax highlighting to use
#
[[ -z "$PS1" ]] && return # <-- If not running interactively, don't
# do anything. Printing any output breaks ssh and various things.
# This is why this line is important at the very top!

# .bashrc: Loaded when the shell is non-interactively started up

# ~/.bashrc: executed by bash(1) for non-login shells.
# see /usr/share/doc/bash/examples/startup-files (in the package bash-doc)
# for examples

# This site also has some useful commands: https://github.com/mrzool/bash-sensible/blob/master/sensible.bash

# ============= Set up things below ONLY if this shell is interactive ===============

function agw_cmd_exists() { # Check if a command exists
    type "$1" &> /dev/null
    # or try: if [[ -n `which exa 2> /dev/null` ]] ... # Usage example: if agw_cmd_exists "exa" && [ "$isMac" == "1" ] ; then ... 
}

if [[ -n "$SSH_CLIENT" || -n "$SSH2_CLIENT" ]] ; then is_sshing=1 ; fi   ## We are connected via SSH or SSH2... ## $SSH_CLIENT and $SSH2_CLIENT are set automatically by SSH.

# If we're in TMUX, then change the screen type to "screen-256color" and export the TERM

export TERM=xterm-256color # override the defaults and always assume 256 colors
case "$TERM" in
    xterm-color|xterm-256color|screen-256color)	color_prompt=1 ;;
    *)	                        ;;
esac
[[ -n "$TMUX" ]] && [[ "$color_prompt" == 1 ]] && export TERM=screen-256color
if [[ "$OSTYPE" == darwin* ]] ; then isMac="1" ; fi

COMPYNAME="$HOSTNAME" # <-- we will have to modify this if it's my home machine / some machine where $HOSTNAME doesn't work

command -v "scutil" > /dev/null # <-- Note: we check the exit code from this ("$?") below
if [[ "0" == $? ]] && [[ $(scutil --get ComputerName) == "Slithereens" ]]; then
    isAgwHomeMachine=1
    COMPYNAME="Slithereens"
fi

if [[ -n "$color_prompt" ]] ; then
    color_prefix="\033" ## <-- \033 works everywhere. \e works on Linux
    a_echo_color="${color_prefix}[1;32m" ## green  ## can also be: 1;33 ## was [3;40m before
    a_status_color="${color_prefix}[1;33m"
    a_warning_color="${color_prefix}[1;31m"
    a_end_color="${color_prefix}[m"
else
    color_prefix='' # Sadly, we do not have color on this terminal.
    a_echo_color=''
    a_status_color=''
    a_warning_color=''
    a_end_color=''
fi

echo -e "${a_echo_color}>>> BASH: Loading .bashrc...${a_end_color}" ## <-- comes after the colors are set up in platform-specific fashion

if [[ -d "$HOME/work" ]]; then
    export BINF_CORE_WORK_DIR="$HOME/work" # <-- set BINF_CORE work directory
elif [[ -d "/work" ]]; then
    export BINF_CORE_WORK_DIR="/work"  # <-- set BINF_CORE work directory
else
    echo "[WARNING in .bashrc: WORK directory not found]"
    export BINF_CORE_WORK_DIR="WORK_DIR_NOT_FOUND"
fi

if [[ -d "$HOME/TimeForScience" ]]; then ## If this is in the home directory, then set it no matter what.
    export TIME_FOR_SCIENCE_DIR="$HOME/TimeForScience"
elif [[ "$USER" == "alexgw" ]]; then ## Location of the TIME FOR SCIENCE directory. This is mostly Alex's code.
    export TIME_FOR_SCIENCE_DIR="$HOME/TimeForScience"
else
    export TIME_FOR_SCIENCE_DIR="/home/alexgw/TimeForScience"
    if [[ -d ${TIME_FOR_SCIENCE_DIR} ]]; then
	true # Don't do anything; it was found, which is fine
    else
	echo "[WARNING in .bashrc: TIME_FOR_SCIENCE_DIR not found]"
	export TIME_FOR_SCIENCE_DIR="COULD_NOT_FIND_TIME_FOR_SCIENCE_DIRECTORY_ON_FILESYSTEM"
    fi
fi

# Alias definitions.
# You may want to put all your additions into a separate file like
# ~/.bash_aliases, instead of adding them here directly.
# See /usr/share/doc/bash-doc/examples in the bash-doc package.
# Sets up the aliases whether or not this shell is interactive
if [[ -f ${TIME_FOR_SCIENCE_DIR}/Config/Alex_Williams/.aliases ]]; then
    source ${TIME_FOR_SCIENCE_DIR}/Config/Alex_Williams/.aliases
elif [[ -f ${BINF_CORE_WORK_DIR}/Code/TimeForScience/Config/Alex_Williams/.aliases ]] ; then
    source ${BINF_CORE_WORK_DIR}/Code/TimeForScience/Config/Alex_Williams/.aliases
elif [[ -f ${HOME}/TimeForScience/Config/Alex_Williams/.aliases ]] ; then
    source ${HOME}/TimeForScience/Config/Alex_Williams/.aliases
fi

function deduped() { # input: one string to de-dupe. Usage: PATH=$(deduped $PATH). May break on things with spaces.
    local X="$1"
    local new=""
    local old="$X:"
    if [ -n "$X" ]; then
	while [ -n "$old" ]; do
	    x=${old%%:*}
	    case "$new": in
		*:"$x":*) ;; # <-- entry EXISTS already
		*)new="$new:$x";; # Does not exist already
	    esac
	    old=${old#*:}
	done
	new=${new#:}
    fi
    new=$(echo "$new" | perl -pe 's/^[:+]//' | perl -pe 's/[:]+/:/g' | perl -pe 's/[:+]\$//') # remove leading/trailing ':'
    # Remove any leading or trailing ':' from the path
    echo "$new"
}

# ============================= PATH STUFF ============================
export BOWTIE_INDEXES=${BINF_CORE_WORK_DIR}/Apps/Bio/bowtie/current-bowtie/indexes/ ## <-- MUST have a trailing "/" after it!

export MYPERLDIR=${TIME_FOR_SCIENCE_DIR}/Lab_Code/Perl/
export R_BINF_CORE=${BINF_CORE_WORK_DIR}/Code/R

export PATH=/usr/bin:/bin:/usr/sbin:/sbin:/usr/X11/bin:/sw/bin ## <-- clear out initial path!!
export PATH=/opt/pbs/default/bin:/usr/local/bin:${BINF_CORE_WORK_DIR}/Apps/bin:/usr/local/sbin:/opt/bin:/projects/bin:${PATH}:${HOME}/.local/bin
export PATH=$(deduped "${HOME}/bin:${HOME}/.linuxbrew/bin:$TIME_FOR_SCIENCE_DIR/Lab_Code/Perl/Tools:$TIME_FOR_SCIENCE_DIR/Lab_Code/Perl/Scientific:$TIME_FOR_SCIENCE_DIR/Lab_Code/Python/Tools:$TIME_FOR_SCIENCE_DIR/Lab_Code/Shell:$TIME_FOR_SCIENCE_DIR/Lab_Code/R:${PATH}:${BINF_CORE_WORK_DIR}/Code/Python:${BINF_CORE_WORK_DIR}/Code/alexgw")   # <-- PRIORITY: programs in your own home directory come FIRST, then System-wide "Science" bin, then other stuff.

# PATH="$(echo $PATH | perl -e 'print join(":", grep { not $seen{$_}++ } split(/:/, scalar <>))')"  # remove dupes

BINFAPPBASE=/data/applications
BINFVERSION=2015_06
# 'export' added to "BINFSWROOT" on August 25, 2016
export BINFSWROOT=$(deduped $BINFAPPBASE/$BINFVERSION)
export BINFBINROOT=$(deduped $BINFSWROOT/bin)
export PERL5LIB=$(deduped $BINFSWROOT/libperl/lib/perl5:$HOME/.linuxbrew/lib/perl5/site_perl:$PERL5LIB)
export PERLLIB=$PERL5LIB
# PYTHONPATH: ------- UPDATE: Added August 25, 2016
export BINFPYROOT=$(deduped $BINFSWROOT/libpython2.7)
export PYTHONPATH=$(deduped $BINFPYROOT:$BINFPYROOT/lib64/python2.7/site-packages:$BINFPYROOT/lib/python2.7/site-packages:$PYTHONPATH)
#$BINFPYROOT:$BINFPYROOT/dist-packages:$BINFPYROOT/lib64/python2.7/site-packages:$PYTHONPATH)
export R_LIBS=$(deduped $BINFSWROOT/libr)
export PATH=$(deduped $BINFBINROOT:$PATH:$BINFPYROOT/bin)
export LD_LIBRARY_PATH=$(deduped "$HOME/.linuxbrew/lib:$BINFSWROOT/lib:$LD_LIBRARY_PATH")
export LIBRARY_PATH=$LD_LIBRARY_PATH
export CPATH=$(deduped $BINFSWROOT/include:$CPATH)

export CXX=/usr/bin/g++ # The location of the c++ compiler. NOT "cpp"--that is the C preprocessor.
export LANG="en_US.UTF-8"    # Set the LANG to UTF-8
#export LANG=C               # Set the LANG to C. Speeds some command line tools up, but *BREAKS* some UTF-8 stuff, like emacs!
#export LANG="en_US.UTF-8" &&  emacs --no-splash -nw ~/Downloads/XDownloads/unicode.txt

stty -ixon  # Disable the useless START/STOP output control (the one that pauses input is you press 'Ctrl-S' and resumes is you press 'Ctrl-Q')
shopt -s globstar # With globstar set (bash 4.0+), bash recurses all the directories. In other words, enables '**'
set   -o ignoreeof  # Prevent Ctrl-D from exiting! Still exits if you press it 10 times.
shopt -s checkwinsize # Check the window size after each command and, if necessary, update the values of LINES and COLUMNS.
shopt -s cdspell ## Fix obvious typeos in directory names when "cd"-ing
shopt -s cmdhist ## Save multi-line pasted commands into one single history command
#shopt -s lithist ##
shopt -s no_empty_cmd_completion ## Don't display ALL commands on an empty-line tab
shopt -s nocaseglob ## Match glob / regexp in case-insensitive fashion
set   -o noclobber    # Prevent file overwrite on stdout redirection. Override with ">|", e.g. echo 'a' >| file_that_exists
shopt -s histappend # Save terminal history between sessions
shopt -s cmdhist   # Save multi-line commands as one command

PROMPT_DIRTRIM=2   # Automatically trim long paths in the prompt (requires Bash 4.x)
bind 'set mark-directories on'           # show a '/' at the end of a directory name
bind 'set mark-symlinked-directories on'
bind "set completion-ignore-case on"     # Perform file completion in a case insensitive fashion
bind "set completion-map-case on"        # Treat hyphens and underscores as equivalent
bind "set show-all-if-ambiguous on"      # Display matches for ambiguous patterns at first tab press

# =========================== TERMINAL HISTORY =========================
export HISTSIZE=500000
export HISTFILESIZE=100000
export HISTCONTROL="erasedups:ignoreboth"  # Avoid duplicate entries
export HISTIGNORE="&:[ ]*:kpk:exit:p:pwd:rr:clear:history:fg:bg" ## Commands that are NOT saved to the history!
export HISTTIMEFORMAT='%F %T '
export PROMPT_COMMAND='history -a' ## save ALL terminal histories
# =========================== TERMINAL HISTORY =========================

# set variable identifying the chroot you work in (used in the prompt below)
# What the heck even is this
#if [[ -z $debian_chroot ]] && [[ -r /etc/debian_chroot ]]; then
#    debian_chroot=$(cat /etc/debian_chroot)
#fi

#PS1="[\D{%e}=\t \h:\W]$ "

highlight_text() { # prints colored text apparently. One argument, which is the color setting. 1 = red, 2 = green... etc
    if [ -x /usr/bin/tput ]; then tput bold; tput setaf $1; fi
    shift
    printf -- "$@"
    if [ -x /usr/bin/tput ]; then tput sgr0 ; fi # sgr0 = "turn off all attributes"
}

print_if_nonzero_exit_code() { # Example of using this in your prompt: PS1='$(highlight_exit_code)...' # <-- note the SINGLE quotes (or put a backslash before '$')
    exit_code=$?
    CTRLC_COLOR=2 # 2 = green
    ERR_COLOR=1 # 1 = red
    if [ $exit_code -eq 0 ]; then return;
    elif [ $exit_code -eq 130 ]; then highlight_text $CTRLC_COLOR "[Ctrl-C ($exit_code)]";
    else highlight_text $ERR_COLOR   "[Exit code $exit_code]"; fi
} 


export PS1="\$(print_if_nonzero_exit_code)\n[\D{%e}~\t~${COMPYNAME:0:3}~\W]$ " ## ${HOSTNAME:0:3} means only show the first 3 characters of the hostname! "\h" is the whole thing, also.
# Note: requires a "\n" after the "print_if_nonzero" or else Ctrl-A / Ctrl-E gets messed up when pasting

case "$COMPYNAME"
    in
    $RIGNODE_IP)
	export PS1="[COMPUTE_NODE] $PS1" # prepend "compute node" to it
	;;
    *) ## Otherwise...
	;;
esac

#if [[ -n "$color_prompt" ]] ; then
#    if [[ -z "$is_sshing" ]] ; then

## Conditional logic IS allowed in the ps1, bizarrely!
#PS1='$(if [[ $USER == alexgw ]]; then echo "$REGULAR_PROMPT"; else echo "$PWD%"; fi)'
#export PS1=\$(date +%k:%M:%S)

if [[ -n "$color_prompt" ]]; then
    #PS1="[\d \t \u@\h:\w]$ "
    #PS1="$(date +[%d]%H:%M) \u@\h:\W]\$ "
    #PS1="\[${a_machine_prompt_main_color}\]\$(date +[%d]%H:%M) \u@\h:\W]\$\[\033[0m\]\[${a_end_color} "
    #alias ls='ls --color=auto'
    #alias dir='dir --color=auto'
    #alias vdir='vdir --color=auto'
    alias grep='grep --color=auto'
    alias fgrep='fgrep --color=auto'
    alias egrep='egrep --color=auto'
else
    echo -n ''
    #PS1="$(date +[%d]%H:%M) \u@\h:\W]\$ "
fi

# enable color support of ls
[[ -x /usr/bin/dircolors ]] && eval "`dircolors -b`"

# Enable programmable completion features. May already be enabled in
# /etc/bash.bashrc or /etc/profile, in which case this would not be necessary.
if [[ -f /etc/bash_completion ]]; then
    . /etc/bash_completion
fi

export CLICOLOR='Yes'
export LS_OPTIONS='--color=auto'
export LS_COLORS='no=00:fi=00:di=01;35:ln=01;36:pi=40;33:so=01;35:do=01;35:bd=40;33;01:cd=40;33;01:or=43;31;03:ex=01;31'  ## LS_COLORS *with* an underscore is for Ubuntu
export LSCOLORS="FxgxCxDxBxegedabagacad"    ## LSCOLORS *without* an underscore is for Mac OS X.
export HOSTNAME  ## <-- Required for tmux / env to see HOSTNAME as a variable!
export    CVS_RSH=ssh
export  CVSEDITOR="emacs -nw"
export SVN_EDITOR="emacs -nw"
export     EDITOR="emacs -nw"

# Keybindings: Add IJKL navigation to supplement/replace the arrow keys
bind "\M-J:backward-word"
bind "\M-L:forward-word"
bind "\M-j:backward-char"
bind "\M-l:forward-char"
bind "\M-i:previous-history"
bind "\M-I:previous-history"
bind "\M-k:next-history"
bind "\M-K:next-history"
#bind "\C-r:history-search-backward"
bind "\C-s:history-search-forward" # This doesn't seem to work for some reason

umask u=rwx,g=rwx,o=rx # <-- give users and groups full access to files I create, and let other users READ and EXECUTE




# Save the local ethernet "en0" MAC address into the variable LOCAL_EN0_MAC. Note the zero.
# Allows per-machine settinsg.
#export LOCAL_EN0_MAC=`ifconfig en0 | grep -i ether | sed 's/.*ether //' | sed 's/[ ]*$//'`

## BETTER DIRECTORY NAVIGATION ##
# Prepend cd to directory names automatically
#shopt -s autocd
# Correct spelling errors during tab-completion
#shopt -s dirspell
# Correct spelling errors in arguments supplied to cd
#shopt -s cdspell

# This defines where cd looks for targets
# Add the directories you want to have fast access to, separated by colon
# Ex: CDPATH=".:~:~/projects" will look for targets in the current working directory, in home and in the ~/projects folder
#CDPATH="."

# This allows you to bookmark your favorite places across the file system
# Define a variable containing a path and you will be able to cd into it regardless of the directory you're in
#shopt -s cdable_vars

# Examples:
# export dotfiles="$HOME/dotfiles"
# export projects="$HOME/projects"
# export documents="$HOME/Documents"
# export dropbox="$HOME/Dropbox"







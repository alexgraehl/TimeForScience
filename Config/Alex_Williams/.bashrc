# -*-Sh-*- <-- tells emacs what kind of syntax highlighting to use
#

# https://www.shellcheck.net/

[[ -z "$PS1" ]] && return # <-- If not running interactively, don't
# do anything. Printing any output breaks ssh and various things.
# This is why this line is important at the very top!

# Source any global bash definitions
[[ -f "/etc/bashrc" ]] &&  . /etc/bashrc

# Enable programmable completion features. May already be enabled in  /etc/bash.bashrc or /etc/profile, in which case this would not be necessary.
[[ -f "/etc/bash_completion" ]] &&  . /etc/bash_completion

# ~/.bashrc: executed by bash(1) for non-login shells (Loaded when the shell is non-interactively started up)
# see /usr/share/doc/bash/examples/startup-files (in the package bash-doc) for examples
# This site also has some useful commands: https://github.com/mrzool/bash-sensible/blob/master/sensible.bash

# ~~~~~~~~~~~~ Set up things below ONLY if this shell is interactive ~~~~~~~

function agw_cmd_exists() { # Check if a command exists
    type "$1" &> /dev/null
    # Alternative approach: command -v "scutil" > /dev/null # <-- Note: we check the exit code from this ("$?") below
    # And then check the exit code:   HAS_PROGRAM=$((1-$?)) # $? = 0 means the above command SUCCEEDED
    # or try: if [[ -n `which exa 2> /dev/null` ]] ... # Usage example: if agw_cmd_exists "exa" && [ "$isMac" == "1" ] ; then ...
}


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

# If we're in TMUX, then change the screen type to "screen-256color" and export the TERM
export TERM=xterm-256color # override the defaults and always assume 256 colors
case "$TERM" in
    xterm-color|xterm-256color|screen-256color)	color_prompt=1 ;;
    *)	                        ;;
esac
[[ -n "${TMUX}" ]] && [[ "${color_prompt}" == 1 ]] && export TERM=screen-256color
if [[ "${OSTYPE}" == darwin* ]] ; then isMac="1" ; fi

# ~~~~~~~~~~~~~~~~~ See what the Mac (if it is one) thinks the computer's name is. This is the setting in the Sharing preference pane. ~~~~~~~~
if agw_cmd_exists "scutil"; then
    export MAC_SHARING_NAME=$(scutil --get ComputerName)
else
    export MAC_SHARING_NAME="PROBABLY_NOT_A_MAC_due_to_absence_of_the_SCUTIL_program"
fi

if [[ "${MAC_SHARING_NAME}" == "Slithereens" || "${MAC_SHARING_NAME}" == "Capsid" ]]; then
    export isAgwHomeMachine=1
    export COMPYNAME="Mac"
    export COMPYSHORT="Mac"
elif [[ "${MAC_SHARING_NAME}" =~ ^[Gg]host.* ]]; then # no quotation marks!
     export COMPYNAME="Ghost"
     export COMPYSHORT="Gho"
else
    export COMPYNAME="$HOSTNAME"
    export COMPYSHORT="$COMPYNAME"
fi

BLACK=$(tput setaf 0)
RED=$(tput setaf 1)
GREEN=$(tput setaf 2)
LIME_YELLOW=$(tput setaf 190)
YELLOW=$(tput setaf 3)
POWDER_BLUE=$(tput setaf 153)
BLUE=$(tput setaf 4)
MAGENTA=$(tput setaf 5)
CYAN=$(tput setaf 6)
WHITE=$(tput setaf 7)
BRIGHT=$(tput bold)
NORMAL=$(tput sgr0)
BLINK=$(tput blink)
REVERSE=$(tput smso)
UNDERLINE=$(tput smul)

if [[ -n "$color_prompt" ]] ; then
    cpre="\033" ## <-- Color prefix. \033 works everywhere. \e works on Linux
    a_echo_color="${cpre}[1;32m" ## green  ## can also be: 1;33 ## was [3;40m before

    # ######################
    # 256 color escapes from https://misc.flogisoft.com/bash/tip_colors_and_formatting
    a_green82="${cpre}[38;5;82m" # intense lime green
    a_orange136="${cpre}[38;5;136m"
    a_yellow226="${cpre}[38;5;226m" # bright yellow
    a_lavender93="${cpre}[38;5;93m"
    a_red="${cpre}[31m"
    a_green="${cpre}[32m"
    a_yellow="${cpre}[33m"
    a_blue="${cpre}[34m"
    a_magenta="${cpre}[35m"
    a_cyan="${cpre}[36m"
    a_lightgray="${cpre}[37m"
    a_darkgray="${cpre}[38m"

    # ######################
    a_red_bold="${cpre}[1;31m"
    a_green_bold="${cpre}[1;32m"
    a_yellow_bold="${cpre}[1;33m"
    a_blue_bold="${cpre}[1;34m"
    a_pink_bold="${cpre}[1;35m"
    a_cyan_bold="${cpre}[1;36m"
    a_white_bold="${cpre}[1;37m"
    a_status_color="${cpre}[1;33m"
    a_warning_color="${cpre}[1;31m"
    a_err_color=${a_red}
    a_ok_color=${a_green82}
    
    a_banner_color="${cpre}[1;45m"
    a_dirinfo_color="${cpre}${POWDER_BLUE}"
    a_end_color="${cpre}[m"
else
    cpre='' # Sadly, we do not have color on this terminal.
    a_echo_color=''
    a_status_color=''
    a_warning_color=''
    a_end_color=''
fi

echo -e "${a_echo_color}[:FYI:] Loading .bashrc...${a_end_color}" ## <-- comes after the colors are set up in platform-specific fashion

if [[ -d "${HOME}/work" ]]; then
    export BINF_CORE_WORK_DIR="${HOME}/work" # <-- set BINF_CORE work directory
elif [[ -d "/work" ]]; then
    export BINF_CORE_WORK_DIR="/work"  # <-- set BINF_CORE work directory
else
    export BINF_CORE_WORK_DIR="NO_BIOINFORMATICS_CORE_WORK_DIR_FOUND"
fi

if [[ -d "${HOME}/TimeForScience" ]]; then ## If this is in the home directory, then set it no matter what.
    export TIME_FOR_SCIENCE_DIR="${HOME}/TimeForScience"
else
    export TIME_FOR_SCIENCE_DIR="/home/alexgw/TimeForScience"
    if [[ -d ${TIME_FOR_SCIENCE_DIR} ]]; then
	true # Don't do anything; it was found, which is fine
    else
	echo "[:HEY:] WARNING in .bashrc: TIME_FOR_SCIENCE_DIR not found]"
	export TIME_FOR_SCIENCE_DIR="COULD_NOT_FIND_TIME_FOR_SCIENCE_DIRECTORY_ON_FILESYSTEM"
    fi
fi

# Alias definitions.
# You may want to put all your additions into a separate file like
# ~/.bash_aliases, instead of adding them here directly.
# See /usr/share/doc/bash-doc/examples in the bash-doc package.
# Sets up the aliases whether or not this shell is interactive
if [[ -f "${HOME}/.aliases" ]]; then
    source "${HOME}/.aliases"
elif [[ -f ${TIME_FOR_SCIENCE_DIR}/Config/Alex_Williams/.aliases ]]; then
    source ${TIME_FOR_SCIENCE_DIR}/Config/Alex_Williams/.aliases
elif [[ -f ${BINF_CORE_WORK_DIR}/Code/TimeForScience/Config/Alex_Williams/.aliases ]] ; then
    source ${BINF_CORE_WORK_DIR}/Code/TimeForScience/Config/Alex_Williams/.aliases
fi


# ~~~~~~~~~~~~~~~~~~ PATH ~~~~~~~~~~~~~~~~~~
export MYPERLDIR=${TIME_FOR_SCIENCE_DIR}/Lab_Code/Perl/

MINICONDA_BIN="${HOME}/miniconda3/bin"
export PATH="/usr/bin:/bin:/usr/sbin:/sbin:/usr/X11/bin:/sw/bin" ## <-- clear out initial path!!
export PATH="/opt/pbs/default/bin:/usr/local/bin:${BINF_CORE_WORK_DIR}/Apps/bin:/usr/local/sbin:/opt/bin:/projects/bin:${HOME}/.local/bin:${PATH}"
export PATH=$(deduped "${HOME}/bin:${HOME}/.local/bin:${HOME}/.linuxbrew/bin:${TIME_FOR_SCIENCE_DIR}/Lab_Code/Perl:${TIME_FOR_SCIENCE_DIR}/Lab_Code/Python:${TIME_FOR_SCIENCE_DIR}/Lab_Code/Python2:${TIME_FOR_SCIENCE_DIR}/Lab_Code/Python3:${TIME_FOR_SCIENCE_DIR}/Lab_Code/Shell:${TIME_FOR_SCIENCE_DIR}/Lab_Code/R:${PATH}:${BINF_CORE_WORK_DIR}/Code/Python:${MINICONDA_BIN}:${HOME}/workspace/0_CODE/programs:/Applications/Postgres.app/Contents/Versions/latest/bin")   # <-- PRIORITY: programs in your own home directory come FIRST, then System-wide "Science" bin, then other stuff.
export BINFSWROOT="" #$(deduped $BINFAPPBASE/$BINFVERSION)
export PERL5LIB=$(deduped $BINFSWROOT/libperl/lib/perl5:${HOME}/.linuxbrew/lib/perl5/site_perl:$PERL5LIB)
export PERLLIB=$PERL5LIB
#export R_LIBS=$(deduped $BINFSWROOT/libr)
export PATH=$(deduped $PATH) #:$BINFPYROOT/bin)
export LD_LIBRARY_PATH=$(deduped "${HOME}/.linuxbrew/lib:$LD_LIBRARY_PATH")
export LIBRARY_PATH=$LD_LIBRARY_PATH
#export CPATH=$(deduped $BINFSWROOT/include:$CPATH)
export CXX=/usr/bin/g++ # The location of the c++ compiler. NOT "cpp"--that is the C preprocessor.
export LANG="en_US.UTF-8"    # Set the LANG to UTF-8
stty -ixon  # Disable the useless START/STOP output control (the one that pauses input is you press 'Ctrl-S' and resumes is you press 'Ctrl-Q')
set   -o pipefail # Report a problem that occurs ANYWHERE in a piped command.
shopt -s globstar # With globstar set (bash 4.0+), bash recurses all the directories. In other words, enables '**'
set   -o ignoreeof  # Prevent Ctrl-D from exiting! Still exits if you press it 10 times.
shopt -s checkwinsize # Check the window size after each command and, if necessary, update the values of LINES and COLUMNS.
shopt -s cdspell ## Fix obvious typeos in directory names when "cd"-ing
#shopt -s lithist ##
#shopt -s no_empty_cmd_completion ## Don't display ALL commands on an empty-line tab
shopt -s nocaseglob ## Match glob / regexp in case-insensitive fashion
set   -o noclobber  # Prevent file overwrite on stdout redirection. Override with ">|", e.g. echo 'a' >| file_that_exists

export PROMPT_DIRTRIM=2   # Automatically trim long paths in the prompt (requires Bash 4.x)
bind 'set mark-directories on'           # show a '/' at the end of a directory name
bind 'set mark-symlinked-directories on'
bind "set completion-ignore-case on"     # Perform file completion in a case insensitive fashion
bind "set completion-map-case on"        # Treat hyphens and underscores as equivalent
bind "set show-all-if-ambiguous on"      # Display matches for ambiguous patterns at first tab press

# ~~~~~~~~~~~~~~~~~~~~ TERMINAL HISTORY ~~~~~~~~~~~~~~~~~~~~~
export HISTSIZE=500000 # num lines WHEN READING A NEW SESSION
export HISTFILESIZE=500000 # num lines ON DISK
export HISTCONTROL="ignoredups" # do not put 'erasedups' in unless you like annoyance
#export HISTIGNORE="&:[ ]*:kpk:exit:p:pwd:rr:clear:history:fg:bg" ## Commands that are NOT saved to the history!
export HISTTIMEFORMAT='%F %T '
#export PROMPT_COMMAND='history -a' ## save ALL terminal histories
shopt -s histappend # Save terminal history between sessions
shopt -s cmdhist ## Save multi-line pasted commands into one single history command
# Save and reload the history after each command finishes
export PROMPT_COMMAND="history -a; history -c; history -r; $PROMPT_COMMAND"


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

#parse_git_branch() { # Shows the current 'git' branch if you're in a git-controlled directory
#     git branch 2> /dev/null | sed -e '/^[^*]/d' -e 's/* \(.*\)/[\1]/'
#}
#RIGNODE_HOSTNAME=${RIGNODE_HOSTNAME:-no_rignode_hostname} # bash, assign a DEFAULT value if the argument is not defined
# ~~~~~~~~~~ Start of code block from http://ezprompt.net/ ~~~~~~~~~
function nonzero_return() {
    RETVAL=$?
    #echo "hello $? $RETVAL ${?#0}"
    [ $RETVAL -ne 0 ] && echo "$RETVAL"
}

# get current branch in git repo
function parse_git_branch() {
	BRANCH=`git branch 2> /dev/null | sed -e '/^[^*]/d' -e 's/* \(.*\)/\1/'`
	if [ ! "${BRANCH}" == "" ]
	then
		STAT=`parse_git_dirty`
		echo "[${BRANCH}${STAT}]"
	else
		echo ""
	fi
}

# get current status of git repo
function parse_git_dirty {
	status=`git status 2>&1 | tee`
	dirty=`echo -n "${status}" 2> /dev/null | grep "modified:" &> /dev/null; echo "$?"`
	untracked=`echo -n "${status}" 2> /dev/null | grep "Untracked files" &> /dev/null; echo "$?"`
	ahead=`echo -n "${status}" 2> /dev/null | grep "Your branch is ahead of" &> /dev/null; echo "$?"`
	newfile=`echo -n "${status}" 2> /dev/null | grep "new file:" &> /dev/null; echo "$?"`
	renamed=`echo -n "${status}" 2> /dev/null | grep "renamed:" &> /dev/null; echo "$?"`
	deleted=`echo -n "${status}" 2> /dev/null | grep "deleted:" &> /dev/null; echo "$?"`
	bits=''
	if [ "${renamed}" == "0" ]; then
		bits=">${bits}"
	fi
	if [ "${ahead}" == "0" ]; then
		bits="*${bits}"
	fi
	if [ "${newfile}" == "0" ]; then
		bits="+${bits}"
	fi
	if [ "${untracked}" == "0" ]; then
		bits="?${bits}"
	fi
	if [ "${deleted}" == "0" ]; then
		bits="x${bits}"
	fi
	if [ "${dirty}" == "0" ]; then
		bits="!${bits}"
	fi
	if [ ! "${bits}" == "" ]; then
		echo " ${bits}"
	else
		echo ""
	fi
}


export LOGLEVEL=DEBUG # python logging level

# \w and PROMPT_DIRTRIM=2 will give you a decent 'sampler' of the path
#export PS1="[\[\033[0;33m\]\$(date +%d%b) \$(date +%T) \[\033[0;34m\]\u@\h:] \[\033[1;31m\]\w\[\033[0m\]>"
#export PS1="\033[0;33m\]\$(date +%d%b) \$(date +%T) \[\033[0;34m\]\u@\h: \[\033[1;31m\]\w\[\033[0m\]$"

#PS1="[\D{%e}=\t \h:\W]$ "
# Note: all color sequence escapes must be surrounded by \[ ... \]. So all instances of \e must have \[ before them, or line wrap will break. See this: https://stackoverflow.com/questions/342093/ps1-line-wrapping-with-colours-problem . More specifically: "terminal escapes should be surrounded by \[ \] in order for bash to be able to correctly count printing characters to the beginning of the line. "
#PRE_PS1="\[${BRIGHT}\]\[${BLUE}\]\t\[${MAGENTA}\]\$(parse_git_branch)\[${NORMAL}\]"
#export PS1="\$(print_if_nonzero_exit_code)\n[\D{%e}~\t~${COMPYNAME:0:3}~\W]$ " ## ${HOSTNAME:0:3} means only show the first 3 characters of the hostname! "\h" is the whole thing, also.
# Note: requires a "\n" after the "print_if_nonzero" or else Ctrl-A / Ctrl-E gets messed up when pasting

a_red_bg_ps1="\[\e[33;41m\]"
a_yellow_bg_ps1="\[\e[31;43m\]"
a_blue_bg_ps1="\[\e[36;44m\]"

COMPCOL="\[\e[37;44m\]"

# For whatever reason, the end color for the PS1 must be THIS form below, or else previous/next line does not work properly
ENDCOL_PS1="\[\e[m\]"
#ENDCOL_ECHO="\e[0m"

test -e "${HOME}/.iterm2_shell_integration.bash" && source "${HOME}/.iterm2_shell_integration.bash" # loads 'it2setcolor' and others

# Color terminal settings / iTerm background color settings / iterm color settings
if [[ "${AFRESH__CLOUD_PROVIDER}" == 'azure' ]]; then
   if [[ "${AFRESH__ENV}" == 'development' ]]; then
       ENVTEXT="${ENDCOL_PS1} [DEV AZ]"
       if agw_cmd_exists "it2setcolor"; then
	   it2setcolor bg 000000  # Set iterm2 background color (terminal background color)
	   it2setcolor fg FFFFFF  # Set iterm2 foreground color (terminal foreground color)
       fi
   elif [[ "${AFRESH__ENV}" == 'staging' ]]; then
       ENVTEXT="${ENDCOL_PS1} [++STAGE AZ++]"
       if agw_cmd_exists "it2setcolor"; then
	   it2setcolor bg 442200  # Set iterm2 background/foreground colors
	   it2setcolor fg FFFFFF  # white
       fi
   elif [[ "${AFRESH__ENV}" == 'production' ]]; then
       ENVTEXT="${ENDCOL_PS1} [**PROD AZ***]"
       if agw_cmd_exists "it2setcolor"; then
	   it2setcolor bg 660000  # Set iterm2 background/foreground colors
	   it2setcolor fg FFBBFF  # white
       fi
   else
       ENVTEXT=""
   fi
   
else
    if [[ "${AFRESH__ENV}" == 'development' ]]; then
       ENVTEXT="${ENDCOL_PS1} [DEV AWS]"
       if agw_cmd_exists "it2setcolor"; then
	   it2setcolor bg 000033  # Set iterm2 background/foreground colors
	   it2setcolor fg FFFFFF
       fi
   elif [[ "${AFRESH__ENV}" == 'staging' ]]; then
       ENVTEXT="${ENDCOL_PS1} [++STAGE AWS++]"
       if agw_cmd_exists "it2setcolor"; then
	   it2setcolor bg 664400  # Set iterm2 background/foreground colors
	   it2setcolor fg FFFFFF  # white
       fi
   elif [[ "${AFRESH__ENV}" == 'production' ]]; then
       ENVTEXT="${ENDCOL_PS1} [**PROD AWS***]"
       if agw_cmd_exists "it2setcolor"; then
	   it2setcolor bg 660066  # Set iterm2 background/foreground colors
	   it2setcolor fg FFBBFF  # pink
       fi
   else
       ENVTEXT=""
   fi
fi
   
case "$COMPYNAME"
in
    #zzz*|*GYFH|Mac*|Slithereens|Capsid)
	# 37;41m is white on a GREEN (42) background. Red would be ;41m
	#export PS1="\[\e[37;42m\]\t\[\e[m\]\[\e[37;41m\]\`parse_git_branch\`\[\e[m\]\[\e[42m\][\[\e[m\]\[\e[37;42m\]${COMPYSHORT}\[\e[m\]\[\e[42m\]]\[\e[m\]\[\e[33;45m\]\`nonzero_return\`\[\e[m\] " # don't use \h for host here; we use COMPYNAME instead, since it was sanitized / redone
	#export PS1="${PRE_PS1}\[${POWDER_BLUE}\].mac\$\[${NORMAL}\] " # we are on a mac laptop probably
#	;;
    *) ## Otherwise...
	# 37;44m is white on blue
	#export PS1="\[\e[37;44m\]\t\[\e[m\]\[\e[37;41m\]\`parse_git_branch\`\[\e[m\]\[\e[44m\][\[\e[m\]\[\e[37;44m\]${COMPYSHORT}\[\e[m\]\[\e[44m\]]\[\e[m\]\[\e[33;45m\]\`nonzero_return\`\[\e[m\] " # don't use \h for host here; we use COMPYNAME instead, since it was sanitized / redone
	export PS1="\[\e[37;44m\]\t\[\e[m\]\[\e[37;41m\]\
\`parse_git_branch\`\
\[\e[m\]\[\e[44m\][\[\e[m\]\[${COMPCOL}${COMPYSHORT}${ENDCOL_PS1}\
\[\e[44m\]]\
${ENVTEXT}\
\[\e[m\]\[\e[33;45m\]\
\`nonzero_return\`\[\e[m\] " # don't use \h for host here; we use COMPYNAME instead, since it was sanitized / redone
	;;
esac

# ~~~~~~~~~~~~~~~~ End of code block from http://ezprompt.net/ ~~~~~~~~~~~~~~~

#if [[ -n "$SSH_CLIENT" || -n "$SSH2_CLIENT" ]] ; then is_sshing=1 ; fi   ## We are connected via SSH or SSH2... ## $SSH_CLIENT and $SSH2_CLIENT are set automatically by SSH.
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
    true
fi

# enable color support of ls
[[ -x "/usr/bin/dircolors" ]] && eval "`dircolors -b`"

export CLICOLOR='Yes'
export LS_OPTIONS='--color=auto'
export LS_COLORS='no=00:fi=00:di=01;35:ln=01;36:pi=40;33:so=01;35:do=01;35:bd=40;33;01:cd=40;33;01:or=43;31;03:ex=01;31'  ## LS_COLORS *with* an underscore is for Ubuntu
export LSCOLORS="FxgxCxDxBxegedabagacad"    ## LSCOLORS *without* an underscore is for Mac OS X.
export HOSTNAME  ## <-- Required for tmux / env to see HOSTNAME as a variable!
export    CVS_RSH="ssh"
export  CVSEDITOR="emacs"
export SVN_EDITOR="emacs"
export     EDITOR="emacs"

export HUE_BASE_IP="10.0.0.9"
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
# export dotfiles="${HOME}/dotfiles"
# export projects="${HOME}/projects"
# export documents="${HOME}/Documents"
# export dropbox="${HOME}/Dropbox"

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# From here: https://gist.github.com/pablete/5871811
# Changing iTerm2 color in MacOSX when SSHing (so you know at a glance that you're no longer in Kansas)
# Adapted from https://gist.github.com/porras/5856906
# 1. Create a theme in your terminal setting with the name "SSH" and the desired colors, background, etc.
# 2. Add this to your .bash_profile (or .bashrc, I always forget the difference ;))
# 3. Optional but useful: in the terminal, go to Settings > Startup and set "New tabs open with" to
#    "default settings" (otherwise, if you open a new tab from the changed one, you get a local tab with
#    the SSH colors)
function tabc() {
  NAME=$1; if [ -z "$NAME" ]; then NAME="Default"; fi # if you have trouble with this, change
                                                      # "Default" to the name of your default theme
  echo -e "\033]50;SetProfile=$NAME\a"
}

function colorssh() {
  tabc SSH
  ssh "$*"
  tabc
}

#alias ssh="colorssh"
# This would be easy to extend to check if a theme with the name of the server exists and set it, and
# fall back to the SSH theme. This way you can have different colors for different remote environments
# (per project, production, staging, etc.)

if [[ ! "${isMac}" == "1" ]]; then
    # Only do this if we are NOT on a Mac!
    if __bp_interactive_mode; then # <-- only run __bp_interactive_mode order to check that its exit status is ZERO
	echo -n ''
	# no problems
    else
	echo "[:HEY:] Unsetting the PROMPT_COMMAND, because this machine doesn't appear to support it. See the ~/.bashrc for details."
	unset PROMPT_COMMAND
    fi
fi

if [[ -f "${HOME}/bin/activate" ]]; then
    echo "[:OK:] .bashrc reporting: Activating virtualenv at ~/bin/activate"
    source "${HOME}/bin/activate"
fi

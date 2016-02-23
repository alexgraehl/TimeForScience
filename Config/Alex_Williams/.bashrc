# -*-Sh-*- <-- tells emacs what kind of syntax highlighting to use
#
[[ -z "$PS1" ]] && return # <-- If not running interactively, don't
# do anything. Printing any output breaks ssh and various things.
# This is why this line is important at the very top!

# .bashrc: Loaded when the shell is non-interactively started up

# ~/.bashrc: executed by bash(1) for non-login shells.
# see /usr/share/doc/bash/examples/startup-files (in the package bash-doc)
# for examples

# ============= Set up things below ONLY if this shell is interactive ===============

if [[ -n "$SSH_CLIENT" || -n "$SSH2_CLIENT" ]] ; then is_sshing=1 ; fi   ## We are connected via SSH or SSH2... ## $SSH_CLIENT and $SSH2_CLIENT are set automatically by SSH.

# If we're in TMUX, then change the screen type to "screen-256color" and export the TERM

export TERM=xterm-256color # override the defaults and always assume 256 colors

case "$TERM" in
    xterm-color|xterm-256color|screen-256color)	color_prompt=1 ;;
    *)	                        ;;
esac
[[ -n "$TMUX" ]] && [[ color_prompt==1 ]] && export TERM=screen-256color
if [[ "$OSTYPE" == darwin* ]] ; then isMac=1 ; fi

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

bind 'set mark-directories on'
bind 'set mark-symlinked-directories on'

# ============================= PATH STUFF ============================
export BOWTIE_INDEXES=/work/Apps/Bio/bowtie/current-bowtie/indexes/ ## <-- MUST have a trailing "/" after it!

if [[ -d "$HOME/TimeForScience" ]]; then
    ## If this is in the home directory, then set it no matter what.
    export TIME_FOR_SCIENCE_DIR="$HOME/TimeForScience"
elif [[ "$USER" == "alexgw" ]]; then
    ## Location of the TIME FOR SCIENCE directory. This is mostly Alex's code.
    export TIME_FOR_SCIENCE_DIR="$HOME/TimeForScience"
else
    export TIME_FOR_SCIENCE_DIR="/home/alexgw/TimeForScience"
    if [[ -d ${TIME_FOR_SCIENCE_DIR} ]]; then
	true # Don't do anything; it was found, which is fine
    else
	export TIME_FOR_SCIENCE_DIR="COULD_NOT_FIND_TIME_FOR_SCIENCE_DIRECTORY_ON_FILESYSTEM"
    fi
fi
export MYPERLDIR=${TIME_FOR_SCIENCE_DIR}/Lab_Code/Perl/ ## <-- Used by Josh Stuart's scripts. Mostly these are enhanced perl versions of UNIX scripts, like "cut.pl" and "join.pl".


# PATH: The FIRST things get run first!
# Low priority: Normal UNIX paths. Clear out the initial path here.
export PATH=/usr/bin:/bin:/usr/sbin:/sbin:/usr/X11/bin:/sw/bin ## <-- clear out initial path!!
# Medium priority: specific installs of tools
export PATH=/usr/local/bin:/work/Apps/bin:/usr/local/sbin:/opt/bin:/projects/bin:${PATH}  # <-- PRIORITY: programs in your own home directory come FIRST, then System-wide "Science" bin, then other stuff.
# Highest priority: things in the HOME directory or TimeForScience or /bioinformatics/bin
export PATH="${HOME}/bin:${HOME}/.linuxbrew/bin:${HOME}/.linuxbrew/bin:$TIME_FOR_SCIENCE_DIR/Lab_Code/Perl/Tools:$TIME_FOR_SCIENCE_DIR/Lab_Code/Perl/Scientific:$TIME_FOR_SCIENCE_DIR/Lab_Code/Python/Tools:$TIME_FOR_SCIENCE_DIR/Lab_Code/Shell:$TIME_FOR_SCIENCE_DIR/Lab_Code/R:/data/bin:${PATH}"
export LD_LIBRARY_PATH="${HOME}/.linuxbrew/lib:/data/lib:$LD_LIBRARY_PATH"

#export PERL5LIB=/usr/local/lib/perl5/site_perl/5.12.1/  # note: sometimes it's in /usr/lib instead of /usr/local!!

### ============== Below: Added Nov 27, 2012 for MRC metagenomics classifier project ======================
### Setup stuff for MRC for Tom Sharpton
#export PKG_CONFIG_PATH="/usr/share/doc/libgenome-1.3-0:$PKG_CONFIG_PATH"
#export PERL5LIB=/home/sharpton/dev/bioperl-live:/home/sharpton/sandbox/bioperl-run/lib:/home/sharpton/sandbox/bioperl-live/:/home/sharpton/src/bioperl-live:/home/sharpton/src/bioperl-live/t/lib:/lib:/home/sharpton/projects/MRC/lib:/home/sharpton/projects/MetaPASSAGE/simPipeModules:$PERL5LIB:/home/sharpton/projects/MRC/scripts/:/home/sharpton/src/netcdf-perl-1.2.4/src/perl/:/home/sharpton/projects/MRC/db_loader/lib:/home/sharpton/projects/matchmaking/scripts/lib:/home/sharpton/scripts/perl_modules/
#export WNHOME=/home/sharpton/projects/s_clustering/data/WNdb-3.0.tar.gz
### ============== Above: Added Nov 27, 2012 for MRC metagenomics classifier project ======================

### ============= For VCFtools, Jan 28, 2016:
#brew install vcftools
#==> Caveats  To use the Perl modules, make sure Vcf.pm, VcfStats.pm, and FaSlice.pm are included in your PERL5LIB environment variable:
export PERL5LIB=/data/home/alexgw/.linuxbrew/lib/perl5/site_perl:${PERL5LIB}
export PERLLIB=${PERL5LIB}
### =============


if [[ "$COMPYNAME" == "Slithereens" ]]; then
    export BINF_CORE_WORK_DIR="/Users/${USER}/work" # <-- set BINF_CORE work directory
elif [[ "$HOSTNAME" == "westway" ]] || [[ "$HOSTNAME" == "bueno" ]]; then
    export BINF_CORE_WORK_DIR="/home/alexgw/work"
elif [[ "$COMPYNAME" == "lighthousewww" ]] || [[ "$COMPYNAME" == "lighthouse" ]]; then
    export BINF_CORE_WORK_DIR="/home/awilliams/work"
else
    export BINF_CORE_WORK_DIR="/work"  # <-- set BINF_CORE work directory
fi

export PATH="${PATH}:$BINF_CORE_WORK_DIR/Common/Code/Python:$BINF_CORE_WORK_DIR/Common/Code/alexgw"

export R_BINF_CORE="$BINF_CORE_WORK_DIR/Common/Code/R_Binf_Core"

# ============================= DONE WITH PATH STUFF ============================

function agw_cmd_exists() {
    # or try: if [[ -n `which exa 2> /dev/null` ]] ...
    type "$1" &> /dev/null ;
    #export -f agw_cmd_exists
    # This returns something that will show up as true
    # in an "if" statement if a command exists, and false otherwise.
    # Note: the exit code is backwards--0 = exists, 1 = nope
}

if agw_cmd_exists "ls"; then
    echo "ls appears to exist"
fi

if agw_cmd_exists "fdsfds"; then
    echo "fdsfds appears to exist"
else
    echo "fdsfds does not exist"
fi


# Alias definitions.
# You may want to put all your additions into a separate file like
# ~/.bash_aliases, instead of adding them here directly.
# See /usr/share/doc/bash-doc/examples in the bash-doc package.
# Sets up the aliases whether or not this shell is interactive
if [[ -f ${TIME_FOR_SCIENCE_DIR}/Config/Alex_Williams/.aliases ]]; then
    source ${TIME_FOR_SCIENCE_DIR}/Config/Alex_Williams/.aliases
elif [[ -f ${BINF_CORE_WORK_DIR}/Common/Code/TimeForScience/Config/Alex_Williams/.aliases ]] ; then
    source ${BINF_CORE_WORK_DIR}/Common/Code/TimeForScience/Config/Alex_Williams/.aliases
elif [[ -f ${HOME}/TimeForScience/Config/Alex_Williams/.aliases ]] ; then
    source ${HOME}/TimeForScience/Config/Alex_Williams/.aliases
fi

# ======== SET THE COMMAND PROMPT COLOR FOR THIS MACHINE ======== #
a_machine_prompt_main_color=''
if [[ -n "$color_prompt" ]] ; then
    if [[ -z "$is_sshing" ]] ; then
	## LOCAL machine gets a specific color...
	#a_machine_prompt_main_color="$color_prefix[44m$color_prefix[3;33m" ## Blue background, white foreground
	a_machine_prompt_main_color="${color_prefix}[1;32m" ## Bold green text
	iterm_bg=000000 # iTerm console window background
	iterm_top="180 180 180"
    else
	## Ok, we are SSHed into a remote machine...
	case "$COMPYNAME"
	    in
	    $BN_IP)
		a_machine_prompt_main_color="${color_prefix}[1;35m" ; ## 1;35m == Bold magenta
		iterm_bg=002200 ;# iTerm console window background
		iterm_top="40 120 40" ; # iTerm window top bar color
		;;
	    $BC_IP)
		a_machine_prompt_main_color="${color_prefix}[1;34m" # ???
		iterm_bg=000022 ;
		iterm_top="80 80 250" ;
		;;
	    $PB_IP)
		a_machine_prompt_main_color="${color_prefix}[44m${color_prefix}[3;36m" # cyan text / blue background
		iterm_bg=220000 ;
		iterm_top="140 40 40" ;
		;;
	    $PL_IP)
		a_machine_prompt_main_color="${color_prefix}[43m${color_prefix}[3;30m" # yellow background
		iterm_bg=302000 ;
		iterm_top="120 100 40" ;
		;;
	    --) ## OTHER unknown machine
		a_machine_prompt_main_color="${color_prefix}[41m${color_prefix}[3;37m"
		iterm_bg=333333 ;
		iterm_top="120 50 120" ;
		;;
	esac
    fi

    #\033[47m\033[3;31m   ## <-- For the second part, where it says [3;31m 1=bold, 2=light, 3=regular
fi
# ======== SET THE COMMAND PROMPT COLOR FOR THIS MACHINE ======== #

#if [[ -z "$STY" ]] && [[ "$HOSTNAME" == 'something' ]] ; then
#    echo -e "${a_warning_color}*\n*\n* Remember to resume the screen session on this machine!!!\n*\n*${a_end_color}"
#else
#    echo -n ''
#fi

## ===============================================
## ====== TERMINAL HISTORY =======================
export HISTSIZE=30123
export HISTCONTROL=ignoredups # Do not save duplicate/blank lines to the history
export HISTIGNORE="ls:kpk:exit:p:pwd:rr:clear:history:fg:bg" ## Commands that are NOT saved to the history!
set revert-all-at-newline on # <-- prevents editing of history lines! If this gets turned off, it is SUPER annoying.
shopt -s histappend # Save terminal history between sessions
PROMPT_COMMAND='history -a' ## save ALL terminal histories
## ====== TERMINAL HISTORY =======================
## ===============================================

stty -ixon  # Disable the totally useless START/STOP output control (enables you to pause input by pressing the Ctrl-S key sequence and resume output by pressing the Ctrl-Q key sequence)

shopt -s globstar # With globstar set (bash 4.0+), bash recurses all the directories. In other words, enables '**'

set   -o ignoreeof  # Prevent Ctrl-D from exiting! Still exits if you press it 10 times.
shopt -s checkwinsize # Check the window size after each command and, if necessary, update the values of LINES and COLUMNS.
shopt -s cdspell ## Fix common mis-spellings in directories to "cd" to
shopt -s cmdhist ## Save multi-line pasted commands into one single history command
#shopt -s lithist ##
shopt -s no_empty_cmd_completion ## Don't display ALL commands on an empty-line tab
shopt -s nocaseglob ## Match glob / regexp in case-insensitive fashion

# set variable identifying the chroot you work in (used in the prompt below)
if [[ -z $debian_chroot ]] && [[ -r /etc/debian_chroot ]]; then
    debian_chroot=$(cat /etc/debian_chroot)
fi

#PS1="[\D{%e}=\t \h:\W]$ "
PS1="[\D{%e}~\t~${COMPYNAME:0:3}~\W]$ " ## ${HOSTNAME:0:3} means only show the first 3 characters of the hostname! "\h" is the whole thing, also.

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

## ===============================================
## ====== LS COLORS ==============================
## COMMAND LINE COLOR / LS COLOR
export CLICOLOR='Yes'
export LS_OPTIONS='--color=auto'

## LS_COLORS *with* an underscore is for Ubuntu
export LS_COLORS='no=00:fi=00:di=01;35:ln=01;36:pi=40;33:so=01;35:do=01;35:bd=40;33;01:cd=40;33;01:or=43;31;03:ex=01;31' 

## LSCOLORS *without* an underscore is for Mac OS X.
export LSCOLORS=FxgxCxDxBxegedabagacad

#if [[ -n "$isMac" ]] ; then
#    AGW_LS_COLOR_OPTION=" -G "  ## -G means "color" on the mac...
#else
#    AGW_LS_COLOR_OPTION=' --color=auto ' 
#fi
## ====== LS COLORS ==============================
## ===============================================

export HOSTNAME  ## <-- Required in order for tmux / env to see HOSTNAME as a variable!
export    CVS_RSH=ssh
export  CVSEDITOR="emacs -nw"
export SVN_EDITOR="emacs -nw"
export     EDITOR="emacs -nw"

export LANG=C    # Set the LANG to C. Speeds some command line tools up, but *BREAKS* some UTF-8 stuff!

export CXX=/usr/bin/g++ # The location of the c++ compiler. NOT "cpp"--that is the C preprocessor.

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

# Save the local ethernet "en0" MAC address into the variable LOCAL_EN0_MAC. Note the zero.
# Allows per-machine settinsg.
#export LOCAL_EN0_MAC=`ifconfig en0 | grep -i ether | sed 's/.*ether //' | sed 's/[ ]*$//'`

umask u=rwx,g=rwx,o=rx # <-- give users and groups full access to files I create, and let other users READ and EXECUTE
#umask 0007 # <-- give users and groups full access to files I create, but give no access to other users.
# 0 = "full access", 7 = "no access"


# =======

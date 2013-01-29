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

if [[ -f ~/TimeForScience/Config/Alex_Williams/bash-platform-specific ]] ; then
    source ~/TimeForScience/Config/Alex_Williams/bash-platform-specific
elif [[ -f /work/Common/Code/TimeForScience/Config/Alex_Williams/bash-platform-specific ]] ; then
    source /work/Common/Code/TimeForScience/Config/Alex_Williams/bash-platform-specific
fi

echo -e "${a_echo_color}>>> BASH: Loading .bashrc...${a_end_color}" ## <-- comes after the colors are set up in platform-specific fashion

if [[ -f ~/TimeForScience/Config/Alex_Williams/bash-path-setup ]] ; then
    source ~/TimeForScience/Config/Alex_Williams/bash-path-setup
elif [[ -f /work/Common/Code/TimeForScience/Config/Alex_Williams/bash-path-setup ]] ; then
    source /work/Common/Code/TimeForScience/Config/Alex_Williams/bash-path-setup
fi

# Alias definitions.
# You may want to put all your additions into a separate file like
# ~/.bash_aliases, instead of adding them here directly.
# See /usr/share/doc/bash-doc/examples in the bash-doc package.
# Sets up the aliases whether or not this shell is interactive
if [[ -f ~/TimeForScience/Config/Alex_Williams/.aliases ]]; then
    source ~/TimeForScience/Config/Alex_Williams/.aliases
elif [[ -f /work/Common/Code/TimeForScience/Config/Alex_Williams/.aliases ]] ; then
    source /work/Common/Code/TimeForScience/Config/Alex_Williams/.aliases
fi

#iterm_bg=000000
#iterm_title_r255=200
#iterm_title_g255=0
#iterm_title_b255=0

# function iTermColor {
# # Set the Mac application "iTerm2" to have a custom title bar or a custom background
# # Argument $1 is either "top" or "bg"
# ## If the argument $1 is "top", then:
# ##   $2/$3/$4 are R/G/B, which can range from 0 to 255
# ##   This will, on the mac, set the iTerm2 title bar's color
# ## If the aargument $1 is "bg", then the input $2 is a single hex color string in 00-FF format.
# ## Example, it could be FF00FF, for magenta.
#     if [[ "$1" == "top" ]] ; then
# 	[[ -n "$2" ]] && echo -e "\033]6;1;bg;red;brightness;$2\a" ; ## title bar of iTerm
# 	[[ -n "$3" ]] && echo -e "\033]6;1;bg;green;brightness;$3\a" ; ## title bar of iTerm
# 	[[ -n "$4" ]] && echo -e "\033]6;1;bg;blue;brightness;$4\a" ## title bar of iTerm
#     elif [[ "$1" == "bg" ]] ; then
# 	[[ -n "$2" ]] && echo -e "\033]Ph$2\033\\" ; ## Set the iTerm background. $2 must be like FF00FF
#     fi
# }

# function agwReColorTerminal {
#     iTermColor "bg" ${iterm_bg} ; iTermColor "top" ${iterm_top}
#     echo "Re-coloring the terminal..."
# }


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
	case "$HOSTNAME"
	    in
	    nausicaa)
		a_machine_prompt_main_color="${color_prefix}[1;35m" ; ## 1;35m == Bold magenta
		iterm_bg=002200 ;# iTerm console window background
		iterm_top="40 120 40" ; # iTerm window top bar color
		;;
	    catbus)
		a_machine_prompt_main_color="${color_prefix}[1;34m"
		iterm_bg=000022 ;
		iterm_top="80 80 250" ;
		;;
	    bueno)
		a_machine_prompt_main_color="${color_prefix}[44m${color_prefix}[3;36m" # cyan text / blue background
		iterm_bg=220000 ;
		iterm_top="140 40 40" ;
		;;
	    lighthouse)
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

#agwReColorTerminal

if [[ -z "$STY" ]] && [[ "$HOSTNAME" == 'nausicaa' ]] ; then
    echo -e "${a_warning_color}*\n*\n* Remember to resume the screen session on this machine!!!\n*\n*${a_end_color}"
else
    echo -n ''
fi

## ===============================================
## ====== TERMINAL HISTORY =======================
export HISTSIZE=30123
export HISTCONTROL=ignoredups # duplicate and blank lines are not saved in the history
export HISTIGNORE="ls:ll:kpk:exit:p:pwd:rr:clear:history:fg:bg" ## Commands that are NOT saved to the history!
set revert-all-at-newline on # <-- prevents editing of history lines! This is SUPER important, and SUPER annoying if it gets set to 'off'
shopt -s histappend # Save terminal history between sessions
PROMPT_COMMAND='history -a' ## save ALL terminal histories
## ====== TERMINAL HISTORY =======================
## ===============================================

stty -ixon -ixoff  # Disable the totally useless START/STOP output control (enables you to pause input by pressing the Ctrl-S key sequence and resume output by pressing the Ctrl-Q key sequence)

set   -o ignoreeof  # Prevent Ctrl-D from exiting!
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
PS1="[\D{%e}/\t:${HOSTNAME:0:3}:\W]$ " ## ${HOSTNAME:0:3} means only show the first 3 characters of the hostname! "\h" is the whole thing, also.

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
## LSCOLORS *without* an underscore is for Mac OS X
export LSCOLORS=FxgxCxDxBxegedabagacad
#if [[ -n "$isMac" ]] ; then
#    AGW_LS_COLOR_OPTION=" -G "  ## -G means "color" on the mac...
#else
#    AGW_LS_COLOR_OPTION=' --color=auto ' 
#fi
## ====== LS COLORS ==============================
## ===============================================

export HOSTNAME  ## <-- Required in order for tmux / env to see HOSTNAME as a variable!
export CVS_RSH=ssh
export CVSEDITOR="emacs -nw"
export SVN_EDITOR="emacs -nw"
export EDITOR="emacs -nw"

# Keybindings: Add IJKL navigation to supplement/replace the arrow keys
bind "\M-J:backward-word"
bind "\M-L:forward-word"
bind "\M-j:backward-char"
bind "\M-l:forward-char"
bind "\M-i:previous-history"
bind "\M-I:previous-history"
bind "\M-k:next-history"
bind "\M-K:next-history"
bind "\C-s:history-search-forward"

# Save the local ethernet "en0" MAC address into the variable LOCAL_EN0_MAC. Note the zero.
# Allows per-machine settinsg.
#export LOCAL_EN0_MAC=`ifconfig en0 | grep -i ether | sed 's/.*ether //' | sed 's/[ ]*$//'`

umask u=rwx,g=rwx,o=rx # <-- give users and groups full access to files I create, and let other users READ and EXECUTE
#umask 0007 # <-- give users and groups full access to files I create, but give no access to other users.
# 0 = "full access", 7 = "no access"


# =======

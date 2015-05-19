#function fish_prompt
#	 set_color $fish_color_cwd
#	 echo -n (prompt_pwd)
#	 set_color normal
#	 echo -n ' > '
#end

#[[ "$color_prompt" && "$isMac" ]] && AGW_LS_OPT=' -F -G -@ ' ## Mac: Color option is -G. Also, on a Mac, show the extended attributes (-@)
#[[ "$color_prompt" && (-z "$isMac") ]] && AGW_LS_OPT=' --indicator-style=slash --color=auto ' ## Ubuntu: color is --color=auto

set -l AGW_LS_OPT=' -F -G -@ '
alias ls='/bin/ls {$AGW_LS_OPT}'
alias l='ls'

alias ll='/bin/ls -l -h -A -F {$AGW_LS_OPT} '
alias lo='ll -r -t' # Depends on "ll" already being defined above. List by TIME so that NEWEST files are at the bottom

alias c='cd'
alias p='pwd -P'
alias tc='randomize_terminal_color.pl -cycle'
#alias cd='pushd'
alias b='popd'
alias ..='cd ..'
#alias ..='pushd ..'
alias recrash='sudo launchctl unload /Library/LaunchDaemons/com.crashplan.engine.plist; sleep 1; sudo launchctl load /Library/LaunchDaemons/com.crashplan.engine.plist'
alias res='source ~/.bashrc ; source ~/.bash_profile'
alias mv='mv -i'
alias cp='cp -i'
alias kpk='exit'
alias diffdir='diff -rq' ## diff on directories


## Mac-specific commands:
set -l SETFILE_LOCATION="/Applications/Xcode*/Contents/Developer/Tools/SetFile"
alias invisible="chflags -h hidden"
alias visible="chflags -h nohidden"
alias clearicon="{$SETFILE_LOCATION} -P -a c"
alias version="lsb_release -a"

#alias rr='tmux attach || tmux new' #"screen -xR -U" ## Reconnect to the previous screen, or make a new one if there isn't one already
alias sc='tmux new-window' #'screen' ## Make a new screen session

alias wcl='wc -l'
alias sortt='sort -t "	"' # sort with tab as separator
alias sortg='sort -g -t "	"' # sort NUMBERS, with tab as separator

alias e="emacsclient -nw  -c --alternate-editor=\"\"" # "As a special exception, if command is the empty string, then emacsclient starts Ema

alias make='make --warn-undefined-variables --print-directory'
#alias mcm='make clean && make'

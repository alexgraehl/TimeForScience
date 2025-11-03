# If you come from bash you might have to change your $PATH.
# export PATH=$HOME/bin:/usr/local/bin:$PATH

# Path to your oh-my-zsh installation.
export ZSH="$HOME/.oh-my-zsh"

MAYBE_LOCAL_CONFIG_ONLY_ZSH_CFG="$HOME/.mocal-config-zsh.sh"
if [ -f $MAYBE_LOCAL_CONFIG_ONLY_ZSH_CFG ]; then
    echo "[:OK:] Loading local-only configs $MAYBE_LOCAL_CONFIG_ONLY_ZSH_CFG..."
    source $MAYBE_LOCAL_CONFIG_ONLY_ZSH_CFG
else
    echo "[:HEY:] NOT Loading nonexistent file $MAYBE_LOCAL_CONFIG_ONLY_ZSH_CFG..."
fi


# Set name of the theme to load --- if set to "random", it will
# load a random theme each time oh-my-zsh is loaded, in which case,
# to know which specific one was loaded, run: echo $RANDOM_THEME
# See https://github.com/ohmyzsh/ohmyzsh/wiki/Themes
ZSH_THEME="robbyrussell"

# Set list of themes to pick from when loading at random
# Setting this variable when ZSH_THEME=random will cause zsh to load
# a theme from this variable instead of looking in $ZSH/themes/
# If set to an empty array, this variable will have no effect.
# ZSH_THEME_RANDOM_CANDIDATES=( "robbyrussell" "agnoster" )

# Uncomment the following line to use case-sensitive completion.
# CASE_SENSITIVE="true"

# Uncomment the following line to use hyphen-insensitive completion.
# Case-sensitive completion must be off. _ and - will be interchangeable.
# HYPHEN_INSENSITIVE="true"

# Uncomment one of the following lines to change the auto-update behavior
# zstyle ':omz:update' mode disabled  # disable automatic updates
# zstyle ':omz:update' mode auto      # update automatically without asking
# zstyle ':omz:update' mode reminder  # just remind me to update when it's time

# Uncomment the following line to change how often to auto-update (in days).
# zstyle ':omz:update' frequency 13

# Uncomment the following line if pasting URLs and other text is messed up.
# DISABLE_MAGIC_FUNCTIONS="true"

# Uncomment the following line to disable colors in ls.
# DISABLE_LS_COLORS="true"

# Uncomment the following line to disable auto-setting terminal title.
# DISABLE_AUTO_TITLE="true"

# Uncomment the following line to enable command auto-correction.
# ENABLE_CORRECTION="true"

# Uncomment the following line to display red dots whilst waiting for completion.
# You can also set it to another string to have that shown instead of the default red dots.
# e.g. COMPLETION_WAITING_DOTS="%F{yellow}waiting...%f"
# Caution: this setting can cause issues with multiline prompts in zsh < 5.7.1 (see #5765)
# COMPLETION_WAITING_DOTS="true"

# Uncomment the following line if you want to disable marking untracked files
# under VCS as dirty. This makes repository status check for large repositories
# much, much faster.
# DISABLE_UNTRACKED_FILES_DIRTY="true"

# Uncomment the following line if you want to change the command execution time
# stamp shown in the history command output.
# You can set one of the optional three formats:
# "mm/dd/yyyy"|"dd.mm.yyyy"|"yyyy-mm-dd"
# or set a custom format using the strftime function format specifications,
# see 'man strftime' for details.
# HIST_STAMPS="mm/dd/yyyy"

# Would you like to use another custom folder than $ZSH/custom?
# ZSH_CUSTOM=/path/to/new-custom-folder

# Which plugins would you like to load?
# Standard plugins can be found in $ZSH/plugins/
# Custom plugins may be added to $ZSH_CUSTOM/plugins/
# Example format: plugins=(rails git textmate ruby lighthouse)
# Add wisely, as too many plugins slow down shell startup.
plugins=(git)


function deduped() { # input: one string to de-dupe. Usage: PATH=$(deduped $PATH). May break on things with spaces.
    local X="$1"
    local new=""
    local old="$X:"
    if [ -n "$X" ]; then
	while [ -n "$old" ]; do
	    x=${old%%:*}
	    case "$new": in
		    *:"$x":*) ;;        # <-- This element exists already (do not append it)
		    *) new="$new:$x" ;; # <-- Does not exist already, so append "new"
	    esac
	    old=${old#*:}
	done
	new=${new#:}
    fi
    new=$(echo "$new" | perl -pe 's/^[:+]//' | perl -pe 's/[:]+/:/g' | perl -pe 's/[:+]\$//') # remove leading/trailing ':'
    # Remove any leading or trailing ':' from the path
    echo "$new"
}


if [[ -d "${HOME}/TimeForScience" ]]; then ## If this is in the home directory, then set it no matter what.
    export TIME_FOR_SCIENCE_DIR="${HOME}/TimeForScience"
else
    export TIME_FOR_SCIENCE_DIR="/home/alexgw/TimeForScience"
    if [[ -d ${TIME_FOR_SCIENCE_DIR} ]]; then
	true # Do not't do anything; it was found, which is fine
    else
	echo "[:HEY:] WARNING in .bashrc: TIME_FOR_SCIENCE_DIR not found]"
	export TIME_FOR_SCIENCE_DIR="COULD_NOT_FIND_TIME_FOR_SCIENCE_DIRECTORY_ON_FILESYSTEM"
    fi
fi


DISABLE_AUTO_UPDATE="true"
source $ZSH/oh-my-zsh.sh
source $HOME/.aliases-zsh

# Add `flutter` to the path, per https://docs.flutter.dev/get-started/install/macos
export PATH="$HOME/bin:/opt/homebrew/bin:/opt/homebrew/sbin:$PATH:/Users/alexgw/Applications/flutter/bin"

# Added by LM Studio CLI (lms)
export PATH="$PATH:/Users/alexgw/.cache/lm-studio/bin"

if [ -f ~/.local-config-no-git.sh ]; then
    source ~/.local-config-no-git.sh
else
    echo -e "[:FYI:] No ~/.local-config-no-git.sh was found or loaded."
fi


# User configuration

# export MANPATH="/usr/local/man:$MANPATH"

# You may need to manually set your language environment
# export LANG=en_US.UTF-8

# Preferred editor for local and remote sessions
# if [[ -n $SSH_CONNECTION ]]; then
#   export EDITOR='vim'
# else
#   export EDITOR='mvim'
# fi

# Compilation flags
# export ARCHFLAGS="-arch x86_64"

# Set personal aliases, overriding those provided by oh-my-zsh libs,
# plugins, and themes. Aliases can be placed here, though oh-my-zsh
# users are encouraged to define aliases within the ZSH_CUSTOM folder.
# For a full list of active aliases, run `alias`.
#
# Example aliases
# alias zshconfig="mate ~/.zshrc"
# alias ohmyzsh="mate ~/.oh-my-zsh"


export NVM_DIR="$HOME/.nvm"
[ -s "$NVM_DIR/nvm.sh" ] && \. "$NVM_DIR/nvm.sh"  # This loads nvm
[ -s "$NVM_DIR/bash_completion" ] && \. "$NVM_DIR/bash_completion"  # This loads nvm bash_completion

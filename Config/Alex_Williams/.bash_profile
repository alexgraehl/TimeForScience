
# This is the ~/.bash_profile that should be in place for any Macintosh system!

source ~/.bashrc




# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/Users/alexgw/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/Users/alexgw/miniconda3/etc/profile.d/conda.sh" ]; then
        . "/Users/alexgw/miniconda3/etc/profile.d/conda.sh"
    else
        export PATH="/Users/alexgw/miniconda3/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<


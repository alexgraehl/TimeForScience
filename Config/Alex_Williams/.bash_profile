
# This is the ~/.bash_profile that should be in place for any Macintosh system!

source ~/.bashrc

# The next line updates PATH for the Google Cloud SDK.
#if [ -f '/Users/alexwilliams/Downloads/google-cloud-sdk/path.bash.inc' ]; then source '/Users/alexwilliams/Downloads/google-cloud-sdk/path.bash.inc'; fi

# The next line enables shell command completion for gcloud.
#if [ -f '/Users/alexwilliams/Downloads/google-cloud-sdk/completion.bash.inc' ]; then source '/Users/alexwilliams/Downloads/google-cloud-sdk/completion.bash.inc'; fi

if [[ -f "$HOME/.modules" ]]; then
    echo "[MODULES]: loading..."
    time source "$HOME/.modules"
    echo "[MODULES]: Done loading"
fi

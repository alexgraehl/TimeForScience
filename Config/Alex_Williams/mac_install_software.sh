#!/usr/bin/env bash -u
set -e
set -o pipefail
set -o xtrace

cd "$HOME"

SHOULD_HOMEBREW=1

if [[ -z "${TIME_FOR_SCIENCE_DIR:-}" ]] || [[ ! -d "$TIME_FOR_SCIENCE_DIR" ]] || [[ ! -e "$TIME_FOR_SCIENCE_DIR/README.md" ]]; then
    echo "[FYI]: TIME_FOR_SCIENCE_DIR is not set. You can set it by following the instructions at https://github.com/alexgraehl/TimeForScience.";
    echo "[FYI]: Continuing without symlinking..."
else
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #            Create symlinks to standard config files
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if [[ $(dirname $TIME_FOR_SCIENCE_DIR) == "$HOME" ]]; then echo "use relative links"; else echo "do not use relative links--apparently the TimeForScience directory is NOT located immediately in your home directory2. Exiting here, as we cannot deal with this."; exit 1; fi

    HOMEDIR_CONFIG_FILES_TO_SYMLINK=(.aliases .bash_profile .bashrc .emacs .inputrc .screenrc .tmux.conf)
    for FILE in ${HOMEDIR_CONFIG_FILES_TO_SYMLINK[@]}; do
	TIME_FOR_SCIENCE_BASENAME=$(basename $TIME_FOR_SCIENCE_DIR)
	CFG="./${TIME_FOR_SCIENCE_BASENAME}/Config/Alex_Williams/$FILE"
	echo "CFG="$CFG
	echo "FILE="$FILE
	if [[ ! -e $CFG ]]; then echo "[ERROR] The config file that we expected to exist, specifically, '$CFG', did not exist."; exit 1; else echo "[OK] Linking $CFG to the home directory (via relative path)"; fi
	ln -f -s $CFG ./
    done

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #                          SSH Config
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if [[ -e "$HOME/.ssh/config" ]]; then
	echo "[:HEY:] [SKIPPING] The ssh config file at .ssh/config ALREADY exists, so we are not overwriting it."
    else
	cp ${TIME_FOR_SCIENCE_DIR}/Config/Alex_Williams/ssh-slash-config $HOME/.ssh/config
	chmod -R go-w   ~/.ssh/
	chmod    go-rwx ~/.ssh/id_rsa
    fi
fi

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                         Homebrew
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if [[ 1 == $SHOULD_HOMEBREW ]]; then
    echo "Installing a ton of stuff via homebrew (which must already be installed in order for this to work)"

    brew cask install java # required for gradle and more
    brew cask install tldr # Better examples, like 'man' or 'info'
    #brew cask install iterm2 # hexyl glances  # hexyl: command line color hex editor. glances: fancier htop
    brew cask install github
    brew install bash bash-completion cairo catimg dialog \
	 emacs ffmpeg fontconfig freetype gcc gdbm gettext git glances glib gmp gnutls gradle gsl \
	 htop imagemagick ipython isl jpeg \
	 lame libevent libffi libgit2 libmpc libpng libtasn1 libtiff libtool libvo-aacenc \
	 mono mpfr nettle numpy openssl pcre pigz pixman pyenv pkg-config pv python qemu \
	 R readline rsync scipy sqlite szip tmux watch wget wxmac wxpython x264 xvid xz

    if [[ ! -e "/Applications/iTerm.app"              ]]; then brew cask install iterm2; fi
    if [[ ! -e "/Applications/Stats.app"              ]]; then brew install stats; fi   # Menu bar stats
    if [[ ! -e "/Applications/Visual Studio Code.app" ]]; then brew cask install visual-studio-code; fi
    if [[ ! -e "/Applications/Spotify.app"            ]]; then brew cask install spotify; fi
    if [[ ! -e "/Applications/BetterTouchTool.app"    ]]; then brew cask install bettertouchtool; fi
    if [[ ! -e "/Applications/Dropbox.app"            ]]; then brew install dropbox; fi
    if [[ ! -e "/Applications/LibreOffice.app"        ]]; then brew cask install libreoffice; fi
    if [[ ! -e "/Applications/Atom.app"               ]]; then brew cask install atom; fi  # markdown editor
fi

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
AGW_BASH_MAJOR_VERSION=$(echo "$BASH_VERSION" | perl -pe "s/(\d+).*/\1/; chomp;")
if (( ${AGW_BASH_MAJOR_VERSION} >= 4 )); then
    echo "[OK] globstar should be supported on bash version $AGW_BASH_MAJOR_VERSION, which is >= 4."
else
    echo "[GOTTA FIX] globstar is NOT supported on bash version $AGW_BASH_MAJOR_VERSION. This is probably because the Apple-included bash is ancient (3.x). Let's now use Homebrew to set a new bash."
    ALLOWED_TO_INSTALL_BASH_4_ON_MAC=1
    if [[ "${ALLOWED_TO_INSTALL_BASH_4_ON_MAC}" == "1" ]]; then
	brew install bash
	NEW_BASH_LOCATION=/usr/local/bin/bash
	echo "Assuming that homebrew has installed bash to the location ${NEW_BASH_LOCATION}."
	if [[ -f ${NEW_BASH_LOCATION} ]]; then
	    #grep $NEW_BASH_LOCATION /etc/shells
	    sudo bash -c "echo '${NEW_BASH_LOCATION}' >> /etc/shells"
	    sudo chsh -s ${NEW_BASH_LOCATION} ${USER}
	else
	    echo "Failure! The bash installation was not located in the expected location. Did homebrew fail to install bash? The expected location (which is seemingly NOT A VALID FILE) was supposed to be --> ${NEW_BASH_LOCATION}" && exit 1
	fi
    fi
fi

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
echo "Making a few required directories, in case they aren't already present."
mkdir -p ${HOME}/bin

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
echo "Making a bunch of uninteresting files/folders hidden in the default Finder view..."
FILES_TO_INVISIFY=("$HOME/Movies" "$HOME/Music" "$HOME/Public" "$HOME/TimeForScience" "$HOME/bin")
for FILE in ${FILES_TO_INVISIFY[@]}; do
    chflags -h hidden ${FILE}
done
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
echo "A bunch of things that should be installed via pip3 (Python 3)"
PIP3_PACKAGES=("matplotlib" "numpy" "scipy" "jupyter" "mypy" "scikit-learn" "docutils" "pandas" "pylint" "sphinx" "sphinx-autobuild")
for PACKAGE in ${PIP3_PACKAGES[@]}; do
    echo "[PYTHON3 (PIP3)] Installing the package $PACKAGE..."
    python3 -m pip install ${PACKAGE}
done
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#echo "Install the amazon command line tools (via python's pip)"
#pip install awscli
#echo "[Installed] the 'aws' tool"

echo "MAC APP STORE: install from the Mac App store:"
echo "      * Divvy"
echo "      * Slack"
echo "      * Pixelmator"
echo "      * XCode"

echo "YOU SHOULD MANUALLY INSTALL THE FOLLOWING PROGRAMS:"
echo "           cool-retro-term           FROM --> brew cask install cool-retro-term (retro-style Terminal)"
echo "           RStudio                   FROM --> https://www.rstudio.com/products/rstudio/download/#download"
echo "           SCROLL REVERSER           FROM --> https://pilotmoon.com/downloads/ScrollReverser-1.7.6.zip"
echo "           CLIPY (replaces ClipMenu) FROM --> https://clipy-app.com/"

echo "Remember to sync the iTerm2 preferences from my config folder on github!"

echo "You will need to MANUALLY run XCode to let it install its weird components"
echo "maybe something like this: sudo /Applications/Xcode.app/Contents/Developer/usr/bin/xcodebuild  -license accept"
#sudo installer -pkg /Applications/Xcode-beta.app/Contents/Resources/Packages/MobileDevice.pkg -target /
#sudo installer -pkg /Applications/Xcode-beta.app/Contents/Resources/Packages/MobileDeviceDevelopment.pkg -target /

echo "You probably want to set up pyenv: e.g., pyenv install 3.9.0 and then pyenv global 3.9.0"

echo "[DONE!] But remember to manually install the software shown above."

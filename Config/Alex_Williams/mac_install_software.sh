#!/usr/bin/env bash -u
set -e
set -o pipefail
set -o xtrace

cd "$HOME"

if [[ -z "$TIME_FOR_SCIENCE_DIR" ]]; then echo "[ERROR: a required variable is NOT set!] You need to set the TIME_FOR_SCIENCE_DIR to the location of the TimeForScience directory (probably in your home directory). You could do this by following the instructions at https://github.com/alexgraehl/TimeForScience."; exit 1; else echo "[OK] The TIME_FOR_SCIENCE_DIR environment variable is '$TIME_FOR_SCIENCE_DIR'"; fi

if [[ ! -d "$TIME_FOR_SCIENCE_DIR" ]]; then echo "[ERROR: the TIME_FOR_SCIENCE_DIR does not exist (or is not a directory)] The TIME_FOR_SCIENCE_DIR environment variable must point to a valid checked-out copy of the TimeForScience repository. It appears to not be a directory. The directory we are looking for (which you should verify is correct) is: '$TIME_FOR_SCIENCE_DIR'"; exit 1; else echo "[OK] Successfully verified the existence of the directory '$TIME_FOR_SCIENCE_DIR'"; fi

if [[ ! -e "$TIME_FOR_SCIENCE_DIR/README.md" ]]; then echo "[ERROR: the TIME_FOR_SCIENCE_DIR appears NOT to contain a 'README.md', which implies that it is not a valid checked-out copy of the TimeForScience repository."; exit 1; else echo "[OK] The $TIME_FOR_SCIENCE_DIR/README.md file was successfully found."; fi

SHOULD_HOMEBREW=1

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#            Create symlinks to standard config files
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if [[ $(dirname $TIME_FOR_SCIENCE_DIR) == "$HOME" ]]; then echo "use relative links"; else echo "do not use relative links--apparently the TimeForScience directory is NOT located immediately in your home directory2. Exiting here, as we cannot deal with this."; exit 1; fi

HOMEDIR_CONFIG_FILES_TO_SYMLINK=(.aliases .bash_profile .bashrc .emacs .inputrc .screenrc .tmux.conf)
for FILE in ${HOMEDIR_CONFIG_FILES_TO_SYMLINK[@]}; do
	TIME_FOR_SCIENCE_BASENAME=$(basename $TIME_FOR_SCIENCE_DIR)
	CFG="./$TIME_FOR_SCIENCE_BASENAME/Config/Alex_Williams/$FILE"
	echo "CFG="$CFG
	echo "FILE="$FILE
	if [[ ! -e $CFG ]]; then echo "[ERROR] The config file that we expected to exist, specifically, '$CFG', did not exist."; exit 1; else echo "[OK] Linking $CFG to the home directory (via relative path)"; fi
	ln -f -s $CFG ./
done

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                          SSH Config
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if [[ -e "$HOME/.ssh/config" ]]; then
    echo "[:ERR:] [SKIPPING] The ssh config file at .ssh/config ALREADY exists, so we are not overwriting it."
else
    cp TimeForScience/Config/Alex_Williams/ssh-slash-config $HOME/.ssh/config
    chmod -R go-w   ~/.ssh/
    chmod    go-rwx ~/.ssh/id_rsa
fi

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                         Homebrew
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if [[ 1 == $SHOULD_HOMEBREW ]]; then
    echo "Installing a ton of stuff via homebrew (which must already be installed in order for this to work)"
    brew cask install java # required for gradle and more
    brew cask install iterm2
    brew cask install slack notion visual-studio-code github spotify transmit bettertouchtool tableplus
    brew install bash bbedit cairo catimg dialog dropbox emacs ffmpeg fontconfig freetype gcc gdbm gettext git glances glib gmp gnutls gradle gsl htop imagemagick ipython isl jpeg lame libevent libffi libgit2 libmpc libpng libtasn1 libtiff libtool libvo-aacenc mono mpfr nettle numpy openssl pcre pigz pixman pkg-config pv python qemu R readline rsync samtools scipy sqlite szip tmux watch wget wxmac wxpython x264 xvid xz
    brew install hexyl # hexyl: command line color hex editor. glances: fancier htop
    brew cask install atom # markdown (.md) viewer/editor
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
echo "A bunch of things that should be installed via pip2 (Python 2)"
PIP2_PACKAGES=("awscli" "pylint" "matplotlib" "numpy" "scipy")
for PACKAGE in ${PIP2_PACKAGES[@]}; do
    echo "[PYTHON2 (PIP2)] Installing the package $PACKAGE..."
    python  -m pip install ${PACKAGE}
done
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
echo "A bunch of things that should be installed via pip3 (Python 3)"
PIP3_PACKAGES=("awscli" "pylint" "matplotlib" "numpy" "scipy" "jupyter")
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
echo "           RStudio                   FROM --> https://www.rstudio.com/products/rstudio/download/#download"
echo "           SCROLL REVERSER           FROM --> https://pilotmoon.com/downloads/ScrollReverser-1.7.6.zip"
echo "           PYCHARM                   FROM --> https://www.jetbrains.com/pycharm/"
echo "           CLIPY (replaces ClipMenu) FROM --> https://clipy-app.com/"
echo "           OMNIGRAFFLE (Paid!)       FROM --> https://www.omnigroup.com/"

echo "UCSC Genome Browser software can be obtained for a Mac with this command: (as documented on the download page: http://hgdownload.cse.ucsc.edu/admin/exe/macOSX.x86_64/ )"
echo "mkdir -p ~/bin/ucsc && cd ~/bin/ucsc && rsync -aP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/macOSX.x86_64/ ./  (requires rsync >= 3.1, so get that from homebrew if you're using the ancient Mac version of rsync)"
echo "Remember to sync the iTerm2 preferences from my config folder on github!"


echo "You will need to MANUALLY run XCode to let it install its weird components"
echo "maybe something like this: sudo /Applications/Xcode.app/Contents/Developer/usr/bin/xcodebuild  -license accept"
#sudo installer -pkg /Applications/Xcode-beta.app/Contents/Resources/Packages/MobileDevice.pkg -target /
#sudo installer -pkg /Applications/Xcode-beta.app/Contents/Resources/Packages/MobileDeviceDevelopment.pkg -target /

echo "[DONE!] But remember to manually install the software shown above."

if [[ 1 == 2 ]]; then
    echo "NOTE: currently this is never run"
    echo "Installing IGV..."
    cd ~/bin
    git clone https://github.com/igvteam/igv.git
    cd igv/
fi

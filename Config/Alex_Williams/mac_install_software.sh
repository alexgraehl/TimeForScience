#!/usr/bin/env bash

cd ~/

if [[ -z "$TIME_FOR_SCIENCE_DIR" ]]; then echo "[ERROR: a required variable is NOT set!] You need to set the TIME_FOR_SCIENCE_DIR to the location of the TimeForScience directory (probably in your home directory). You could do this by following the instructions at https://github.com/alexgraehl/TimeForScience."; exit 1; else echo "[OK] The TIME_FOR_SCIENCE_DIR environment variable is '$TIME_FOR_SCIENCE_DIR'"; fi

if [[ ! -d "$TIME_FOR_SCIENCE_DIR" ]]; then echo "[ERROR: the TIME_FOR_SCIENCE_DIR does not exist (or is not a directory)] The TIME_FOR_SCIENCE_DIR environment variable must point to a valid checked-out copy of the TimeForScience repository. It appears to not be a directory. The directory we are looking for (which you should verify is correct) is: '$TIME_FOR_SCIENCE_DIR'"; exit 1; else echo "[OK] Successfully verified the existence of the directory '$TIME_FOR_SCIENCE_DIR'"; fi

if [[ ! -e "$TIME_FOR_SCIENCE_DIR/README.md" ]]; then echo "[ERROR: the TIME_FOR_SCIENCE_DIR appears NOT to contain a 'README.md', which implies that it is not a valid checked-out copy of the TimeForScience repository."; exit 1; else echo "[OK] The $TIME_FOR_SCIENCE_DIR/README.md file was successfully found."; fi

cd "$HOME"


SHOULD_HOMEBREW=1

# ==============================================================================
#            Create symlinks to standard config files
# ==============================================================================
if [[ $(dirname $TIME_FOR_SCIENCE_DIR) == "$HOME" ]]; then echo "use relative links"; else echo "do not use relative links--apparently the TimeForScience directory is NOT located immediately in your home directory2. Exiting here, as we cannot deal with this."; exit 1; fi
cd $HOME
HOMEDIR_CONFIG_FILES_TO_SYMLINK=(.aliases .bash_profile .bashrc .emacs .inputrc .screenrc .tmux.conf)
for FILE in ${HOMEDIR_CONFIG_FILES_TO_SYMLINK[@]}; do
	TIME_FOR_SCIENCE_BASENAME=$(basename $TIME_FOR_SCIENCE_DIR)
	CFG="./$TIME_FOR_SCIENCE_BASENAME/Config/Alex_Williams/$FILE"
	echo "CFG="$CFG
	echo "FILE="$FILE
	if [[ ! -e $CFG ]]; then echo "[ERROR] The config file that we expected to exist, specifically, '$CFG', did not exist."; exit 1; else echo "[OK] Linking $CFG to the home directory (via relative path)"; fi
	ln -s $CFG ./
done
# ==============================================================================

if [[ -e "$HOME/.ssh/config" ]]; then
    echo "[SKIPPING] The ssh config file at .ssh/config ALREADY exists, so we are not overwriting it."
else
    cp TimeForScience/Config/Alex_Williams/ssh-slash-config $HOME/.ssh/config
    chmod -R go-w   ~/.ssh/
    chmod    go-rwx ~/.ssh/id_rsa
fi

echo "install from the Mac App store: Slack Divvy Pixelmator"
echo "install iTerm2 from ____________"
echo "Sync the iTerm2 preferences from my config folder on github"
echo "install OmniGraffle (Paid!) from ____________"
echo "install Transmit (Paid!) (File transfer) from ____________"
echo "install Evernote from ____________"
echo "install Sourcetree (git/mercurial) from ____________"
echo "install Hermes (Pandora) from ____________"
echo "install Clipy from https://clipy-app.com/ (successor to ClipMenu)"

echo "Install XCode"
echo "Run XCode to let it install its weird components"
echo "maybe something like this: sudo /Applications/Xcode.app/Contents/Developer/usr/bin/xcodebuild  -license accept"
#sudo installer -pkg /Applications/Xcode-beta.app/Contents/Resources/Packages/MobileDevice.pkg -target /
#sudo installer -pkg /Applications/Xcode-beta.app/Contents/Resources/Packages/MobileDeviceDevelopment.pkg -target /

# ==============================================================================
#                         Homebrew
# ==============================================================================


if [[ 1 == $SHOULD_HOMEBREW ]]; then
    brew tap homebrew/science
    brew tap homebrew/versions
    brew update
    brew cask install java # required for gradle and more
    brew install bamutil bash bcftools bedtools blast boost cairo dialog emacs ess fastqc ffmpeg fontconfig fqzcomp freetype gcc gdbm gettext git glib gmp gnutls gradle gsl hdf5 htslib imagemagick isl jpeg lame libevent libffi libmpc libpng libtasn1 libtiff libtool libvo-aacenc mercurial mono mpfr nettle openssl pcre pixman pkg-config poretools pv qemu readline samtools sqlite szip tmux watch wget wxmac wxpython x264 xvid xz
    brew install docker
    brew install Caskroom/cask/gpgtools
    brew install homebrew/dupes/unzip htop-osx
fi

# ==============================================================================

AGW_BASH_MAJOR_VERSION=$(echo "$BASH_VERSION" | perl -pe "s/(\d+).*/\1/; chomp;")

if (( $AGW_BASH_MAJOR_VERSION >= 4 )); then
	echo "[OK] globstar should be supported on bash version $AGW_BASH_MAJOR_VERSION"
else
	echo "[GOTTA FIX] globstar is NOT supported on bash version $AGW_BASH_MAJOR_VERSION, so let's add a new bash"
fi

if [[ 'some bash thing' == "yep" ]]; then
    NEW_BASH_LOCATION=/usr/local/bin/bash
    chsh -s $NEW_BASH_LOCATION
    grep $NEW_BASH_LOCATION /etc/shells
    echo "$NEW_BASH_LOCATION" >> /etc/shells
    #sudo bash -c 'echo "/usr/local/bin/bash" >> /etc/shells'
fi

# ==============================================================================
echo "Making a bunch of uninteresting files/folders hidden in the default Finder view..."
FILES_TO_INVISIFY=("$HOME/Movies" "$HOME/Music" "$HOME/Public" "$HOME/TimeForScience" "$HOME/bin")
for FILE in ${FILES_TO_INVISIFY[@]}; do
    chflags -h hidden $FILE
done
# ==============================================================================

# ==============================================================================
echo "A bunch of things that should be installed via pip (2 vs 3?)"
PIP2_PACKAGES=("awscli", "pylint", "matplotlib", "numpy")
for PACKAGE in ${PIP2_PACKAGES[@]}; do
    print "Installing the package $PACKAGE via pip..."
    pip2 install $PACKAGE
done
# ==============================================================================

#echo "Install the amazon command line tools (via python's pip)"
#pip install awscli
#echo "[Installed] the 'aws' tool"

#echo "Installing pylint for Python 2.x"
#pip2 install pylint

echo "[DONE]"



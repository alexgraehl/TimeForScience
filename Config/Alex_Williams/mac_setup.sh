#!/usr/bin/env bash

cd ~/

if [[ -z "$TIME_FOR_SCIENCE_DIR" ]]; then echo "[ERROR: a required variable is NOT set!] You need to set the TIME_FOR_SCIENCE_DIR to the location of the TimeForScience directory (probably in your home directory). You could do this by following the instructions at https://github.com/alexgraehl/TimeForScience."; exit 1; else echo "[OK] The TIME_FOR_SCIENCE_DIR environment variable is '$TIME_FOR_SCIENCE_DIR'"; fi

if [[ ! -d "$TIME_FOR_SCIENCE_DIR" ]]; then echo "[ERROR: the TIME_FOR_SCIENCE_DIR does not exist (or is not a directory)] The TIME_FOR_SCIENCE_DIR environment variable must point to a valid checked-out copy of the TimeForScience repository. It appears to not be a directory. The directory we are looking for (which you should verify is correct) is: '$TIME_FOR_SCIENCE_DIR'"; exit 1; else echo "[OK] Successfully verified the existence of the directory '$TIME_FOR_SCIENCE_DIR'"; fi

if [[ ! -e "$TIME_FOR_SCIENCE_DIR/README.md" ]]; then echo "[ERROR: the TIME_FOR_SCIENCE_DIR appears NOT to contain a 'README.md', which implies that it is not a valid checked-out copy of the TimeForScience repository."; exit 1; else echo "[OK] The $TIME_FOR_SCIENCE_DIR/README.md file was successfully found."; fi

cd "$HOME"


if [[ $(dirname $TIME_FOR_SCIENCE_DIR) == "$HOME" ]]; then echo "use relative links"; else echo "do not use relative links. Exiting here, as we cannot deal with this."; exit 1; fi


cd $HOME

HOMEDIR_CONFIG_FILES_TO_SYMLINK=(.aliases .bash_profile .bashrc .emacs .inputrc .screenrc .tmux.conf)
for XFN in ${HOMEDIR_CONFIG_FILES_TO_SYMLINK[@]}; do
	TIME_FOR_SCIENCE_BASENAME=$(basename $TIME_FOR_SCIENCE_DIR)
	CFG="./$TIME_FOR_SCIENCE_BASENAME/Config/Alex_Williams/$XFN"
	echo "CFG="$CFG
	echo "XFN="$XFN
	if [[ ! -e $CFG ]]; then echo "[ERROR] The config file that we expected to exist, specifically, '$CFG', did not exist."; exit 1; else echo "[OK] Linking $CFG to the home directory (via relative path)"; fi
	ln -s $CFG ./
done


cp TimeForScience/Config/Alex_Williams/ssh-slash-config .ssh/config
chmod -R go-w   ~/.ssh/
chmod    go-rwx ~/.ssh/id_rsa


echo "wget some common applications if they aren't present"
echo "get evernote"
echo "get omnigraffle"
echo "tell the user to download slack and divvy from the app store"

echo "symlink the .whateverrc files"

echo "get iterm2"
echo "symlink the settings for iTerm2"

echo "brew install whatever"

echo "installing new bash since old one does not use globstar"
echo "brew install bash git emacs vim imagemagick wget htop"

echo "make some files invisible"

if [[ "homebrew" == "yep do it" ]]; then
    brew tap homebrew/science
    brew tap homebrew/versions
    brew install bamutil bash bedtools blast boost cairo clustal-w dialog emacs ess ffmpeg fontconfig fqzcomp freetype gcc gdbm gettext git glib gmp gnutls gradle gsl hdf5 htslib imagemagick isl jpeg lame libevent libffi libmpc libpng libtasn1 libtiff libtool libvo-aacenc mono mpfr nettle openssl pcre pixman pkg-config poretools pv qemu readline samtools sqlite szip tmux watch wget wxmac wxpython x264 xvid xz
fi

AGW_BASH_MAJOR_VERSION=$(echo "$BASH_VERSION" | perl -pe "s/(\d+).*/\1/; chomp;")

if (( $AGW_BASH_MAJOR_VERSION >= 4 ]]; then
	echo "[OK] globstar should be supported on bash version $AGW_BASH_MAJOR_VERSION"
else
	echo "[GOTTA FIX] globstar is NOT supported on bash version $AGW_BASH_MAJOR_VERSION, so let's add a new bash"
fi

NEW_BASH_LOCATION=/usr/local/bin/bash
chsh -s $NEW_BASH_LOCATION
grep $NEW_BASH_LOCATION /etc/shells
echo "$NEW_BASH_LOCATION" >> /etc/shells
#sudo bash -c 'echo "/usr/local/bin/bash" >> /etc/shells'



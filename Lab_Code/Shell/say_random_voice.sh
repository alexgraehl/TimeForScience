#!/bin/bash
### say_random_voice.sh: Run the Mac 'say' script with a randomly-chosen voice
###
### Usage:
###   say_random_voice.sh 'Hello, I am a computer'
###
### Options:
###   (none)         (No options yet)
###   -h or --help:  Show this help message

set -u
set -o pipefail
# Don't set -e, which EXITS if any return code is non-zero!

# ~~~~~~~~~ Method for printing usage based on the '###' at very top of file ~~~~
if [[ $# == 0 ]] || [[ "$1" == "-h" ]] || [[ "$1" == "--help" ]]; then
    # Cool way to print help messages from https://samizdat.dev/help-message-for-shell-scripts/
    # Any line ANYWHERE in this file that starts with '###' will be printed.
    awk -F'### ' '/^###/ { print $2 }' "$0"
    exit 1
fi
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

echoerr() { echo -e "$@" 1>&2; }


nonempty_file_exists() { # returns 1, but does not exit with an error
    # How to use it:   if nonempty_file_exists "myfile" ; then echo "Exists!"; fi
    if [[ "$#" != 1 ]]; then echoerr "[:ERR:]] requires EXACTLY ONE file/directory name!"; exit 1; fi # Bad arguments to the function
    if [[ ! -f "$1" ]]; then echoerr "[:HEY:] File does not exist: '$1'"; return 1; fi
    if [[ ! -s "$1" ]]; then echoerr "[:HEY:] File exists, but is EMPTY: '$1'"; return 1; fi
    return 0 # 0 is the 'good' status
}

require_nonempty_file() {
    if [[ "$#" != 1 ]]; then echoerr "[:ERR:]] requires EXACTLY ONE file/directory name!"; exit 1; fi # Bad arguments to the function
    if ! nonempty_file_exists "$1"; then exit 1; fi
    return 0 # 0 is the 'good' status
}

require_executable() {
    if [[ "$#" != 1 ]]; then echoerr "[:ERR:]] requires EXACTLY ONE file/directory name!"; exit 1; fi
    if [[ ! -f "$1" ]]; then echoerr "[:ERR:] Cannot find required executable: '$1'"; exit 1; fi
    if [[ ! -x "$1" ]]; then echoerr "[:ERR:] Program was found, but is NOT EXECUTABLE at: '$1'"; exit 1; fi
    if [[ ! -s "$1" ]]; then echoerr "[:ERR:] The executable exists, but is EMPTY, which is probably a mistake!: '$1'"; exit 1; fi
    return 0 # 0 is the 'good' status
}

require_executable $(which 'say')

say_this="$1"  # the thing to say

voice_str=$(say -v '?' | cut -d ' ' -f 1)
IFS=$'\n'  read -r -d '' -a voice_arr <<< "$voice_str"

#for x in "${voice_arr[@]}"
#do
#    echo "Voice option --> $v"
#done

RANDOM=$$$(date +%s)
# echo $RANDOM
rand_idx=$(expr $RANDOM % ${#voice_arr[@]})  # random valid indx in the voice_arr (% is modulo)

rand_voice="${voice_arr[$rand_idx]}"

#echo "(The voice is ${rand_voice})"
(set -x ; say -v "${rand_voice}" "${say_this}")


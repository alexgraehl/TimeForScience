#!/bin/bash
### Test_script.sh: DESCRIPTION SHOULD GO HERE
###
### Usage:
###   HOW_TO_USE_IT
###
### Options:
###   THING 1        A thing
###   THING 2        Another thing
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

#if ! type "command_name_here" > /dev/null; then
#    echo "Executable / command is missing"
#fi

# This does not work
#is_numeric() {
#    if [[ "$#" != 1 ]]; then echoerr "[:ERR:]] is_numeric requires EXACTLY ONE ARGUMENT!"; exit 1; fi
#    if [[ $1 =~ ^[+-]?[0-9]+([.][0-9]+)?$ ]] ; then
#	return 0 # is a number (0 is 'success' in bash)
#    else
#	return 1 # not a number (1 is 'not success' in bash)
#    fi
#    # Remember that you can't capture the return value with x=$(is_numeric ...): you must check $?
#}

#if [[ -z $(is_numeric $1) ]]; then
#    echo "Is numeric"
#else
#    echo "Not numeric"
#fi

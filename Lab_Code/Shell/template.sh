#!/bin/bash
set -u
set -o pipefail
# Don't set -e, which EXITS if any return code is non-zero!

echoerr() { echo -e "$@" 1>&2; }


nonempty_file_exists() { # returns 1, but does not exit with an error
    if [[ "$#" != 1 ]]; then echoerr "[:ERR:]] requires EXACTLY ONE file/directory name!"; exit 1; fi # Bad arguments to the function
    if [[ ! -f "$1" ]]; then echoerr "[:HEY:] File does not exist: '$1'"; return 1; fi
    if [[ ! -s "$1" ]]; then echoerr "[:HEY:] File exists, but is EMPTY: '$1'"; return 1; fi
    return 0 # 0 is the 'good' status
}

require_nonempty_file() {
    if [[ "$#" != 1 ]]; then echoerr "[:ERR:]] requires EXACTLY ONE file/directory name!"; exit 1; fi # Bad arguments to the function
    if [[ ! $(nonempty_file_exists "$1") ]]; then exit 1; fi
    return 0 # 0 is the 'good' status
}

require_executable() {
    if [[ "$#" != 1 ]]; then echoerr "[:ERR:]] requires EXACTLY ONE file/directory name!"; exit 1; fi
    if [[ ! -f "$1" ]]; then echoerr "[:ERR:] Cannot find required executable: '$1'"; exit 1; fi
    if [[ ! -x "$1" ]]; then echoerr "[:ERR:] Program was found, but is NOT EXECUTABLE at: '$1'"; exit 1; fi
    if [[ ! -s "$1" ]]; then echoerr "[:ERR:] The executable exists, but is EMPTY, which is probably a mistake!: '$1'"; exit 1; fi
    return 0 # 0 is the 'good' status
}

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

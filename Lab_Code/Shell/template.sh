#!/bin/bash
set -u
set -o pipefail
# Don't set -e, which EXITS if any return code is non-zero!

echoerr() { echo "$@" 1>&2; }

we_need_file() {
    if [[ "$#" != 1 ]] ; then echo -e "[:ERR:]] requires EXACTLY ONE file/directory name!"; exit 1; fi
    if [[ ! -f "$1" ]]; then echo "[:ERR:] Cannot find required file: '$1'"; exit 1; fi
    if [[ ! -s "$1" ]]; then echo "[:ERR:] The file exists, but is EMPTY, which is probably a mistake!: '$1'"; exit 1; fi
    return 0
}

need_executable() {
    if [[ "$#" != 1 ]] ; then echo -e "[:ERR:]] requires EXACTLY ONE file/directory name!"; exit 1; fi
    if [[ ! -f "$1" ]]; then echo "[:ERR:] Cannot find required executable: '$1'"; exit 1; fi
    if [[ ! -x "$1" ]]; then echo "[:ERR:] Program was found, but is NOT EXECUTABLE at: '$1'"; exit 1; fi
    if [[ ! -s "$1" ]]; then echo "[:ERR:] The file exists, but is EMPTY, which is probably a mistake!: '$1'"; exit 1; fi
    return 0
}

# This does not work
#is_numeric() {
#    if [[ $1 =~ ^[+-]?[0-9]+([.][0-9]+)?$ ]] ; then
#	return 0 # is a number (0 means 'success' in bash)
#    else
#	return 1 # not a number (1 means 'not success' in bash)
#    fi
#    # Remember that you can't capture the return value with x=$(is_numeric ...): you must check $?
#}

#if [[ -z $(is_numeric $1) ]]; then
#    echo "Is numeric"
#else
#    echo "Not numeric"
#fi

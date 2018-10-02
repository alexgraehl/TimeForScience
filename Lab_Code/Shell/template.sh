#!/bin/bash
set -ueo pipefail

we_need() {
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

echoerr() { echo "$@" 1>&2; }


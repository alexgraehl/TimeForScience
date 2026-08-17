#!/usr/bin/env bash
# By Alex, last edited in 2026.

### say_random_voice.sh: Runs the Mac built-in 'say' utility with a randomly-chosen
###                      voice (optionally in a specified language/locale). Only works on the Mac.
###
### Usage:
###   say_random_voice.sh              "Hello, I am a computer"         # Implicitly English voices only.
###   say_random_voice.sh --lang=en_gb "Hello, I am a British computer" # A British-locale English voice
###   say_random_voice.sh --lang=es    "Hola, soy computadora"          # Any locale of Spanish (es_...)
###   say_random_voice.sh --lang=any   "Any voice is allowed here"
###   say_random_voice.sh --lang=zh    "你好"                           # (Any locale of Chinese (zh_...)
###   say_random_voice.sh --lang=en --list                             # Just list all English voices
###   say_random_voice.sh --lang=en_us --list                          # Just list all U.S. English voices
###
### Options:
###   --lang=PATT: Only use voices for this LOCALE. Case-insensitive. Default: "en_" (any English).
###                Format is xx_YY (language_COUNTRY: run say -v '?' on the command line to see the full list)
###                Examples: "--lang=fr" matches any French, "--lang=en_gb" matches British English only.
###                Set to --lang="" or --lang=ANY to allow ALL languages.
###   --list:      List voices matching the '--lang' filter and then exit.
###   --help (-h): Show this help message

set -u
set -o pipefail
# Don't set -e, which EXITS if any return code is non-zero!

# ~~~~~~~~~ Method for printing usage based on the '###' at very top of file ~~~~
# Detect 1) no input arguments, 2) '-h' or 3) '--help' and print all lines beginning with '###'.
if [[ $# == 0 ]] || [[ "$1" == "-h" ]] || [[ "$1" == "--help" ]]; then
    # Cool way to print help messages from https://samizdat.dev/help-message-for-shell-scripts/
    # Any line ANYWHERE in this file that starts with '###' will be printed.
    awk -F'### ' '/^###/ { print $2 }' "$0"
    exit 1
fi
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Handle the input args:
lang="en_"   # Used to select which langages to support.
onlylist=0   # Effectively a boolean (0 or 1). Only
say_this=""
for arg in "$@"; do
    case "$arg" in
        --lang=*)  lang="${arg#--lang=}" ;;
        --list)    onlylist=1 ;;
        *)         say_this="${say_this:+$say_this }$arg" ;;  # APPEND unquoed multi-word input
    esac
done

if [[ "$(echo "$lang" | tr '[:upper:]' '[:lower:]')" == "any" ]]; then
    lang="" # Treat "--lang=any" the same as "--lang=" (no filtering at all, ALL voices matched)
fi

if [[ "$onlylist" == 1 ]] && [[ -n "$say_this" ]]; then
    echoerr "[:ERR:] You can't '--list' AND say something at the same time. Try '--help'."
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

get_matching_voices() {
    # Args:
    #    $1: A search pattern for the locale to allow, matched case-insensitively (e.g. "^en_" for english-only)
    #        Expected to already start with "^" before being passed in, or something has gone wrong.
    # Returns:
    #    All the matching voices (if any).
    # Annoyingly, the "say -v '?'" output is not delimited in a normal fashion, so this function is
    # WAY more complicated than it seems like it needs to be. (It's space-delimited...ish, to
    # a relatively small maximum column width). Voice names can include spaces & parens, which causes woes.
    # Example output of say -v '?'). It's *NOT* tab or fixed-column-size delimited, even though it looks like it should be.
    #      Grandpa (English (UK)) en_GB    # Hello! My name is Grandpa.
    #      Grandpa (English (US)) en_US    # Hello! My name is Grandpa.
    #      Jester               en_US    # Hello! My name is Jester.
    #      Junior               en_US    # Hello! My name is Junior.
    #      Karen                en_AU    # Hello! My name is Karen.
    say -v '?' | awk -v \
        searchpat="$(echo "$1" | tr '[:upper:]' '[:lower:]')" '
        match($0, /[A-Za-z][A-Za-z]_[A-Za-z][A-Za-z][[:space:]]/) {
            locale = substr($0, RSTART, 5)
            voicename = substr($0, 1, RSTART - 1)
            gsub(/[[:space:]]+$/, "", voicename)
            if (tolower(locale) ~ searchpat) { print voicename }
        }'
}

# Read all the voice options (the first "column" of the say -v '?' output):
voice_arr=()
while IFS= read -r line; do
    voice_arr+=("$line")
done < <(get_matching_voices "^$lang")  # Note that we add '^" to the search pattern to REQUIRE matching from the start of the xx_YY locale string

# Abort if we didn't find a voice. This should never happen.
if [[ "${#voice_arr[@]}" == 0 ]]; then echoerr "[:ERR:] Nothing matched the locale filter: '${lang}'"; exit 1; fi

if [[ "$onlylist" == 1 ]]; then
    # List all voices that survived the filter, then exit normally.
    printf '%s\n' "${voice_arr[@]}"
    exit 0
fi

RANDOM=$$$(date +%s)                         # Random seed based on current date/time
rand_idx=$(expr $RANDOM % ${#voice_arr[@]})  # Random index in our voice array...
rand_voice="${voice_arr[$rand_idx]}"         # Finally get the random voice name (e.g. "Fred")

# The 'set -x' will echo the command to the command line, so we can see what voice was used, etc.)
(set -x ; say -v "${rand_voice}" "${say_this}")


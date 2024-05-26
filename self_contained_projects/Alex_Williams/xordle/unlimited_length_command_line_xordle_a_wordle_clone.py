#!/usr/bin/env python3
"""
A word-guessing game inspired by the game "Wordle" / "Quordle" / "Octordle", etc.
Similar to the word-guessing hacking mini-game in the Fallout series, which is also
similar to a word version of the game "Mastermind".

By default, after you make a guess:
 * A yellow letter means "this character is present in this word, but not at this position"
 * A green letter means "this character is present at this exact position."
 * Anything else means the letter does not appear in that word at all.

FYI: for testing, here's a one-liner to pull words out of the dictionary of a particular length:
  Length=8 --> cat /usr/share/dict/words | perl -ne 'print if /^.{8}$/' | less -S

Example command line arguments:
 SCRIPTNAME -n=4 -l=5 --daily --debug

Show a NON-INTERACTIVE demo of a set of guesses for six 12-letter words:
 SCRIPTNAME -n=6 -l=12 --daily --debug

Deterministically always pick the same words:
 SCRIPTNAME -n=5 -l=4 --seed=1234

Script by Alex Williams, Feb 2022.
"""

# pylint: disable=unnecessary-pass,useless-return,line-too-long
import argparse
import datetime
import random
import re
from ssl import HAS_TLSv1_1
import sys
from enum import Enum
from typing import Iterable


class Verdict(Enum):
    ABSENT = 0  # Letter is not in this word at all
    ELSEWHERE = 1  # Letter is present but not at this position
    CORRECT = 2  # Letter is present AND at this position
    FULL_WORD_CORRECT = 3  # Not a per-letter verdict. Indicates that it's part of a COMPLETE guess.


class Redaction(Enum):
    VISIBLE = 0  # Visibile, no redaction.
    ASCII = 1  # Ascii text
    EMOJI = 2  # Redact with emojis


RESET_COLOR = "\x1b[0m"  # Should reset the terminal to "normal" text color.
ANSWER_WORD_COLOR = "\x1b[0;31;40m"
BORDER_COLOR = "\x1b[0;34m"
GUESS_COUNT_WITHIN_LIMIT_COLOR = "\x1b[0;32;40m"
GUESS_COUNT_EXCEEDED_LIMIT_COLOR = "\x1b[0;31;40m"
ALREADY_GUESSED_LETTER_COLOR = "\x1b[0;34;40m"  # A dim color that isn't too eye-catching
NOT_YET_GUESSED_LETTER_COLOR = "\x1b[0;30;42m"

REDACT_CHAR_MAPPING: dict[Verdict, str] = {
    # Unix color formatting
    Verdict.ABSENT: " ",
    Verdict.ELSEWHERE: "-",
    Verdict.CORRECT: "=",
    Verdict.FULL_WORD_CORRECT: "#",
}

REDACT_EMOJI_MAPPING: dict[Verdict, str] = {
    # Emoji mapping
    Verdict.ABSENT: "⬛",
    Verdict.ELSEWHERE: "🟡",
    Verdict.CORRECT: "🟩",
    Verdict.FULL_WORD_CORRECT: "🔹",
}

COLOR_FORMAT_DEFAULT: dict[Verdict, str] = {
    # Unix color formatting
    Verdict.ABSENT: "\x1b[0;37;40m",  # White on black
    Verdict.ELSEWHERE: "\x1b[0;30;43m",  # Black text, yellow background
    Verdict.CORRECT: "\x1b[0;30;42m",  # Black text, green background
    Verdict.FULL_WORD_CORRECT: "\x1b[1;30;46m",  # Note the '1;' for "bold" colors
}

COLOR_FORMAT_COLORBLIND_1: dict[Verdict, str] = {
    # Unix color formatting
    Verdict.ABSENT: "\x1b[1;37;40m",  # White on black
    Verdict.ELSEWHERE: "\x1b[0;33;41m",  # Yellow on red
    Verdict.CORRECT: "\x1b[0;33;44m",  # Yellow on blue
    Verdict.FULL_WORD_CORRECT: "\x1b[1;30;46m",  # Note the '1;' for "bold" colors
}


def get_color_formatter(name: str = "") -> dict[Verdict, str]:
    if name in ("default", ""):
        return COLOR_FORMAT_DEFAULT
    elif name in ("colorblind", "colorblind_1", "colorblind1"):
        return COLOR_FORMAT_COLORBLIND_1


# The set of all possible words you'd find in a dictionary. Note that it's not a *python* dict.
WORDS: set[str] = set()

# The actually-picked words. One word per N "boards"
ANSWERS: list[str] = list()

UNSOLVED_VALUE = -1

# int showing which guess index each answer was found on, starting from index = 0 (first guess). -1 means "not correct yet"
SOLVED_AT: list[int] = list()

GUESSES: list[str] = list()  # your previous guesses


def is_guess_valid(guess: str, valid_words: Iterable[str], expected_len: int) -> bool:
    if len(guess) != expected_len:
        raise ValueError(
            f"Your guess must be a word of length {expected_len}, but your guess of `{guess}` was of length {len(guess)}."
        )
    if guess not in valid_words:
        raise ValueError(f"Guesses must be in the dictionary: your word (`{guess}`) was not.")
    return True


def manually_colorized(word: str, color: str) -> str:
    return f"""{color}{word}{RESET_COLOR}"""


def colorized(letter: str, v: Verdict) -> str:
    assert len(letter) == 1, "This should be a SINGLE character we are evaluating/colorizing."
    return f"""{get_color_formatter()[v]}{letter}{RESET_COLOR}"""


def evaluate_guess(guess_word: str, answer_word: str) -> list[Verdict]:
    """Evaluates a guessed WORD versus the answer WORD."""
    result: list[Verdict] = list()
    assert len(guess_word) == len(answer_word), "Programming error: length mismatch!"
    for i in range(len(guess_word)):
        if guess_word[i] == answer_word[i]:
            result.append(Verdict.CORRECT)  # Exactly correct
        elif guess_word[i] in answer_word:
            result.append(Verdict.ELSEWHERE)  # Letter is in the word, but not at this position
        else:
            result.append(Verdict.ABSENT)
    return result


def print_keyboard(guesses: Iterable[str]) -> None:
    char_between_keys: str = " "  # Horizontal space to print between keys (a delimiter)
    xxx: str = "EMPTY"  # Just used for layout / spacing
    h1: str = "HALF"  # Half width, for offsetting the keyboard by HALF as much as 'xxx' would. This is how we offset the keyboard to get the 'staggered' row look.
    h2: str = "HALF_RIGHT"  # This is used to make the right border line up, in case it doesn't with the default spacing settings.
    kb_layout: list[list[str]] = [  # This gets printed out essentially verbatim
        ["q", "w", "e", "r", "t", "y", "u", "i", "o", "p"],
        [xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx],  # blank line
        [h1, "a", "s", "d", "f", "g", "h", "j", "k", "l", h2],
        [xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx, xxx],  # blank line
        [xxx, "z", "x", "c", "v", "b", "n", "m", xxx, xxx],
    ]
    all_letters_guessed: set[str] = set("".join(guesses))  # Letters already guessed

    key_width = 1  # Each character is only one character wide (unless we had something like '[A]' instead of just 'A')

    l_border = "║   "
    r_border = "  ║"

    # How many monospaced characters wide is the keyboard?
    kb_width: int = (key_width + len(char_between_keys)) * len(kb_layout[0]) + len(l_border) + len(r_border)

    t_border = "╔" + "═" * (kb_width - 2) + "╗" + "\n"
    b_border = "╚" + "═" * (kb_width - 2) + "╝" + "\n"
    # per_row_delim = f"\n{l_border}\n"  # Two newlines means we get a blank line between rows
    sys.stdout.write(manually_colorized(t_border, BORDER_COLOR))
    row: list[str]
    for row in kb_layout:
        sys.stdout.write(manually_colorized(l_border, BORDER_COLOR))  # Left border for the "keyboard"
        k: str
        for k in row:
            if k == xxx:  # If it's the special "full spacer" string
                sys.stdout.write(f" {char_between_keys}")
                continue
            if k == h1:  # If it's the special "half spacer" string
                sys.stdout.write(f"{char_between_keys}")  # <-- less space than the full spacer
                continue
            if k == h2:  # If it's the special "half spacer" string
                sys.stdout.write(f"{char_between_keys}")  # <-- less space than the full spacer
                continue

            color: str = (
                ALREADY_GUESSED_LETTER_COLOR if (k in all_letters_guessed) else NOT_YET_GUESSED_LETTER_COLOR
            )
            sys.stdout.write(f"{color}{k}{RESET_COLOR}{char_between_keys}")  # Print this key
            pass  # End of this row.

        sys.stdout.write(
            f"{manually_colorized(r_border, BORDER_COLOR)}\n"
        )  # Just a newline, for the last keyboard row only.

        pass
    sys.stdout.write(manually_colorized(b_border, BORDER_COLOR))
    return


def print_board(
    guesses: Iterable[str],
    answers: Iterable[str],
    max_guesses: int = -1,  # <-- Not currently used
    redact_type: Redaction = Redaction.VISIBLE,
    toupper: bool = False,
    show_answers: bool = False,
    print_guess_num: bool = True,
    print_in_color: bool = True,
) -> None:
    n_boards: int = len(answers)

    ANSWER_HORIZONTAL_DELIM = "  "  # The delimiter between answers (columns)
    if n_boards == 0:
        raise Exception("Zero words? Questionable.")

    if print_guess_num:
        N_GUESS_DIGITS: int = 5
        LEFT_PAD_BEFORE_GUESS: str = " "
    else:
        N_GUESS_DIGITS: int = 0
        LEFT_PAD_BEFORE_GUESS: str = ""
        pass

    if show_answers:
        sys.stdout.write(" " * N_GUESS_DIGITS)  # Left-padding (guess numbers will go here later
        sys.stdout.write(LEFT_PAD_BEFORE_GUESS)
        for i, a in enumerate(answers):
            if i > 0:
                sys.stdout.write(ANSWER_HORIZONTAL_DELIM)
            sys.stdout.write(manually_colorized(a, ANSWER_WORD_COLOR))
        sys.stdout.write("\n")
        pass

    # Keep track of when we discover that each answer is correct.
    # Note that this is different from the global variable version, which doesn't actually store the position.
    guess: str
    for guess_num, guess in enumerate(guesses):
        # Each guess is a ROW
        padded_guess_num = "{0: >{width}}".format(guess_num + 1, width=N_GUESS_DIGITS)
        guess_color = GUESS_COUNT_WITHIN_LIMIT_COLOR
        if print_guess_num:
            sys.stdout.write(manually_colorized(padded_guess_num, color=guess_color))
            sys.stdout.write(LEFT_PAD_BEFORE_GUESS)
            pass

        answer: str
        for a_idx, answer in enumerate(answers):
            # Each answer is a COLUMN.
            verdict_this_word: list[Verdict] = evaluate_guess(guess, answer)
            assert len(guess) == len(verdict_this_word)

            if guess == answer:  # Got this entire word totally correct
                # print(f"Found this word --> {guess} was {answer} at index {a_idx}")
                SOLVED_AT[a_idx] = guess_num

            for letter, l_verdict in zip(guess, verdict_this_word):
                if toupper:
                    letter = letter.capitalize()

                if (SOLVED_AT[a_idx] != UNSOLVED_VALUE) and (SOLVED_AT[a_idx] < guess_num):
                    # AFTER we've solve the word, print it in a different color from that guess number onward
                    # (SO we print it in 'all solved' color for one line only)
                    l_verdict = Verdict.FULL_WORD_CORRECT
                    letter = " "  # Even if we aren't redacting, we still don't keep printing guesses for already-solved words.

                if redact_type == Redaction.ASCII:
                    # If we're in "redaction" mode, then don't print the guesses
                    letter = REDACT_CHAR_MAPPING[l_verdict]
                elif redact_type == Redaction.EMOJI:
                    letter = REDACT_EMOJI_MAPPING[l_verdict]
                    pass

                formatted_letter = colorized(letter=letter, v=l_verdict) if print_in_color else letter
                sys.stdout.write(formatted_letter)
                pass  # end of printing each letter

            sys.stdout.write(ANSWER_HORIZONTAL_DELIM)
            # end of printing all the 'boards' on this row

        sys.stdout.write("\n")  # <-- End of this row.
        pass
    return


def argErrorAndExit(msg="(No additional information given)"):
    raise SystemExit("[ERROR] in arguments to this script: " + msg)


def main():
    parser = argparse.ArgumentParser(
        description="%(prog)s: Word-guessing game.",
        epilog="""Example usage: %(prog)s -w 6 -l 5.""",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-d",
        "--dict",
        dest="dictionary",
        type=str,
        default="/usr/share/dict/words",
        metavar="FILE",
        help="The UNIX dictionary file. Default on the Mac is /usr/share/dict/words.",
    )
    parser.add_argument(
        "-n",
        "--n-boards",
        dest="n_boards",
        type=int,
        default=4,
        help="Number of `boards` you are playing on (one word per board). Anywhere from 1 to 16 is `reasonable`, but more is OK too.",
    )

    parser.add_argument("-l", "--letters", dest="word_len", type=int, default=5, help="Letters per word.")

    parser.add_argument(
        "--debug",
        dest="debug_force_input",
        action="store_true",
        help="'Automatically' runs through some guesses. Deterministic, so you may want to also use --seed=####.",
    )

    parser.add_argument(
        "--no-keyboard",
        dest="omit_keyboard",
        action="store_true",
        help="Don't print the keyboard.",
    )
    parser.add_argument(
        "-s",
        "--seed",
        dest="rand_seed",
        type=int,
        default=None,
        help="Random seed, for deterministic games. Otherwise",
    )

    parser.add_argument(
        "--daily-challenge",
        dest="is_daily_challenge",
        action="store_true",
        help="Should we use today's date (UTC) as the random seed? The default is just to make a random board.",
    )

    parser.add_argument(
        "--exclude-upper-case",
        dest="exclude_upper_case",
        action="store_false",
        help="Exclude words with ANY capital letters.",
    )

    parser.add_argument(
        "--exclude-non-english-alphabet",
        dest="exclude_non_english_alphabet",
        action="store_false",
        help="Exclude words with ANY letters that aren't in [A-Za-z]. E.g. no accent marks or apostrophes.",
    )

    parser.add_argument(
        "-g",
        "--num-guesses",
        dest="n_guesses",
        type=int,
        default=4,
        help="Number of guesses before you 'lose'. Doesn't really mean anything.",
    )
    parser.add_argument(
        "-q", "--quiet", dest="verbose", action="store_false", help="don't print status messages to stdout"
    )
    parser.add_argument(
        "remainder", nargs=argparse.REMAINDER
    )  # get the REMAINING un-parsed arguments (for example, a bunch of filenames)

    args = parser.parse_args()

    if not len(args.remainder) == 0:
        argErrorAndExit("Looks like you have an extra command line argument?")

    if args.n_boards == 0:
        argErrorAndExit("You can't have ZERO words to guess. Try 1 or more.")

    if args.n_boards > 600:
        argErrorAndExit("The terminal can't currently handle this many boards. Try fewer.")

    try:
        with open(args.dictionary, str("r")) as dfile:
            for linenum, line in enumerate(dfile):
                word: str = line.strip()  # Strip newline
                if len(word) != args.word_len:
                    continue  # Word is not the desired length
                if args.exclude_upper_case and re.search(r"[A-Z]", word):
                    # print(f'''Line {linenum+1}: skipped word with upper case: {word}''')
                    continue  # Skip upper-case word
                if args.exclude_non_english_alphabet and re.search(r"[^A-Za-z]", word):
                    # Skip this word with a non-A-through-Z letter in it
                    # print(f'''Line {linenum+1}: skipped word non-"English A-Z" char: {word}''')
                    continue  # Skip word with non-"A-Z" symbol in it.
                WORDS.add(word)
            pass
    except Exception as e:
        argErrorAndExit(f"Failed to read from alleged dictionary file `{args.dictionary}: {e}")

    # print(f"Read this many words of length {args.word_len}: {len(WORDS)}")

    if len(WORDS) < args.n_boards:
        argErrorAndExit(
            f"""Error: you requsted {args.n_boards} words,
         but the supplied dictionary fileonly has {len(WORDS)} words of a suitable length.
         Double check the dictionary (at `{args.dictionary}`), or request fewer words."""
        )

    if args.is_daily_challenge and (args.rand_seed is not None):
        argErrorAndExit("You cannot specify BOTH `--daily` and `--seed` simultaneously.")

    seed: int | None
    if args.is_daily_challenge:
        # The random seed is "today's" (UTC timezone) days since Jan 1, 2000.
        seed = (datetime.datetime.utcnow().date() - datetime.date(2000, 1, 1)).days
    elif args.rand_seed is None:
        seed = None  # Randomize the seed based on the current time
        print(f"This is a random run, with seed={seed}")
    else:
        seed = args.rand_seed
        print(f"This is a NON-RANDOM run, with seed={args.rand_seed}")

    random.seed(seed)

    # Sorted list avoids nondeterminism. set->list is nondeterministic.
    ANSWERS.extend(random.sample(sorted(list(WORDS)), args.n_boards))
    SOLVED_AT.extend([UNSOLVED_VALUE] * len(ANSWERS))

    print(f"There are {args.n_boards} word(s) to be solved:")
    # print("FYI, the answers are: ", " ".join(ANSWERS))

    # For testing purposes, '--debug' will force a specific set of automatic TEST_GUESSES.
    TEST_GUESSES: list[str] = list()
    if args.debug_force_input:
        print("DEBUGGING is on: we're generating DETERMINISTIC test input.")
        print("DEBUG guesses are: N random gueses, 1 correct guess, repeat…")
        answers_not_yet_guessed: list[str] = random.sample(ANSWERS, len(ANSWERS))
        sorted_wordlist = sorted(list(WORDS))  # Set -> list is apparently nondeterministic.
        while len(answers_not_yet_guessed) > 0:
            n_rand_guesses_between_correct: int = 2
            TEST_GUESSES.extend(random.sample(population=sorted_wordlist, k=n_rand_guesses_between_correct))
            TEST_GUESSES.append(answers_not_yet_guessed.pop())
            pass
        pass

    finished: bool = False
    while not finished:
        if TEST_GUESSES:
            guess = TEST_GUESSES.pop(0)
        else:
            guess = input("🟢 Guess a word: ")

        try:
            _ = is_guess_valid(guess, valid_words=WORDS, expected_len=args.word_len)
        except ValueError as e:
            sys.stderr.write(f"🔴 {e}\n")
            sys.stderr.write("🔴 Try another word...\n")
            continue

        GUESSES.append(guess)
        print(f"Guess #{len(GUESSES)}: “{manually_colorized(guess, color=GUESS_COUNT_WITHIN_LIMIT_COLOR)}”")
        print_board(guesses=GUESSES, answers=ANSWERS, show_answers=False)
        if not args.omit_keyboard:
            print_keyboard(guesses=GUESSES)
        print("\n")
        if all([x != UNSOLVED_VALUE for x in SOLVED_AT]):
            print(f"🟢 You got all {len(ANSWERS)} word(s) in this many guesses: {len(GUESSES)}")
            finished = True
        pass

    print()  # newline
    print_board(guesses=GUESSES, answers=ANSWERS, show_answers=True)
    print("\n\nSpoiler-free text version...\n")
    print_board(guesses=GUESSES, answers=ANSWERS, redact_type=Redaction.ASCII, show_answers=False)
    print("\n\nSpoiler-free emoji version:\n")
    print_board(
        guesses=GUESSES,
        answers=ANSWERS,
        redact_type=Redaction.EMOJI,
        print_guess_num=False,
        show_answers=False,
        print_in_color=False,
    )
    return  # end of 'main'


if __name__ == "__main__":
    main()
    pass

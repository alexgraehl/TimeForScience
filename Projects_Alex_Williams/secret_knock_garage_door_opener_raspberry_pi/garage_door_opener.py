#!/usr/bin/env python3

# Example:
#  python garage_door_opener.py --knock=15 --code=1234 --scale=3 --verbose

# Venv cheat sheet:
# 1. One time only:
#    python3 -m venv venv  # The second 'venv' is the directory name!
# 2. To install the prereqs:
#    1. Install the `pyaudio` prerequisites, which can be done with:
#       A) If you're on a Mac (in 2026):     brew install portaudio  ### Or install `portaudio` some other way
#       B) On a Raspberry Pi (Linux, 2026):  sudo apt install libasound2-dev portaudio19-dev libportaudio2 libportaudiocpp0 python3-pyaudio
#                                 and then:  sudo apt install python3-gpiozero
#    2. Then activate `venv`:  source venv/bin/activate  # Must be in the same directory: note that `venv` is the directory name here
#    3. ...and then install: pip install pyaudio numpy gpiozero
#       (Note that the pip installation only takes effect in that specific directory with the `venv/` subdirectory in it.
# 3. Then, every time you want to run this script:
#    source venv/bin/activate   # (To pick up `pyaudio` and `numpy`)

import sys
import time
import argparse

import gpiozero
import pyaudio
import numpy as np

class color:
    g = '\033[92m' # GREEN
    y = '\033[93m' # YELLOW
    r = '\033[91m' # RED
    reset = '\033[0m'

BUTTON_PRESS_TIME_SEC: int = 1  # How long to press the garage door button
PULL_PIN_LOW_TO_ACTIVATE: bool = True  # Hard-coded. Right now we hard-code that we set the pin to LOW to activate it. Could be refactored.

FORMAT                    = pyaudio.paFloat32
AUDIO_DTYPE               = np.float32
CHANNELS: int             = 1        # (Mono)
RATE: int                 = 44100    # 44.1kHz sampling rate
CHECKS_PER_SECOND: int    = 10  # How many times to check for the current audio volume every second. Can be pretty low.
CHUNK: int                = int(RATE / (1+CHECKS_PER_SECOND))   # Amount of data to read in one "chunk". Might be slightly off (not sure if this division is totally correct)
MAX_VOLUME_BAR_WIDTH: int = 50     # Max number of characters for the volume bar. If this is larger than the terminal width, then the bar won't overwrite itself (and the terminal will scroll)


# These are chosen somewhat empirically and probably should be command line args.
MIN_KNOCK_LENGTH: float      = 0.05  # Absolute minimum number of seconds to wait before counting another knock. If you set it too HIGH, then two knocks will be counted as one.
MIN_TIME_BETWEEN_SETS: float = 0.5   # Minimum required gap (in sec.) between sets of knocks. E.g. {X,X,X} (gap) {X,X} would be 3, 2 (and not "5"). If you set it too low, it will be hard to reliably create "sets" of knocks instead of a bunch of single knocks.
FULL_RESET_TIMEOUT: int      = 5     # Reset the code after this many seconds of silence. Possibly not actually needed.

FINAL_SILENCE_REQUIRED: int = 1  # Register a code only after silence at least THIS long

IS_DRY_RUN: int = False
VERBOSE: bool = False

"""Converts a percent (number) into a color. Loud = red, medium = yellow, quiet = green."""
def color_for_pct(p: int | float):
    if p <= 5:  # Arbitrary cutoffs here
        return color.g
    if p <= 20:
        return color.y
    return color.r

BUFFER_SIZE: int = 100  # Keep track of the last N volumes. Only matters when `--verbose` is on.
RECALCULATE_MAX_EVERY: int = 10 # Recalculate the max voulume every ~N runs through the loop, just to be every so slightly less wasteful. Only matters when `--verbose` is on.

def int_percent(x):
    try:
        x = int(x)
    except ValueError:
        raise argparse.ArgumentTypeError(f"{x} doesn't seem to be castable to an int!")
    if 1 <= x <= 100:
        return x
    raise argparse.ArgumentTypeError(f"{x} is not a valid value for this argument: it must be in the range [1, 100]")

def str_with_digits_from_1_to_9_only(x: str):
    if not x:
        raise argparse.ArgumentTypeError(f"Hey! You MUST specify a value for this argument: it cannot be `None` or the empty string.")
    if str(x).isdigit() and ('0' not in str(x)):
        return x
    raise argparse.ArgumentTypeError(f"'{x}' must consist ONLY of the digits 1-9 and nothing else.")

# Check to see if secret_code is a suffix of entire_knock_sequence.
# E.g. if secret_code = [1,2] and entire_knock_sequence = [9,9,1,4,2,1,2], then it's true (ends in '1,2')
def secret_code_found_in_suffix(secret_code: list[int], entire_knock_sequence: list[int]):
    if not entire_knock_sequence or len(entire_knock_sequence) < len(secret_code):
        return False  # Not long enough!
    return entire_knock_sequence[-len(secret_code):] == secret_code  # Check just the suffix


def handle_audio(rescale_volume: float, knock_pct_threshold: float, secret_knock_code: list[int], pin_obj: gpiozero.OutputDevice):
    assert knock_pct_threshold >= 0 and knock_pct_threshold <= 100, "knock_pct_threshold must be in [0, 100]"
    assert isinstance(secret_knock_code, list) and all(isinstance(x, int) for x in secret_knock_code), "secret_knock_code should be an array of ints"

    last_n_volumes = np.zeros(BUFFER_SIZE, dtype=np.int8)
    buffer_index: int = 0
    max_recent_volume: int = 0

    currently_in_a_knock: bool = False
    last_knock_time = 0
    current_knock_count = 0
    current_sequence = []

    max_knock_sets_to_record_before_resetting: int = len(secret_knock_code) * 20  # Just make sure we don't keep track of knock sets forever if there's a woodpecker or something

    last_code_we_printed = []  # Just remember the last code we updated the user about, so we don't keep spamming the logs

    audio_stream = init_audio_stream()
    while True:
        try:
            data = audio_stream.read(CHUNK, exception_on_overflow=False)
        except KeyboardInterrupt:
            raise  # Abort if there's a keyboard interrupt
        except: # Tolerate anything ELSE...
            # If it fails due to some kind of technical issue with the audio system, try to re-initialize the stream.
            # May or may not actually help, but allegedly this can occur intermittently.
            print("Re-initializing the audio stream...")
            audio_stream = init_audio_stream()
            continue  # Back to the top of the loop
                
        # Volume level is `root mean squared` (rms), NOT just the average of the sound wave. Note that the (naive) mean can be negative.
        audio_data = np.frombuffer(data, dtype=AUDIO_DTYPE) # (Make a numpy array
        rms = np.sqrt(np.mean(audio_data.astype(float)**2)) * rescale_volume  # rms should be in [0, 1]
        # Convert to percent (clamped at 100%). This is a holdover from when I was doing everything in percent space.
        pct: float = 100.0 * min(1, rms)
        
        now: float = time.time()
        time_since_last_knock: float = now - last_knock_time
        if (pct > knock_pct_threshold) and not currently_in_a_knock:
            if time_since_last_knock < MIN_KNOCK_LENGTH:
                pass # We're still in the current knock! Don't double-count it.
            else:
                currently_in_a_knock = True
                current_knock_count += 1
                print(f"KNOCK: {current_knock_count}")
                last_knock_time = now

        # See if the current knock is now "done" being loud (volume < 50% of threshold AND it hasn't been a ludicrously short amount of time) (and prepare us to count the next one)
        if (pct < knock_pct_threshold * 0.5 and time_since_last_knock >= MIN_KNOCK_LENGTH):
            currently_in_a_knock = False

        if current_knock_count > 0:
            silence_duration = now - last_knock_time  # Reset the silence duration so we don't terminate early.
            # Finished a knock set
            if silence_duration > MIN_TIME_BETWEEN_SETS:
                current_sequence.append(current_knock_count)
                print(f"Added this digit: {current_knock_count}. Knock sequence is now: {current_sequence}")
                current_knock_count = 0
                    
        if current_sequence and time_since_last_knock > FINAL_SILENCE_REQUIRED:
            if secret_code_found_in_suffix(secret_code=secret_knock_code, entire_knock_sequence=current_sequence):
                print(f"Found the code ({secret_knock_code}) at the end of this current sequence of knocks ({current_sequence})")
                press_garage_door_button(pin_obj, press_time_sec=BUTTON_PRESS_TIME_SEC)
                current_sequence = [] # Reset for next attempt
            else:
                if (last_code_we_printed != current_sequence):
                    print(f"Rejected this code: {current_sequence}.    We expected this one: {secret_knock_code}\n")
                    last_code_we_printed = current_sequence # Just so we don't spam the logs with "rejected this code..." over and over
            
        if current_sequence and ((time_since_last_knock > FULL_RESET_TIMEOUT) or len(current_sequence) > max_knock_sets_to_record_before_resetting):
            current_sequence = []
        
        # ---------------------- PRINT THE PERCANTAGE NICELY ----------------
        if VERBOSE:
            audio_color = color_for_pct(pct)
            bar_chars: str = "#" * int(pct / (100 / MAX_VOLUME_BAR_WIDTH))
            blanks: str = " " * (MAX_VOLUME_BAR_WIDTH - len(bar_chars))

            right_aligned_rms = f"{rms:>8.6f}"  # The ":>8" should right-justify the result, I think it's (whole number + period + 6 fraction digits)

            # Here we will super inefficiently recalculate the max value every N runs through the loop, instead of on every loop
            if buffer_index % RECALCULATE_MAX_EVERY == 0:
                max_recent_volume = int(np.max(last_n_volumes))

            right_aligned_max = f"{max_recent_volume:>3}%"  # The ":>3" right-justifies the result
            max_color = color_for_pct(max_recent_volume)

            # Note the leading '\r\ to overwrite the same line in the terminal. Still makes a newline if the terminal is narrow.
            # ">5" right-aligns the percentage in a 5-digit-wide field.
            print(f"\rLevel: {audio_color}{pct:>5.1f}% {bar_chars}{blanks}{right_aligned_rms}"
                + f" {max_color}{right_aligned_max}{color.reset}", flush=True, end="")
            last_n_volumes[buffer_index] = int(pct)
            buffer_index = (buffer_index + 1) % BUFFER_SIZE  # Point to the next element...
            pass
        # -------------------- DONE with printing and printed-related bookkeeping -----

        pass # End of infinite loop
    return

def unit_tests():
    assert secret_code_found_in_suffix([1,2,3], [1,2,3]), "Exact should match"
    assert secret_code_found_in_suffix([1,2,3], [4, 1,2,3]), "Suffix should match"
    assert not secret_code_found_in_suffix([1,2,3], [1,2,3,4]), "Not a suffix: should not match"
    assert not secret_code_found_in_suffix([], [1,2,3,4]), "Empty secret code should not match"
    assert not secret_code_found_in_suffix([], []), "Empty/empty should not match"
    print("✅✅ Ad-hoc unit tests have passed")
    return

"""Press the button.

Args:
    gpio_pin_name: E.g. "GPIO4". Note that the GPIO pin name is usually not the same as the (physical) BOARD pin name.
                   The meaning of `on` and `off` are determined at the pin creation time: note that 'on' may actually represent
                   toggling the pin to 0 volts (see the `active_is_low` parameter in the initialization of the OutputDevice.)
    press_time_sec: How long to keep the states switched.
"""
def press_garage_door_button(pin_obj: gpiozero.OutputDevice | None, press_time_sec: int):
    if IS_DRY_RUN:
        print("⭐️⚠️⭐️ DRY RUN: NOT PUSHING THE BUTTON")
        return
    if pin_obj is None:
        print("⭐️⚠️⭐️ NOT PUSHING THE BUTTON: Could not initialize the GPIO pin. Check the logs.")
        return
    print("⭐️✅⭐️ ATTEMPTING TO OPEN THE GARAGE DOOR")
    pin_obj.on()  # Note that "on" MAY actually be the "low' voltage setting.
    time.sleep(press_time_sec)
    pin_obj.off()
    return

def init_audio_stream():
    return pyaudio.PyAudio().open(format=FORMAT, channels=CHANNELS, rate=RATE, input=True, frames_per_buffer=CHUNK)

def main():
    parser = argparse.ArgumentParser(description="Audio Handleroo")
    parser.add_argument("--code",                 type=str_with_digits_from_1_to_9_only, default="", help="Required. Secret code in single-digits of knocks. E.g. 123 = {1, 2, 3 knocks}.")
    parser.add_argument("--scale",                type=float,       default="1.0", help="[-Inf, Inf] Rescale the volume by this amount. Useful if the numbers seem 'off'")
    parser.add_argument("--knock",                type=int_percent, default="35", help="[1, 100] Threshold for detecting a 'knock' event, from 1 (percent) to 100 (percent).")
    parser.add_argument("--gpio_pin_to_pull_low", type=str,         default="GPIO4", help=f"The pin to pull LOW for {BUTTON_PRESS_TIME_SEC} second(s) to open/close the garage door.")

    parser.add_argument("--debug_test_button_now", action='store_true', help="Debug option. Press the button and then exit.")
    parser.add_argument("--run_unit_tests"       , action='store_true', help="Run our rather informal unit tests")
    parser.add_argument("--dry"                  , action='store_true', help="If true, do NOT actually activate the garage door opener")
    parser.add_argument("--verbose"              , action='store_true', help="If true, print the volume bar. This does A LOT of printing.")

    args = parser.parse_args()

    if args.run_unit_tests:
        unit_tests()

    global IS_DRY_RUN
    IS_DRY_RUN = args.dry
    global VERBOSE
    VERBOSE = args.verbose

    secret_code_as_list: list[int] = [int(x) for x in args.code]

    print("Listening for audio!")
    print(f"""  * Scaling factor is:   {args.scale}""")
    print(f"""  * Knock threshold is:  {args.knock}%""")
    print(f"""  * Knock code is:       {secret_code_as_list} ("--code={args.code}")""")
    print(f"""  * GPIO pin (pull low): {args.gpio_pin_to_pull_low}""")
    if IS_DRY_RUN:
        print(f"  * DRY RUN")
        pass

    try:
        # Note: If this is set up wrong, then this initialization may also open the garage door. Which we do not want.
        # The `initial_value` should always be False even if the pin is active on "high"
        pin_obj = gpiozero.OutputDevice(args.gpio_pin_to_pull_low, initial_value=False, active_high=(not PULL_PIN_LOW_TO_ACTIVATE))
    except gpiozero.exc.BadPinFactory:
        print("❌ Not able to actually get the real GPIO pin. This is OK if you're testing on a Mac, but it means we CANNOT actually control the garage door.")
        pin_obj = None

    if (args.debug_test_button_now):
        press_garage_door_button(pin_obj, BUTTON_PRESS_TIME_SEC)
        print("Exiting early because `--debug_test_button_now` was specified")
        return # Exit early

    try:
        handle_audio(rescale_volume=args.scale, knock_pct_threshold=args.knock, secret_knock_code=secret_code_as_list, pin_obj=pin_obj)
    except KeyboardInterrupt:
        print("\n\n⚠️ User pressed Ctrl-C, so we're quitting!")
    finally:
        pass

if __name__ == "__main__":
    main()

#!/usr/bin/env python3

# Example:
#  python garage_door_opener.py --knock=15 --code=123 --pause_time=0.70 --scale=3 --verbose --log=/tmp/garage_log.txt

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

import asyncio
import re
import sys
import time
import argparse
import datetime
from typing import Collection

import bleak  # For bluetooth device filtering
import gpiozero
import pyaudio
import numpy as np

IS_DRY_RUN: int = False
VERBOSE: bool   = False

LOG_FILE_PATH: str = ""

# Characters that we replace with '.' when checking for Bluetooth device names. Anything matching here will be replaced by an `UNSAFE_CHAR_REPLACEMENT` char.
UNSAFE_CHAR_MATCHER: str = "[^-~A-Za-z0-9_ .]"
SAFER_REPLACEMENT: str = "."

"""Escape sequence for terminal colors."""
class color:
    g     = '\033[92m' # GREEN
    y     = '\033[93m' # YELLOW
    r     = '\033[91m' # RED
    reset = '\033[0m'

BUTTON_PRESS_TIME_SEC:   float = 0.5  # How long to press the garage door button
PULL_PIN_LOW_TO_ACTIVATE: bool = True  # Hard-coded. Right now we hard-code that we set the pin to LOW to activate it. Could be refactored.

FORMAT                    = pyaudio.paFloat32
AUDIO_DTYPE               = np.float32
CHANNELS: int             = 1        # (Mono)
RATE: int                 = 44100    # 44.1kHz sampling rate
CHECKS_PER_SECOND: int    = 10  # How many times to check for the current audio volume every second. Can be pretty low.
CHUNK: int                = int(RATE / (1+CHECKS_PER_SECOND))   # Amount of data to read in one "chunk". Might be slightly off (not sure if this division is totally correct)
MAX_VOLUME_BAR_WIDTH: int = 50     # Max number of characters for the volume bar. If this is larger than the terminal width, then the bar won't overwrite itself (and the terminal will scroll)

# These are chosen somewhat empirically and probably should be command line args.
MIN_KNOCK_LENGTH: float        = 0.05  # Absolute minimum number of seconds to wait before counting another knock. If you set it too HIGH, then two knocks will be counted as one.

"""Prints iff VERBOSE (global variable) is True."""
def verbose_print(*args, **kwargs):
    if VERBOSE:
        print(*args, **kwargs)

"""Converts a percent (number) into a color. Loud = red, medium = yellow, quiet = green."""
def color_for_pct(p: int | float):
    if p <= 5:  # Arbitrary cutoffs here
        return color.g
    if p <= 20:
        return color.y
    return color.r

"""Turns "a ,  bb,ccc,d,  e" into ["a", "bb", "ccc", "d", "e"]. Note the stripped whitespace."""
def type_csv(in_arg) -> list[str]:
    if not in_arg:
        return [] # Empty list (instead of [''])
    items: list[str] = [x.strip() for x in str(in_arg).split(",")]
    return items

def replace_potentially_unsafe_chars(s: str, unsafe_char_regexp: str, replacement_char: str) -> str:
    return re.sub(unsafe_char_regexp, replacement_char, s)


def int_percent(x) -> int:
    try:
        x = int(x)
    except ValueError:
        raise argparse.ArgumentTypeError(f"{x} doesn't seem to be castable to an int!")
    if 1 <= x <= 100:
        return x
    raise argparse.ArgumentTypeError(f"{x} is not a valid value for this argument: it must be in the range [1, 100]")

def str_with_digits_from_1_to_9_only(x: str) -> str:
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


def handle_audio(rescale_volume: float, knock_pct_threshold: float, secret_knock_code: list[int], pin_obj: gpiozero.OutputDevice, required_pause_time: float):
    assert knock_pct_threshold >= 0 and knock_pct_threshold <= 100, "knock_pct_threshold must be in [0, 100]"
    assert isinstance(secret_knock_code, list) and all(isinstance(x, int) for x in secret_knock_code), "secret_knock_code should be an array of ints"

    currently_in_a_knock: bool  = False
    last_knock_time: float      = 0.0
    current_knock_count: int    = 0
    knock_seq: list[int] = []

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
        pct: float = 100.0 * min(1, rms)  # Convert to percent (clamped at 100%). This is a holdover from when I was doing everything in percent space.
        
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
                time_since_last_knock = now - last_knock_time  # Recompute time_since_last_knock for the benefit of the code below. (Probably there's a more elegant way to do this)

        # See if the current knock is now "done" being loud (volume < 50% of threshold AND it hasn't been a ludicrously short amount of time) (and prepare us to count the next one)
        if (pct < knock_pct_threshold * 0.5 and time_since_last_knock >= MIN_KNOCK_LENGTH):
            currently_in_a_knock = False

        if current_knock_count > 0:
            if time_since_last_knock > required_pause_time:  # Finished a set of knocks! In a perfect world, maybe we'd look at AVERAGE time between knocks and actually set this to (average * 2) or something.
                knock_seq.append(current_knock_count)
                print(f"{time_since_last_knock=} Added [{current_knock_count}], so now {knock_seq=}")
                current_knock_count = 0
        

        final_silence_required_before_activating_garage = required_pause_time  # Same as the normal pause time
        if knock_seq and time_since_last_knock > final_silence_required_before_activating_garage:
            if secret_code_found_in_suffix(secret_code=secret_knock_code, entire_knock_sequence=knock_seq):
                print(f"Found {secret_knock_code=}) at the end of {knock_seq=})")
                press_garage_door_button(pin_obj, press_time_sec=BUTTON_PRESS_TIME_SEC)
                knock_seq = [] # Reset for next attempt
            else:
                if (last_code_we_printed != knock_seq):
                    print(f"Rejected this code: {knock_seq}.    We expected this one: {secret_knock_code}\n")
                    last_code_we_printed = knock_seq # Just so we don't spam the logs with "rejected this code..." over and over
        
        full_reset_timeout: float = required_pause_time * 10   # Reset the code after this many seconds of silence. Possibly not actually needed.
        if knock_seq and ((time_since_last_knock > full_reset_timeout) or len(knock_seq) > max_knock_sets_to_record_before_resetting):
            knock_seq = []
        
        # ---------------------- PRINT THE VOLUME BAR NICELY ----------------
        if VERBOSE:
            audio_color = color_for_pct(pct)
            bar_chars: str = "#" * int(pct / (100 / MAX_VOLUME_BAR_WIDTH))
            blanks: str = " " * (MAX_VOLUME_BAR_WIDTH - len(bar_chars))
            right_aligned_rms = f"{rms:>8.6f}"  # The ":>8" should right-justify the result, I think it's (whole number + period + 6 fraction digits)
            # Note the leading '\r\ to overwrite the same line in the terminal. Still makes a newline if the terminal is narrow.
            # ">5" right-aligns the percentage in a 5-digit-wide field.
            print(f"\rLevel: {audio_color}{pct:>5.1f}% {bar_chars}{blanks}{right_aligned_rms} {color.reset}", flush=True, end="")
            pass
        # -------------------- DONE with printing and printed-related bookkeeping -----

        pass # End of infinite loop
    return

"""Press the button.

Args:
    gpio_pin_name: E.g. "GPIO4". Note that the GPIO pin name is usually not the same as the (physical) BOARD pin name.
                   The meaning of `on` and `off` are determined at the pin creation time: note that 'on' may actually represent
                   toggling the pin to 0 volts (see the `active_is_low` parameter in the initialization of the OutputDevice.)
    press_time_sec: How long to keep the states switched.
"""
def press_garage_door_button(pin_obj: gpiozero.OutputDevice | None, press_time_sec: float):
    try:
        if IS_DRY_RUN:
            print("⭐️⚠️⭐️ DRY RUN: NOT PUSHING THE BUTTON")
            open_and_write_to_log(LOG_FILE_PATH, "(Dry run, no action performed)")
        elif pin_obj is None:
            print("⭐️⚠️⭐️ NOT PUSHING THE BUTTON: Could not initialize the GPIO pin. Check the logs.")
            open_and_write_to_log(LOG_FILE_PATH, "(NOT PUSHING THE BUTTON: Could not initialize the GPIO pin.)")
        else:
            print("⭐️✅⭐️ ATTEMPTING TO OPEN THE GARAGE DOOR")
            pin_obj.on()  # Note that "on" MAY actually be the "low' voltage setting.
            time.sleep(press_time_sec)
            pin_obj.off()
            open_and_write_to_log(LOG_FILE_PATH, f"Pressed the garage door button.")
            pass
    except Exception as e:
        # Catch all exceptions and write another log line.
        open_and_write_to_log(LOG_FILE_PATH, f"Failure! Exeption was: {e}")
    return

def init_audio_stream():
    return pyaudio.PyAudio().open(format=FORMAT, channels=CHANNELS, rate=RATE, input=True, frames_per_buffer=CHUNK)

def open_and_write_to_log(log_file_path: str, message: str) -> None:
    if not log_file_path:
        verbose_print("Not writing to a log file, since none was specified")
        return

    try:
        with open(log_file_path, "a") as f:
            when: str = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            f.write(f"[{when}] {message}\n")  # E.g. "[1998-01-30] Snakes and cakes"
    except Exception as e:
        print(f"Couldn't append a log entry to {log_file_path=}: {e}")


async def print_all_nearby_bluetooth_devices() -> None:
    scan_time_sec: float = 5
    print(f"Spending {scan_time_sec} seconds to scan for nearby BLUETOOTH devices...")
    try:
        devices = await bleak.BleakScanner.discover(return_adv=True, timeout=scan_time_sec)
    except bleak.exc.BleakError as e:
        print(f"BLUETOOTH is unavailable for some reason. Maybe it was turned off? The exception: {e}")
        raise e
    
    if not devices:
        print("No Bluetooth devices were found in range! This is unusual.")
        return
    
    for address, (device, adv_data) in devices.items():
        print(f"ID: {address}  |   RSSI_SIGNAL: {adv_data.rssi:>4}  |   NAME: {device.name}")
        #                                                     ^ (right-aligned text)

    return


"""Gets the 'friendly' user-facing (and user-configurable) IDs for nearby Bluetooth devices."""
async def get_nearby_bluetooth_device_names(sanitize_output: bool=True) -> set[str]:
    scan_time_sec: float = 1
    names_found: set[str] = set()
    print(f"Spending {scan_time_sec} seconds to scan for nearby BLUETOOTH devices...")
    try:
        devices = await bleak.BleakScanner.discover(return_adv=True, timeout=scan_time_sec)
    except bleak.exc.BleakError as e:
        print(f"BLUETOOTH is unavailable for some reason. Maybe it was turned off? The exception: {e}")
        raise e
    
    if not devices:
        print("No Bluetooth devices were found in range! This is unusual, but isn't an error")
    else:
        for address, (device, adv_data) in devices.items():
            # FYI: 'address' is in a format like 'AA7AD434-B8DD-19AA-44DC-623E444A718A'
            print(f"ID: {address}  |   RSSI_SIGNAL: {adv_data.rssi:>4}  |   NAME: {device.name}")
            #                                                     ^ (right-aligned text)
            if device.name():
                names_found.add(device.name())  # Add non-blank and non-"None" names only. Theoretically might add a name that was just several spaces.

    if sanitize_output:
       return {replace_potentially_unsafe_chars(x, unsafe_char_regexp=UNSAFE_CHAR_MATCHER, replacement_char=SAFER_REPLACEMENT) for x in names_found}
    else:
       return names_found

"""Returns whether or not our 'nearby Bluetooth device name' requirement is met.

required_names: If empty, then this function ALWAYS returns true (always matches). Otherwise, returns whether any regex in 'required_name_regexes' matches any of the actual names.
"""
def bluetooth_match_satisfied(required_name_regexes: Collection[str], actual_names: Collection[str]) -> bool:
    if not required_name_regexes:
        return True # No requirement, so always return true
    for pattern in required_name_regexes:
        for name in actual_names:
            if re.fullmatch(pattern, name):
                verbose_print(f"""The nearby Bluetooth device named '{name}' matched the required regex pattern '{pattern}'""")
                return True
    verbose_print(f"""No nearby Bluetooth device matched any of our naming requirements!""")
    return False


def unit_tests() -> None:
    assert secret_code_found_in_suffix([1,2,3], [1,2,3]), "Exact should match"
    assert secret_code_found_in_suffix([1,2,3], [4, 1,2,3]), "Suffix should match"
    assert not secret_code_found_in_suffix([1,2,3], [1,2,3,4]), "Not a suffix: should not match"
    assert not secret_code_found_in_suffix([], [1,2,3,4]), "Empty secret code should not match"
    assert not secret_code_found_in_suffix([], []), "Empty/empty should not match"
    assert type_csv("") == [], "Empty string should make an empty list, not ['']"
    assert type_csv("  11,  22  ,33  ") == ["11","22","33"], "Normal splitting"
    assert type_csv("  aaa a ") == ["aaa a"], "Normal splitting"
    assert bluetooth_match_satisfied(required_name_regexes=[".*"],     actual_names=[""]), "Dot-star should match the empty string."
    assert bluetooth_match_satisfied(required_name_regexes=[],         actual_names=[]), "Empty required name regexes should even match an empty list of actual names."
    assert not bluetooth_match_satisfied(required_name_regexes=[".*"], actual_names=[]), "Empty actual names doesn't match a non-empty required regex."
    assert bluetooth_match_satisfied(required_name_regexes=["zzz", ".*Cool.*"], actual_names=["ACoolPhone"])
    assert bluetooth_match_satisfied(required_name_regexes=["zzz", "Cool.*"], actual_names=["CoolPhone"])
    assert not bluetooth_match_satisfied(required_name_regexes=["zzz", ".*Cool"], actual_names=["CoolPhone"])
    assert not bluetooth_match_satisfied(required_name_regexes=["a"], actual_names=["aa"]), "'a' is an exact match only"
    assert not bluetooth_match_satisfied(required_name_regexes=["A"], actual_names=["a"]), "'Case matters"
    assert not bluetooth_match_satisfied(required_name_regexes=["zzz", ".*q"], actual_names=["CoolPhone", "zzzz", "zz", "z"]), "No match for three 'z's"
    assert bluetooth_match_satisfied(required_name_regexes=["Z", ".*Cool Phone"], actual_names=["Joe's Cool Phone"]), "Should match"
    assert bluetooth_match_satisfied(required_name_regexes=["Z", ".*Cool[ ]?Phone"], actual_names=["MyCoolPhone"]), "Should match due to the optional space"
    assert replace_potentially_unsafe_chars("Snåkeßß🟨a b c", unsafe_char_regexp="[^a-z ]", replacement_char="#") == "#n#ke###a b c"
    assert replace_potentially_unsafe_chars("", unsafe_char_regexp="[^a-z]", replacement_char="#") == ""
    print("✅✅ Ad-hoc unit tests have passed")
    return


def main():
    parser = argparse.ArgumentParser(description="Audio Handleroo")
    parser.add_argument("--code",                 type=str_with_digits_from_1_to_9_only, default="", help="Required. Secret code in single-digits of knocks. E.g. 123 = {1, 2, 3 knocks}.")
    parser.add_argument("--scale",                type=float,       default="1.0",   help="[-Inf, Inf] Rescale the volume by this amount. Useful if the numbers seem 'off'")
    parser.add_argument("--knock",                type=int_percent, default="35",    help="[1, 100] Threshold for detecting a 'knock' event, from 1 (percent) to 100 (percent).")
    parser.add_argument("--pause_time",           type=float,       default="0.70",  help="Required pause time between 'sets' of knocks, in seconds. 0.5 is typically too short, and 0.8 is typically too long. Default is 0.70.")
    parser.add_argument("--gpio_pin_to_pull_low", type=str,         default="GPIO4", help=f"The pin to pull LOW for {BUTTON_PRESS_TIME_SEC} second(s) to open/close the garage door.")
    parser.add_argument("--log",                  type=str,         default="", help="Default: no log. If specified, log garage-door button-press actions to this file. Note that opening/closing cannot be distinguished.")

    parser.add_argument("--scan_for_bluetooth",             action='store_true', help="If true, runs a scanner ONCE for 5 seconds and then prints the results and exits. This is for helping you figure out the IDs associated with your specific Bluetooth devices of interest.")
    parser.add_argument("--required_bluetooth_names",       type=type_csv,     help="One or more FULL MATCH regexes (e.g. '.*Cool.*' matches 'ACoolPhone', but 'Cool' only matches 'Cool' verbatim). If non-empty, require that these device NAMES (not IDs!) be nearby. WARNING: this uses the easily-detected-and-spoofed COMMON names, not the Bluetooth ID. (This is because Apple devices apparently rotate their Bluetooth ID, and I want something that is static.)")
    parser.add_argument("--allow_nonascii_bluetooth_names", action='store_true', help=f"If true, allow bluetooth names to contain 'surprising' characters. By default, we replace anything that matches '{UNSAFE_CHAR_MATCHER}' with '{SAFER_REPLACEMENT}'.")

    parser.add_argument("--debug_test_button_now", action='store_true', help="Debug option. Press the button and then exit.")
    parser.add_argument("--run_unit_tests"       , action='store_true', help="Run our rather informal unit tests")
    parser.add_argument("--dry"                  , action='store_true', help="If true, do NOT actually activate the garage door opener")
    parser.add_argument("--verbose"              , action='store_true', help="If true, print the volume bar. This does A LOT of printing.")

    args = parser.parse_args()

    if args.run_unit_tests:
        unit_tests()

    if args.scan_for_bluetooth:
        asyncio.run(print_all_nearby_bluetooth_devices())
        print("DONE listing nearby bluetooth devices. Exiting WITHOUT running anything else.")
        sys.exit(0)

    assert args.pause_time >= 0.2 and args.pause_time <= 10, f"Your --pause_time must be some 'normal' number of seconds between 0.2 and 10 seconds, not `{args.pause_time=}`"

    global IS_DRY_RUN
    IS_DRY_RUN = args.dry
    global VERBOSE
    VERBOSE = args.verbose
    global LOG_FILE_PATH
    LOG_FILE_PATH = args.log

    secret_code_as_list: list[int] = [int(x) for x in args.code]

    print("Listening for audio!")
    print(f"""  * Scaling factor is:      {args.scale}""")
    print(f"""  * Knock threshold is:     {args.knock}% (of max volume)""")
    print(f"""  * Min pause between sets: {args.pause_time} seconds""")
    print(f"""  * Knock code is:          {secret_code_as_list} ("--code={args.code}")""")
    print(f"""  * GPIO pin (pull low):    {args.gpio_pin_to_pull_low}""")
    print(f"""  * Log file path (if any): {args.log}""")
    print(f"""  * Required Bluetooth names: {args.required_bluetooth_names}""")
    if IS_DRY_RUN:
        print(f"  * DRY RUN")
        pass

    # Verify that we can write to the log, so we don't fail the first time we ACTUALLY want to write to the log.
    open_and_write_to_log(LOG_FILE_PATH, "<Starting program>")


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
        handle_audio(rescale_volume=args.scale, knock_pct_threshold=args.knock, secret_knock_code=secret_code_as_list, pin_obj=pin_obj, required_pause_time=args.pause_time)
    except KeyboardInterrupt:
        print("\n\n⚠️ User pressed Ctrl-C, so we're quitting!")
    finally:
        open_and_write_to_log(LOG_FILE_PATH, "<Exiting program after running the normal detection loop>")
        pass

if __name__ == "__main__":
    main()

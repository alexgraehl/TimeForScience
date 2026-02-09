#!/usr/bin/env python3

# Example:
#  python garage_door_opener.py --knock=15 --code=123 --pause_time=0.70 --scale=3 --verbose --verbose_bluetooth --required_bluetooth_names=".*CoolPhone.*,.*MyPhone.*" --log=/tmp/garage_log.txt --run_unit_tests

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
import enum
import re
import sys
import time
import argparse
import datetime
from dataclasses import dataclass
from typing import Collection, NoReturn, Sequence, Generator

import bleak  # For bluetooth device filtering
import gpiozero
import pyaudio
import numpy as np

class BluetoothIdType(enum.StrEnum):
    NAME = "BT_Name"  # Indicates that this ID is a "common" name for a Bluetooth device, like "Joe's iPhone 17"
    ADDR = "BT_Addr"  # Indicates that this ID is a Bluetooth address, with a format like "AA9DBB91-1442-ABCD-EF12-AAAABBBE8CCC"

@dataclass
class BluetoothDeviceInfo:
    id_type: BluetoothIdType
    time_last_seen: float

DeviceDictType = dict[str, BluetoothDeviceInfo]

# This will be modified as a GLOBAL variable by the bluetooth scanner task.
RECENT_BLUETOOTH_DEVICES      : DeviceDictType = {}
MAX_BT_DEVICES_TO_REMEMBER    : int = 250  # If `RECENT_BLUETOOTH_DEVICES` is larger than this, then we'll prune old entries. Just to keep the device list from accumulating forever.
DEVICE_PRUNING_THRESHOLD_SEC  : int = 120  # If we exceed the `MAX_BT_DEVICES_TO_REMEMBER`, then prune any entries older than this amount.
BLUETOOTH_REQUIRED_RECENCY_SEC: int = 60   # Must have seen the device in the last N seconds to count it as "recently" enough to open the garage door.

# Time to give up control in the audio 'main' loop. Must be quite LOW in order to reliably pick up knocks.
# Note that this CAN be zero if you want, which means "do it as fast as possible."
SLEEP_TIME_IN_AUDIO_LOOP           : float = 0.001  
SLEEP_TIME_IN_BLUETOOTH_FINDER_LOOP: float = 1  # Check nearby bluetooth devices every this many seconds. A good value is probably like 10-20.

IS_DRY_RUN       : int  = False
VERBOSE          : bool = False
VERBOSE_BLUETOOTH: bool = False  # Super verbose logging specifically for Bluetooth
LOG_FILE_PATH    : str  = ""     # (If this is the empty string, then don't log)

# Characters that we replace with '.' when checking for Bluetooth device names. Anything matching here will be replaced by an `UNSAFE_CHAR_REPLACEMENT` char.
UNSAFE_CHAR_MATCHER: str = "[^-~A-Za-z0-9_ .]"
SAFER_REPLACEMENT  : str = "."

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

"""Prints iff VERBOSE_BLUETOOTH (global variable) is True."""
def verbose_bluetooth_print(*args, **kwargs):
    if VERBOSE_BLUETOOTH:
        print(*args, **kwargs)


"""Prunes the RECENT_BLUETOOTH_DEVICES dictionary, in case it gets too large."""
def pruned_copy_only_entries_after_time(device_dict: DeviceDictType, cutoff_time: float) -> DeviceDictType:
    return {key:val for key, val in device_dict.items() if val.time_last_seen >= cutoff_time}

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

def type_tuple_of_digits_1_to_9_from_str(in_arg: str) -> tuple[int, ...]:  # Get tuple of digits from 1 through 9 (excluding zero!). E.g. input of "123" becomes [1, 2, 3]. CANNOT contain a 0.
    if not in_arg or (not in_arg.isdigit()) or ('0' in in_arg):
        raise argparse.ArgumentTypeError(f"This argument requires a string input consisting of the digits 1-9 only (not zero). The INCORRECT argument we got was \"{in_arg}\".")
    return tuple(int(d) for d in in_arg)  # Return a tuple of ints


"""Check to see if suffix is both non-empty AND is a suffix of `full_seq`.

E.g. if suffix = [1,2] and full_seq = [9,9,1,4,2,1,2], then it's true (ends in '1,2')
Notably always returns False if the input suffix is zero-length, which may be surprising. (i.e. `[]` is not a valid suffix of `[1,2]`).
"""
def suffix_found(suffix: Sequence[int], full_seq: Sequence[int]):
    if not suffix or len(suffix) > len(full_seq):
        return False  # Not long enough!
    return tuple(full_seq[-len(suffix):]) == tuple(suffix)  # Check just the suffix. Tuple-ify both args to ensure equality works even if one input is a list and the other is a tuple.


async def infinite_loop_bluetooth_scanner(required_bluetooth_regexps: Collection[str], sanitize_bt_device_names: bool):
    """When the scanners sees a particular Bluetooth device, it indicates when it was last seen.
    Every so often, we'll clear out anything in the RECENT_DEVICES map that hasn't been seen in a while.
    """
    def bluetooth_detected_device_callback(device, adv_data):
        global RECENT_BLUETOOTH_DEVICES  # <-- Required since we may be overwriting RECENT_BLUETOOTH_DEVICES
        verbose_bluetooth_print(f"Saw this Bluetooth info:      {device.address=}")
        verbose_bluetooth_print(f"                 ...and:         {device.name=}")
        verbose_bluetooth_print(f"                 ...and: {adv_data.local_name=}")
        common_name: str = device.name or adv_data.local_name or ""

        if common_name:
            if sanitize_bt_device_names:
                common_name = replace_potentially_unsafe_chars(common_name, unsafe_char_regexp=UNSAFE_CHAR_MATCHER, replacement_char=SAFER_REPLACEMENT)
                pass
            now: float = time.time()
            RECENT_BLUETOOTH_DEVICES[common_name] = BluetoothDeviceInfo(id_type=BluetoothIdType.NAME, time_last_seen=now)
            is_recent_devices_too_big: bool = (len(RECENT_BLUETOOTH_DEVICES) > MAX_BT_DEVICES_TO_REMEMBER)
            if is_recent_devices_too_big:  # We 'delete' the old items by just making a new dictionary with only the passing-the-filter items.
                RECENT_BLUETOOTH_DEVICES = pruned_copy_only_entries_after_time(RECENT_BLUETOOTH_DEVICES, cutoff_time=(now - DEVICE_PRUNING_THRESHOLD_SEC))

            if VERBOSE_BLUETOOTH:
                print_device_dict(RECENT_BLUETOOTH_DEVICES)
                # Check to see if any of the devices is in our special whitelist
                recency_cutoff: float = now - BLUETOOTH_REQUIRED_RECENCY_SEC
                matched_device_name: str | None = first_matched_bt_device(d=RECENT_BLUETOOTH_DEVICES, allowed_regexps=required_bluetooth_regexps, time_cutoff=recency_cutoff)
                if not matched_device_name:
                    verbose_bluetooth_print("NO MATCHING DEVICES")
                else:
                    device_info = RECENT_BLUETOOTH_DEVICES[matched_device_name]
                    time_since_seen: float = device_info.time_last_seen - now
                    verbose_bluetooth_print(f"""Found a match: device "{matched_device_name}" was seen at {device_info.time_last_seen:.2f} ({time_since_seen:.2f} seconds ago)""")

        return
        
    continuous_bluetooth_scanner = bleak.BleakScanner(detection_callback=bluetooth_detected_device_callback)
    await continuous_bluetooth_scanner.start()
    try:
        while True:
            await asyncio.sleep(delay=SLEEP_TIME_IN_BLUETOOTH_FINDER_LOOP)
    finally:
        verbose_bluetooth_print("\nShutting down the Bluetooth scanner...")
        await continuous_bluetooth_scanner.stop()  # Prevents segfault on exit

"""Note that 'infinite_loop_audio_listener' is basically the main loop. See the "while True" in it."""
async def infinite_loop_audio_listener(rescale_volume: float, knock_pct_threshold: float, secret_knock_code: Sequence[int], pin_obj: gpiozero.OutputDevice, required_pause_time: float, required_bluetooth_regexps: Collection[str]):
    assert knock_pct_threshold >= 0 and knock_pct_threshold <= 100, "knock_pct_threshold must be in [0, 100]"

    currently_in_a_knock: bool  = False
    last_knock_time: float      = 0.0
    current_knock_count: int    = 0
    knock_seq: list[int] = []

    max_knock_sets_to_record_before_resetting: int = len(secret_knock_code) * 20  # Just make sure we don't keep track of knock sets forever if there's a woodpecker or something
    last_code_we_printed = []  # Just remember the last code we updated the user about, so we don't keep spamming the logs

    audio_stream = init_audio_stream()
    try:
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
                    print(f"{time_since_last_knock=} Added [{current_knock_count}], so now {knock_seq=} (we will look for {secret_knock_code})")
                    current_knock_count = 0
            
            final_silence_required_before_activating_garage = required_pause_time  # Same as the normal pause time
            if knock_seq and time_since_last_knock > final_silence_required_before_activating_garage:
                if suffix_found(suffix=secret_knock_code, full_seq=knock_seq):
                    # Check to see if a nearby bluetooth device is found, if applicable.
                    print(f"Found {secret_knock_code=}) at the end of {knock_seq=})")
                    
                    msg: str = ""
                    allowed_device: str | None = None
                    if not required_bluetooth_regexps:
                        msg = "OK to open the door (we are not checking any specific nearby Bluetooth devices)"
                    else:
                        allowed_device = first_matched_bt_device(d=RECENT_BLUETOOTH_DEVICES, allowed_regexps=required_bluetooth_regexps, time_cutoff=now - BLUETOOTH_REQUIRED_RECENCY_SEC)                        
                        if not allowed_device:
                            msg = "NOT opening the garage door (despite the correct code), because we failed to find a mandatory nearby Bluetooth device."
                        else:
                            msg = f"OK to open the door, because we found this eligible nearby Bluetooth device: {allowed_device}"
                            pass
                    
                    open_and_write_to_log(LOG_FILE_PATH, msg)

                    if (not required_bluetooth_regexps) or allowed_device:
                        # We're good to open the garage door
                        press_garage_door_button(pin_obj, press_time_sec=BUTTON_PRESS_TIME_SEC)
                    
                    # Regardless of the Bluetooth device presence above, reset the code for the next attempt.
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

            await asyncio.sleep(0.05) # End of infinite loop
            pass # End of "while True"
    finally:
        audio_stream.stop_stream()
        audio_stream.close()  # Prevents segfault on exit

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
            msg: str = f"[{when}] {message}"
            f.write(msg + "\n")  # E.g. "[1998-01-30] Snakes and cakes"
            verbose_print(msg)
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

"""Returns the name of first nearby device that matches the name requirements AND the time cutoff requirement (i.e. "recently_seen-ness").

Returns None if no nearby device has a matching name and was also seen 'recently' enough.

Args:
    d: The dictionary of Bluetooth devices that we've seen semi-recently.
    allowed_regexps: The devices names that are 'ok'.
    time_cutoff: Only allow a device to be returned if it was seen at this time or later. Don't return 'stale' devices that were seen hours ago.
"""
def first_matched_bt_device(d: DeviceDictType, allowed_regexps: Collection[str], time_cutoff: float) -> str | None:
    for device_name, device_data in d.items():
        if device_data.time_last_seen < time_cutoff:
            verbose_bluetooth_print(f"""Not checking the Bluetooth device '{device_name}', because it wasn't seen recently enough.""")
            continue # Next device please!
        for pattern in allowed_regexps:
            if re.fullmatch(pattern, device_name):
                verbose_bluetooth_print(f"""🟦✅🟦 The nearby Bluetooth device named '{device_name}' matched the required regex pattern '{pattern}' (and was seen on/after '{time_cutoff:.1f}')""")
                return device_name
    verbose_bluetooth_print(f"""🟦❌🟦 No nearby Bluetooth device matched any of our naming requirements and also {time_cutoff=}!""")
    return None # No match

"""Helper method to print a DeviceDictTable in tabular form to STDOUT."""
def print_device_dict(d: DeviceDictType):
    now: float = time.time()
    for name,v in d.items():
        print(f"""{v.time_last_seen::>13.1f}  | {(now - v.time_last_seen):>8.1f} sec ago  |  {v.id_type}  |  {name}""")


def unit_tests() -> None:
    # Note that we're dubiously using the current time (time.time()) in the unit tests, thus making them non-hermetic. Woe!
    now = time.time()
    fake_devices: DeviceDictType = {
        "Joe's iPhone 17": BluetoothDeviceInfo(
            id_type=BluetoothIdType.NAME,
            time_last_seen=(now - 10.1) # Las seen 10.1 seconds ago (recent)
        ),
        "AA:BB:CC:DD": BluetoothDeviceInfo(
            id_type=BluetoothIdType.ADDR,
            time_last_seen=(now - 3700) # Last seen over 1 hour ago (old)
        ),
        "Ghost Pepper": BluetoothDeviceInfo(
            id_type=BluetoothIdType.NAME,
            time_last_seen=(now - 0.5)  # Last seen 0.5 seconds ago (recnt)
        )
    }

    assert suffix_found([1,2,3], [1,2,3]), "Exact should match"
    assert suffix_found(list([2,3]), tuple([1,2,3])), "List and tuple should be interchangeable"
    assert suffix_found(tuple([2,3]), list([1,2,3])), "List and tuple should be interchangeable"
    assert suffix_found([1,2,3], [4, 1,2,3]), "Suffix should match"
    assert not suffix_found([1,2,3], [1,2,3,4]), "Not a suffix: should not match"
    assert not suffix_found([], [1,2,3,4]), "Empty secret code should not match"
    assert not suffix_found([], []), "Empty/empty should not match"
    assert type_csv("") == [], "Empty string should make an empty list, not ['']"
    assert type_csv("  11,  22  ,33  ") == ["11","22","33"], "Normal splitting"
    assert type_csv("  aaa a ") == ["aaa a"], "Normal splitting"
    assert type_tuple_of_digits_1_to_9_from_str("123")  == tuple((1,2,3))
    assert type_tuple_of_digits_1_to_9_from_str("4321") == tuple((4,3,2,1))
    assert type_tuple_of_digits_1_to_9_from_str("3")    == tuple((3,))
    assert bool(first_matched_bt_device(fake_devices, allowed_regexps=[".*"], time_cutoff=now-9999)), "Should find a device (we don't care which one)"
    assert not (first_matched_bt_device(fake_devices, allowed_regexps=[".*"], time_cutoff=now+4444)), "Should NOT find any device, due to the too-strict time cutoff (in the future)"
    assert "Joe's iPhone 17" == first_matched_bt_device(fake_devices, allowed_regexps=["Joe.s iPhone.*"], time_cutoff=now-100), "Should find a device"
    assert "AA:BB:CC:DD" == first_matched_bt_device(fake_devices, allowed_regexps=["ZZZ","AA.*"], time_cutoff=now - 9999), "Finds the `AA:BB:CC:DD` device (even though it's the ADDRESS and not the NAME). In practice this won't actually be found since we are not adding the addresses to the device dict currently."
    assert "Ghost Pepper" == first_matched_bt_device(fake_devices, allowed_regexps=[".*joe.*", ".*Pep.*"], time_cutoff=now - 100), "Finds Ghost Pepper"
    assert "Joe's iPhone 17" == first_matched_bt_device(fake_devices, allowed_regexps=["qq", "zz", "Joe.s iPhone.*"], time_cutoff=now - 100), "Simple match of a recently-seen device"
    assert not first_matched_bt_device(fake_devices, allowed_regexps=["Joe.s iPhone.*"], time_cutoff=now), "Disqualified by time cutoff"
    assert not first_matched_bt_device(fake_devices, allowed_regexps=[".*joe.*", ".*pep.*"], time_cutoff=now - 100), "Case matters, so no match"
    assert not first_matched_bt_device({},           allowed_regexps=[".*"], time_cutoff=now - 9999), "No devices, so no match"
    assert not first_matched_bt_device(fake_devices, allowed_regexps=[], time_cutoff=now - 9999), "No allowed regexps, so no match"

    assert replace_potentially_unsafe_chars("Snåkeßß🟨a b c", unsafe_char_regexp="[^a-z ]", replacement_char="#") == "#n#ke###a b c"
    assert replace_potentially_unsafe_chars("", unsafe_char_regexp="[^a-z]", replacement_char="#") == ""
    print("✅✅ Ad-hoc unit tests have passed")
    return


async def thread_runner(rescale_volume: float, knock_pct_threshold: int, secret_knock_code: Sequence[int], pin_obj, required_pause_time: float, required_bluetooth_regexps: Collection[str], sanitize_bt_device_names: bool):
    bluetooth_scanner_task: asyncio.Task[NoReturn] | None = None
    if required_bluetooth_regexps:
        # Only run the scanner if we are actually going to check the nearby bluetooth devices at all
        bluetooth_scanner_task = asyncio.create_task(infinite_loop_bluetooth_scanner(required_bluetooth_regexps=required_bluetooth_regexps, sanitize_bt_device_names=sanitize_bt_device_names))  # Effectively runs in parallel
    else:
        print("Note: not running the bluetooth scanner, since no `required_bluetooth_names` were specified.")
    
    try:
        # Audio handler is treated as the 'main' loop. It normally never returns.
        await infinite_loop_audio_listener(  
            rescale_volume=rescale_volume, 
            knock_pct_threshold=knock_pct_threshold, 
            secret_knock_code=secret_knock_code, 
            pin_obj=pin_obj, 
            required_pause_time=required_pause_time,
            required_bluetooth_regexps=required_bluetooth_regexps)
    finally:
        if bluetooth_scanner_task:
            bluetooth_scanner_task.cancel()  # Stop the bluetooth scanner task if/when the audio "main" loop exits


def main():
    parser = argparse.ArgumentParser(description="Audio Handleroo")
    parser.add_argument("--code",                 type=type_tuple_of_digits_1_to_9_from_str, default=(), help="Required. Secret code in single-digits of knocks. E.g. 123 = {1, 2, 3 knocks}.")
    parser.add_argument("--scale",                type=float,       default="1.0",   help="[-Inf, Inf] Rescale the volume by this amount. Useful if the numbers seem 'off'")
    parser.add_argument("--knock",                type=int_percent, default="35",    help="[1, 100] Threshold for detecting a 'knock' event, from 1 (percent) to 100 (percent).")
    parser.add_argument("--pause_time",           type=float,       default="0.70",  help="Required pause time between 'sets' of knocks, in seconds. 0.5 is typically too short, and 0.8 is typically too long. Default is 0.70.")
    parser.add_argument("--gpio_pin_to_pull_low", type=str,         default="GPIO4", help=f"The pin to pull LOW for {BUTTON_PRESS_TIME_SEC} second(s) to open/close the garage door.")
    parser.add_argument("--log",                  type=str,         default="", help="Default: no log. If specified, log garage-door button-press actions to this file. Note that opening/closing cannot be distinguished.")

    parser.add_argument("--scan_for_bluetooth",             action='store_true', help="If true, runs a scanner ONCE for 5 seconds and then prints the results and exits. This is for helping you figure out the IDs associated with your specific Bluetooth devices of interest.")
    parser.add_argument("--required_bluetooth_names",       type=type_csv      , default=[], help="If empty (default), then we don't check Bluetooth device names. One or more case-sensitive comma-separated FULL MATCH regexps (e.g. '.*Cool.*' matches 'ACoolPhone', but 'Cool' only matches 'Cool' verbatim). If non-empty, require that these device NAMES (not IDs!) be nearby. WARNING: this uses the easily-detected-and-spoofed COMMON names, not the Bluetooth ID. (This is because Apple devices apparently rotate their Bluetooth ID, and I want something that is static.). A good example that would match 'Jane's iPhone 12' and 'Joe's iPhone SE 8' would be: 'Jane.s.iPhone.*,Joe.s.iPhone.*'. If you need to see what devices are around, run this script with --scan_for_bluetooth (which will print nearby device details and then exit).")
    parser.add_argument("--allow_nonascii_bluetooth_names", action='store_true', help=f"If true, allow bluetooth names to contain 'surprising' characters. By default, we replace anything that matches '{UNSAFE_CHAR_MATCHER}' with '{SAFER_REPLACEMENT}'.")

    parser.add_argument("--debug_test_button_now", action='store_true', help="Debug option. Press the button and then exit.")
    parser.add_argument("--run_unit_tests"       , action='store_true', help="Run our rather informal unit tests")
    parser.add_argument("--dry"                  , action='store_true', help="If true, do NOT actually activate the garage door opener")
    parser.add_argument("--verbose"              , action='store_true', help="If true, print the volume bar. This does A LOT of printing.")
    parser.add_argument("--verbose_bluetooth"    , action='store_true', help="If true, add ultra-verbose printing for Bluetooth devices detected locally. Does a ton of printing.")

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
    global VERBOSE_BLUETOOTH
    VERBOSE_BLUETOOTH = args.verbose_bluetooth
    global LOG_FILE_PATH
    LOG_FILE_PATH = args.log

    print("Listening for audio!")
    print(f"""  * Scaling factor is:      {args.scale}""")
    print(f"""  * Knock threshold is:     {args.knock}% (of max volume)""")
    print(f"""  * Min pause between sets: {args.pause_time} seconds""")
    print(f"""  * Knock code is:          {args.code}""")
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
        asyncio.run(thread_runner(rescale_volume=args.scale, knock_pct_threshold=args.knock, secret_knock_code=args.code, pin_obj=pin_obj, required_pause_time=args.pause_time, required_bluetooth_regexps=args.required_bluetooth_names, sanitize_bt_device_names=(not args.allow_nonascii_bluetooth_names)))
    except KeyboardInterrupt:
        print("\n\n⚠️ User pressed Ctrl-C, so we're quitting!")
    finally:
        open_and_write_to_log(LOG_FILE_PATH, "<Exiting program after running the normal detection loop>")
        if args.log:
            print(f"📝 FYI: You may want to check the logs at\n       tail -n 50 {args.log}")
        pass

if __name__ == "__main__":
    main()

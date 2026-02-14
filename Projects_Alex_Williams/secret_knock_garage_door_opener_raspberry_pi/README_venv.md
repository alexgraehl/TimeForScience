# Readme about the `garage_door_opener.py`

## About `venv/`

Note that running this requires a virtual environment (or all the libraries in your system environment).

Check the header of the python program for details on how to set this up.

Then, normally every time you have to run this:

```shell
   source venv/bin/activate   # (To pick up `pyaudio` and `numpy`)
```

## Bluetooth authentication: How to set up an iPhone specifically

iPhones don't broadcast their "friendly" (user facing) names by default, apparently. This is a security thing.

If you are trying to use an iPhone with the Bluetooth authentication process, you have to first PAIR and TRUST the iPhone with your device via Bluetooth at least ONE time. From that time onward, the phone name ("Joe's iPhone") should show up.

### Pairing and trusting an iPhone from a Raspberry Pi

1. On the Pi, run `bluetoothctl` in a terminal.
2. From inside the `bluetoothctl` command line utility, type each of these commands:
3. `power on`
4. `agent on`
5. `default agent`
6. `discoverable on`
7. `pairable on`
8. Then go to your iPhone and find the Raspberry Pi in `Settings -> Bluetooth`. It might take a minute to find it.
9. On the iPhone, connect to the Raspberry Pi. On the iPhone, it will say something like "Connect to this device (Code: 98989)?" Tap 'yes'.
10. On the Raspberry Pi, the terminal should say something like 'Trust this device? (Code: 98989)'. Type 'yes', assuming the code matches.
11. You might have to type 'yes' a few times for some reason.
12. Find the MAC address of your iPhone. It will be something like AA:BB:11:22:33 or whatever. It should have scrolled by in the terminal above.
13. Type `trust $THAT_MAC_ADDRESS`, e.g. `trust AA:BB:11:22:33`.

Now when you run the bluetooth scanner (--scan), it should print both the MAC address AND the iPhone name. The iPhone "friendly" name can now also be used in the `--required_bluetooth_names` argument.

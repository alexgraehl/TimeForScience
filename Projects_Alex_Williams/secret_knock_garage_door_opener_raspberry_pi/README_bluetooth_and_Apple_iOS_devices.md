# Readme about using the Bluetooth identification system with Apple iOS devices (e.g. the iPhone)


## Bluetooth authentication: How to set up an iPhone specifically

It turns out that iPhones don't broadcast their "friendly" (user facing) names by default. Additionally, you can't use the MAC address reliably either, because it "rotates" (changes) after a few weeks. This is a security thing.

If you are trying to use an iPhone with the Bluetooth authentication process, you have to first PAIR and TRUST the iPhone with your device via Bluetooth at least ONE time. From that time onward, the phone name ("Joe's iPhone") should show up.


### Pairing and trusting an iPhone from a Raspberry Pi

1. On the Raspberry Pi, run `bluetoothctl` in a terminal.
2. From inside the `bluetoothctl` command line utility, type each of these commands:
   * `power on`
   * `agent on`
   * `default agent`
   * `discoverable on`
   * `pairable on`
3. Then go to your iPhone and find the Raspberry Pi in `Settings -> Bluetooth`. It might take a minute to find it.
4. On the iPhone, connect to the Raspberry Pi. On the iPhone, it will say something like "Connect to this device (Code: 98989)?" Tap 'yes'.
5. On the Raspberry Pi, the terminal should say something like 'Trust this device? (Code: 98989)'. Type 'yes', assuming the code matches.
6. You might have to type 'yes' a few times for some reason.
7. Find the MAC address of your iPhone. It will be something like AA:BB:11:22:33 or whatever. It should have scrolled by in the terminal above.
8. Type `trust $THAT_MAC_ADDRESS`, e.g. `trust AA:BB:11:22:33`.

Now when you run the bluetooth scanner (--scan), it should print both the MAC address AND the iPhone name. The iPhone "friendly" name can now also be used in the `--required_bluetooth_names` argument.

### Then you can run the `garage_door_opener.py` like this:

```shell
python3 garage_door_opener.py --code=23 --verbose --knock=50 --log=log.txt --dry --required_bluetooth_names=YourPhonePrefix,Christopher.*iPhone,Jane.*[Aa]pple.Watch --pause=1.0
```

### One last caveat

It turns out that the iPhone does not just constantly broadcast its presence when it doesn't have any specific reason to already be paired to the Raspberry Pi. So it's likely that authentication will fail even if the iPhone is nearby, if you haven't specifically paired with the Pi from the iPhone bluetooth settings page.
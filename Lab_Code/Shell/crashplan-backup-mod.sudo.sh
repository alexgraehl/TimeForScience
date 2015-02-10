#!/bin/bash


## Requires sudo access.
## Make sure this file is executable as root by:

## sudo chown root:staff crashplan-backup-mod.sudo.sh
## sudo chmod 2000 crashplan-backup-mod.sudo.sh
## sudo chmod g+rx crashplan-backup-mod.sudo.sh
## sudo chmod u+rwx crashplan-backup-mod.sudo.sh
## sudo chmod o-rwx crashplan-backup-mod.sudo.sh

if [[ $1 = "off" ]]; then
     echo "Turns off CrashPlan remote backups on the Mac."
     echo "Use the command <backon> to re-enable backups."
     echo ""
     launchctl unload /Library/LaunchDaemons/com.crashplan.engine.plist && echo "CrashPlan Mac remote backups have been disabled for now. Use <backon> to re-enable them. Restarting the computer may also re-enable them."
elif [[ $1 = "on" ]]; then
     launchctl load /Library/LaunchDaemons/com.crashplan.engine.plist && echo "CrashPlan Mac remote backups have been re-enabled."
else
    echo "You must pass in either the argument 'on' or 'off' to this script, to tell it whether to turn CrashPlan backups off or on."
fi


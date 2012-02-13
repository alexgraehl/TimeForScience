# -*-Sh-*- <-- tells emacs what kind of syntax highlighting to use


STATUSFILE="${HOME}/backup-status.txt"
echo ''
echo "About to check the backups..."


case "$HOSTNAME" in
    lighthouse)
	echo "Checking the backups on the machine LIGHTHOUSE..."
	ionice -c3 rsync --dry-run --delete -vr /db/ /it/lighthouse-backup/data/db/
	;;

    bueno)
	echo "Checking the backups on the machine BUENO..."
	ionice -c3 rsync --dry-run --delete -vr /home/ /it/bueno-backup/data/home/
	;;
    *)
	echo 'We do not know how to check backups on any machine except bueno and lighthouse...'
	;;
esac


echo 'Finished checking the backups'
echo '------------------------------------------'
echo '=========================================='
echo '------------------------------------------'


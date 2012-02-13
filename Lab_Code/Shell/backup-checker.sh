# -*-Sh-*- <-- tells emacs what kind of syntax highlighting to use


STATUSFILE="${HOME}/backup-status-${HOSTNAME}.txt"
CHECKSYNC=ionice -c3 rsync --dry-run --delete -vr

echo ''
echo "About to check the backups..."


case "$HOSTNAME" in
    lighthouse)
	echo "Checking the backups on the machine LIGHTHOUSE..."
	${CHECKSYNC} /db/ /it/lighthouse-backup/data/db/ > ${STATUSFILE}
	;;

    bueno)
	echo "Checking the backups on the machine BUENO..."
	${CHECKSYNC} /home/ /it/bueno-backup/data/home/ > ${STATUSFILE}
	;;
    *)
	echo 'We do not know how to check backups on any machine except bueno and lighthouse...'
	;;
esac


echo 'Finished checking the backups'
echo '------------------------------------------'
echo '=========================================='
echo '------------------------------------------'


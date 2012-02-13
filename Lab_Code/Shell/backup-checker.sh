# -*-Sh-*- <-- tells emacs what kind of syntax highlighting to use


STATUSFILE="${HOME}/backup-status-${HOSTNAME}.txt"
CHECKSYNC=sudo ionice -c3 rsync --dry-run --delete -vr

echo ''
echo "About to check the backups on the machine ${HOSTNAME}..."


case "$HOSTNAME" in
    lighthouse)
	${CHECKSYNC} /db/ /it/lighthouse-backup/data/db/ >> ${STATUSFILE}
	;;

    bueno)
	sudo ionice -c3 rsync --dry-run --delete -vr /home/ /it/bueno-backup/data/home/
	;;
    *)
	echo 'We do not know how to check backups on any machine except bueno and lighthouse...'
	;;
esac


echo 'Finished checking the backups'
echo '------------------------------------------'
echo '=========================================='
echo '------------------------------------------'


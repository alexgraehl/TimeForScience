# -*-Sh-*- <-- tells emacs what kind of syntax highlighting to use
#sudo ionice -c3 rsync --dry-run --delete -vr /home/ /it/bueno-backup/data/home/

#echo "ARARA"

STATUSFILE="${HOME}/backup-status-${HOSTNAME}.txt"
CHECKSYNC="sudo ionice -c3 rsync --dry-run --delete -vr"

echo ''
echo "About to check the backups on the machine ${HOSTNAME}..."
date
echo ''

echo '===================================' >> ${STATUSFILE}
echo 'Starting new backup check:' >> ${STATUSFILE}
date >> ${STATUSFILE}
echo '' >> ${STATUSFILE}

case "$HOSTNAME" in
    lighthouse)
	${CHECKSYNC} /db/ /it/lighthouse-backup/data/db/ >> ${STATUSFILE}
	;;
    bueno)
	${CHECKSYNC} /home/ /it/bueno-backup/data/home/ >> ${STATUSFILE}
	;;
    *)
	echo 'We do not know how to check backups on any machine except bueno and lighthouse...'
	;;
esac


echo 'Finished checking the backups'
echo '------------------------------------------'
echo '=========================================='
echo '------------------------------------------'


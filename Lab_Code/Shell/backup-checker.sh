# -*-Sh-*- <-- tells emacs what kind of syntax highlighting to use
#sudo ionice -c3 rsync --dry-run --delete -vr /home/ /it/bueno-backup/data/home/

#echo "ARARA"

STATUSFILE="${HOME}/backup-status-${HOSTNAME}.txt"
CHECKSYNC="sudo ionice -c3 rsync --dry-run --delete -vr \
 --exclude '*.tmp' \
 --cvs-exclude \
 --exclude '*~' \
 --exclude '*.bak' \
 --exclude '\.DS_Store' "

echo $NOSYMWARN

echo ''
echo "About to check the backups on the machine ${HOSTNAME}..."
echo "Note: REQUIRES SUDO ACCESS. You will have to type your password in order for rsync to run here."
date
echo ''

echo '===================================' >> ${STATUSFILE}
echo 'Starting new backup check:' >> ${STATUSFILE}
date >> ${STATUSFILE}
echo '' >> ${STATUSFILE}


case "$HOSTNAME" in
    lighthouse)
	${CHECKSYNC} --exclude "no_backup/*" --exclude "mirror_download/*" /db/ /it/lighthouse-backup/data/db/ | grep -v 'non-regular' >> ${STATUSFILE}
	;;
    bueno)
	${CHECKSYNC}  /home/ /it/bueno-backup/data/home/ | grep -v 'non-regular' >> ${STATUSFILE}
	;;
    *)
	echo 'We do not know how to check backups on any machine except bueno and lighthouse...'
	;;
esac


echo 'Finished checking the backups'
echo '------------------------------------------'
echo '=========================================='
echo '------------------------------------------'


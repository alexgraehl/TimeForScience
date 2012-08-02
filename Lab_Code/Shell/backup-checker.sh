# -*-Sh-*- <-- tells emacs what kind of syntax highlighting to use
#sudo ionice -c3 rsync --dry-run --delete -vr /home/ /it/bueno-backup/data/home/

#echo "ARARA"

HN=`hostname`

STATUSFILE="${HOME}/backup-status-${HN}.txt"
CHECKSYNC="sudo ionice -c3 rsync --dry-run --delete -vr \
 --exclude '*.tmp' \
 --cvs-exclude \
 --exclude '*~' \
 --exclude '*.bak' \
 --exclude '\.DS_Store' "

echo $NOSYMWARN

echo ''
echo "About to check the backups on the machine ${HN}..."
echo "Note: REQUIRES SUDO ACCESS. You will have to type your password in order for rsync to run here."
date
echo ''

echo '===================================' >> ${STATUSFILE}
echo 'Starting new backup check:' >> ${STATUSFILE}
date >> ${STATUSFILE}
echo '' >> ${STATUSFILE}


case "$HN" in
    lighthouse)
	${CHECKSYNC} --exclude "no_backup/*" --exclude "mirror_download/*" /db/ /it/${HN}-backup/data/db/ | grep -v 'non-regular' >> ${STATUSFILE}
	;;
    bueno)
	${CHECKSYNC}  /home/ /it/${HN}-backup/data/home/ | grep -v 'non-regular' >> ${STATUSFILE}
	;;
    nausicaa)
	${CHECKSYNC}  /home/ /it/${HN}-backup/data/home/ | grep -v 'non-regular' >> ${STATUSFILE}.home.txt
	${CHECKSYNC}  /work/ /it/${HN}-backup/data/work/ | grep -v 'non-regular' >> ${STATUSFILE}.work.txt
	;;
    catbus)
	${CHECKSYNC}  /bigdata/ /it/${HN}-backup/data/bigdata/ | grep -v 'non-regular' >> ${STATUSFILE}
	;;
    *)
	"We do not know how to check backups on any machine except <bueno>, <catbus>, <nausicaa> and <lighthouse>. This machine was <$HN>."
	;;
esac


echo 'Finished checking the backups'
echo '------------------------------------------'
echo '=========================================='
echo '------------------------------------------'


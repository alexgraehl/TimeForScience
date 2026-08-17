#!/usr/bin/env bash

# This is a script that was very specifically set up for RSync-ing bioinformatics files, as you can
# tell from the file exclusion list.

REMOTE_USERNAME=*************         # (<-- Fill this in)
WHICH_REMOTE_HOST=*****************   # (<-- Fill this in)
SRC_BACKUP_PATH=/work/Common/          # <-- Change this. Expeted to end in a slash
BKDIR=/Backup_Local                   # (<-- Fill this in)


echo ''
echo "About to copy remote code files from ${WHICH_REMOTE_HOST} to this local machine as a backup..."
echo ''
echo '------------------------------------------'
echo '=========================================='
echo '------------------------------------------'
echo ''


mkdir -p ${BKDIR}

DEL_ARG=
if [[ -z $1 ]]; then
	echo "Not deleting any local files."
	echo "If you want to delete local files that are not on the remote server, run this script with --delete."
	echo ''
else
	if [[ $1 == "--delete" ]] ; then
		echo ''
		echo "WE ARE GOING TO DELETE ANY LOCAL FILES IN ${BKDIR} THAT DO NOT EXIST ON THE SERVER!!! BEWARE!" ;
		echo ''
		DEL_ARG="--delete-after"
	fi
fi

echo 'ENTER YOUR PASSWORD WHEN PROMPTED BY RSYNC:'

rsync --verbose --progress --stats --compress --rsh=/usr/bin/ssh   \
		--human-readable --recursive --times               \
		--archive \
		--delete-excluded \
		--compress \
		--perms --links \
		--exclude "*~"  \
		--exclude "\.DS_Store"    \
		--cvs-exclude             \
		--exclude "apt[-_]out/*"  \
		--exclude "bowtieOut/"    \
		--exclude "chipDat/test"  \
		--exclude "test/Data"     \
		--exclude "*.bowtie"      \
		--exclude "*.wig"         \
		--exclude "*.RData"       \
		--exclude "*.RImage"      \
		--exclude "*.temp"      \
		--exclude "*.tmp"      \
		--exclude "tmp/*"	\
		--exclude "temp/*"	\
		--exclude "TempDir*"	\
		--exclude "*.bak"      \
		--exclude "*.CEL"   \
		--exclude "*.fq"      \
		--exclude "*.fastq"      \
		--exclude "*.aln"      \
		--exclude "*.bar"      \
		--exclude "*.sam" \
		--exclude "*.bam" \
		--exclude "binf-core-internal-*"   \
		--exclude "Data" \
		${DEL_ARG} \
		${REMOTE_USERNAME}@${WHICH_REMOTE_HOST}:${SRC_BACKUP_PATH} ${BKDIR}

echo ''
echo "Finished copying backups to <${BKDIR}>"
echo ''
echo '------------------------------------------'
echo '=========================================='
echo '------------------------------------------'


# --delete-after

#		--delete-after \

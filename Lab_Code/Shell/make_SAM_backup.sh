#!/bin/sh

# This script backs up all of the SAM_T05 materials to the location specified in
# the variable $ARCHIVE_FINAL_LOCATION

ARCHIVE_FINAL_LOCATION=./

CP_FLAGS="-vRp"
# R: recursive--must be upper-case to copy THROUGH symlinks instead of copying the links themselves
# p: preserve modification times
# v: verbose. Maybe want to turn this off if you don't like the scrolling text indicating progress.

BACKUP_TMP_FOLDER_PATH="/tmp/SAM-T05-project-backup" # where to put the temporary .tar file
BACKUP_ARCHIVE_NAME="SAM-T05-archive"

# Below: these are the folders to back up into an archive
SCRIPTS_FOLDER=/projects/compbio/experiments/protein-predict/SAM_T05test/
CGI_BIN_FOLDER=/cse/guests/farmer/.html/cgi-bin/SAM_T05test/
HTML_FOLDER=/projects/compbiousr/.html/SAM_T05test/

mkdir $BACKUP_TMP_FOLDER_PATH

tar -cvf ${BACKUP_TMP_FOLDER_PATH}/${BACKUP_ARCHIVE_NAME}.tar $SCRIPTS_FOLDER $CGI_BIN_FOLDER $HTML_FOLDER
gzip     ${BACKUP_TMP_FOLDER_PATH}/${BACKUP_ARCHIVE_NAME}.tar

mv -i    ${BACKUP_TMP_FOLDER_PATH}/${BACKUP_ARCHIVE_NAME}.tar.gz ${ARCHIVE_FINAL_LOCATION}

rm -rf ${BACKUP_TMP_FOLDER_PATH}



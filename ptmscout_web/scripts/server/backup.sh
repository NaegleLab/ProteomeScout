#!/bin/bash

set -e

DB_BACKUP="/data/ptmscout/ptmscout_web/data/backup/ptmscout.`date +%Y-%m-%d`.sql"
DB_USER="root"
DB_PASS=""

BACKUP_1=""
BACKUP_2=""

mysqldump --user=$DB_USER --password=$DB_PASS ptmscout > $DB_BACKUP
echo "put $DB_BACKUP" | sftp $BACKUP_1
echo "put $DB_BACKUP" | sftp $BACKUP_2

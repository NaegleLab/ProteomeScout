#!/bin/bash

set -e

echo "Enter database username:"
read username

mysqldump -u $username -p --no-create-db --no-data ptmscout_test > scripts/database/database_defs.sql
mysqldump -u $username -p --no-create-db --no-create-info ptmscout PTM PTM_keywords PTM_taxonomy > scripts/database/ptm_defs.sql
mysqldump -u $username -p --no-create-db --no-create-info ptmscout taxonomy > scripts/database/ncbi_taxonomy.sql
mysqldump -u $username -p ptmscout_test > scripts/database/ptmscout_test.sql

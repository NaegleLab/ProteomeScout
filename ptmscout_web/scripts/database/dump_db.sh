#!/bin/bash

mysqldump -u mkmatlock -p --no-create-db --no-data ptmscout > scripts/database/database_defs.sql
mysqldump -u mkmatlock -p --no-create-db --no-create-info ptmscout PTM PTM_keywords PTM_taxonomy > scripts/database/ptm_defs.sql
mysqldump -u mkmatlock -p --no-create-db --no-create-info ptmscout taxonomy > scripts/database/ncbi_taxonomy.sql
mysqldump -u mkmatlock -p ptmscout_test > scripts/database/ptmscout_test.sql

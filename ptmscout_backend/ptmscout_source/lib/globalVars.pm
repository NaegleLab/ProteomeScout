package globalVars;

use strict;
use warnings;

BEGIN{
    require Exporter;
    @globals::testGlobals::ISA = qw(Exporter);
    @globals::testGlobals::EXPROT = qw();

    our ($GO_PATH, $MYSQL_USER, $MYSQL_PWD, $MYSQL_DB, $PFAM_PATH, $LOG_PATH, $REFSEQ_DB, $REFSEQ_PATH, $BACKUP_PATH, $TESTDB_USER, $TESTDB_PWD, $TEST_DB, $DATA_PATH, $MYSQL_UPDATE, @ACC_TYPES, $CPU_DAY, $CPU_NIGHT);
    
    $GO_PATH = '<your-path>/ptmscout_backend/GO/';
    $PFAM_PATH = "<your-path>/ptmscout_backend/pfam/";
    $DATA_PATH = "<your-path>/ptmscout_backend/datasets/";
    
    $MYSQL_USER = ''; #name of your production database account, must have write and read accessibility
    $MYSQL_DB = ''; #name of database, primary here - "production"
    $MYSQL_PWD = ''; #password here

    #$OTHER_DB = ''; #here you can add as many alternate db names..make sure to export and configure DBTools::dbIO to reference a short code. 

    $TEST_DB = 'test'; #if you want a test database here
    $TESTDB_USER = '';
    $TESTDB_PWD = ''; #put password here

    $MYSQL_UPDATE = ''; #this is a good idea for when you want to make major updates to domains, GO terms, etc... create an update database (same user and pwd as production). Once filled, copy this over production for rollout.
    
    $REFSEQ_DB = 'refseq'; #make this the database for your local refseq
    $REFSEQ_PATH = '<your-path>/ptmscout_backend/NCBI_DB/';

    $LOG_PATH = "<your-path>/ptmscout_backend/log";

    $BACKUP_PATH = "<your-path>/ptmscout_backend/backups/";

    $CPU_DAY = 6; #number of processors to use for threaded hmmpfam on weekdays between 8am and 4pm
    $CPU_NIGHT = 8; #numbe of processors to use at all other times. 
    # change weekday hour definitions in commonTools.pm

    @ACC_TYPES = ('refseq', 'UniprotKB', 'GenPept', 'IPI'); 

    @globals::testGlobals::EXPORT_OK = qw($GO_PATH $PFAM_PATH $MYSQL_USER $MYSQL_DB $MYSQL_PWD $LOG_PATH $REFSEQ_DB $REFSEQ_PATH $BACKUP_PATH $TEST_DB $TESTDB_USER $TESTDB_PWD $DATA_PATH $MYSQL_UPDATE \@ACC_TYPES $CPU_DAY $CPU_NIGHT);
#    print "Hi there\n";
}
END {}

1;

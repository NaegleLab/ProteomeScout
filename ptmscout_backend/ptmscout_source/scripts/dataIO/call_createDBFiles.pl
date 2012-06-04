#! /usr/bin/perl 
use dataIO;
use DBTools::dbIO;


if(scalar(@ARGV < 3)){
    print "Usage: call_loadDBFiles.pl DB_CODE DIR \@EXPS\n";
    exit;
}
my ($DB, $dir, @exps) = @ARGV;
my $dbh = returnDBHNOCommit($DB);

createDBFiles($dbh, $dir, \@exps);

$dbh->rollback();

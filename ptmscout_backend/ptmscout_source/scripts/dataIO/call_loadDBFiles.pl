#! /usr/bin/perl
use dataIO;
#use DBI;
use DBTools::dbIO;

if(scalar(@ARGV !=2)){
    print "Usage: call_loadDBFiles.pl DB_CODE DIR\n";
    exit;
}
my ($DB, $dir) = @ARGV;
my $dbh = returnDBHNOCommit($DB);


loadDBFiles($dbh, $dir);

$dbh->commit();

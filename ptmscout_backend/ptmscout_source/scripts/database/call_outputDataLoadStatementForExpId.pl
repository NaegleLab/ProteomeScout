#! /usr/bin/perl
use strict;
use warnings;
use DBTools::dbIO;
use DBTools::dbEntry;

if(scalar(@ARGV) != 2){
    print "Usage call_outputDataLoadStatementForExpId.pl EXP_ID DB\n";
    exit;
}
my ($expId, $DB) = @ARGV;
my $dbh=returnDBHNOCommit($DB);

my $str = outputDataLoadStatementForExpId($dbh, $expId);
print $str."\n";

$dbh->rollback();

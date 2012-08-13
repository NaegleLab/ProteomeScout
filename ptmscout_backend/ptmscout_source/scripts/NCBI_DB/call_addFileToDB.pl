#! /usr/bin/perl
use strict;
use warnings;
use ncbiDBTools;
use DBTools::dbIO;

if(scalar(@ARGV) != 1){

    print "Usage: call_getNCBIFile.pl FILE_NUMBER\n";
    exit;
}
my $number = $ARGV[0];
getNCBIFile($number);

my $dbh = returnRefSeqDBHNOCommit();
my $date = 'March 2008';
my $fileName = returnFileNameMV($number);
addFileToDB($dbh, $fileName, $date, 0);
$dbh->commit();
#$dbh->rollback();


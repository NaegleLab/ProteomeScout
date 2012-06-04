#! /usr/bin/perl
use strict;
use warnings;
use ncbiDBTools;
use DBTools::dbIO;

if(scalar(@ARGV) != 1){
    print "Usage: call_getNCBIFile.pl BOOL_CHECK_FOR_GI\n";
    exit;

}


my $CHECK = $ARGV[0];
my $dbh = returnRefSeqDBHNOCommit();


addAllFilesToDB($dbh, $CHECK);

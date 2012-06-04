#! /usr/bin/perl
use strict;
use warnings;
use DBTools::dbEntry;
use DBTools::dbIO;

if (scalar(@ARGV) < 1){
    print "Usage: call_hanldeProteinByAcc.pl ACC\n";
    print "INTO TEST DATABASE\n";
    exit;

}

my $dbh = returnTestDBH();
my $acc = $ARGV[0];
my($proteinId, $iProteinFlag, $iAccFlag);

($proteinId, $iProteinFlag, $iAccFlag) = handleProteinByAcc($dbh, $acc);
if($proteinId && $iProteinFlag && $iAccFlag){
    print "PASSED\n";
}
else{
    print "FAILED\n";
    print "ProteinID = $proteinId\n INSERT PROTEIN: $iProteinFlag\n INSERT ACC: $iAccFlag\n";
}


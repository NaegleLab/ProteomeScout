#! /usr/bin/perl
use strict;
use warnings;
use DBTools::dbIO;
use DBTools::maintenance;


if(scalar(@ARGV) != 2){
    print "USAGE: call_maintainAmbiguousExperiments.pl LOCAL REFSEQ\n";
    print "Maintains all experiments with ambiguity bit by force, comparing to just LOCAL or REFSEQ or both\n";
    exit();
}

my ($LOCAL, $REFSEQ)= @ARGV;
if(!$LOCAL and !$REFSEQ){
    print "ERROR! You must either insert from one or both databases, both arguments cannot be 0\n";
    exit();
}

my $dbh = returnProductionDBHNOCommit();
my $rdbh = returnRefSeqDBHNOCommit();


maintainAmbiguousExperiments($dbh, $rdbh, $LOCAL, $REFSEQ);

$dbh->commit();

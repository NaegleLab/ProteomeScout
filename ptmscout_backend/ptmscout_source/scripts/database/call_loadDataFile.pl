#! /usr/bin/perl 
use strict;
use warnings;
use DBTools::dbIO;
use DBTools::dbEntry;
use errorHandling;
use HMMTools;

if(scalar(@ARGV) != 18 ){
    print "Usage: call_loadDataFile.pl DATA_FILE NAME AUTHOR DESCRIPTION CONTACT PMID URL PUBLISHED AMBGUITY EXPORT EXPID_LINK SUBMISSION_EMAIL PRIMARY_MOD(S) JOURNAL PUB_DATE VOLUME PAGES [P|U|T]\n";
    print "WHERE database argument: P-production, U-update, T-test\n";
    exit;
}

my $dataFile = $ARGV[0];
my $name = $ARGV[1];
my $author = $ARGV[2];
my $description = $ARGV[3];
my $contact = $ARGV[4];
my $PMID = $ARGV[5];
my $URL = $ARGV[6];
my $published = $ARGV[7];
my $AMBIGUITY = $ARGV[8];
my $EXPORT = $ARGV[9];
my $expId_link = $ARGV[10];
my $submissionEmail = $ARGV[11];
my $primaryMods = $ARGV[12];
my $journal = $ARGV[13];
my $pub_date = $ARGV[14];
my $volume = $ARGV[15];
my $pages = $ARGV[16];
my $DB = $ARGV[17];

my $dbh = returnDBHNOCommit($DB);

#`perl -i~ -pe 's/\cM/\n/g' $dataFile`;

#xsmy $dbh = returnWebDevDBHNOCommit();
#my $dbh = returnTestDBHNOCommit();
flushLog();

my ($expId, $accepted, $rejected);

($expId, $accepted, $rejected) = call_loadDataFile($dbh, $dataFile,$name, $author, $description, $contact, $PMID, $URL, $published, $AMBIGUITY, $EXPORT, $expId_link, $submissionEmail, $primaryMods, $journal, $pub_date, $volume, $pages);
print "Experiment ID: $expId\n";
print "Accepted: $accepted Rejected:$rejected\n";



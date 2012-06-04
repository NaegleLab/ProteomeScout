#! /usr/bin/perl
use strict;
use warnings;
use DBTools::dbIO;
use DBTools::maintenance;


my $dbh = returnUpdateDBHNOCommit();

createNewVersion($dbh);

$dbh->rollback();

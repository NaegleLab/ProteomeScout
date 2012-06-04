#! /usr/bin/perl
use strict;
use warnings;
use DBTools::dbIO;
use DBTools::deleteFunctions;

my $dbh = returnProductionDBHNOCommit();

# find all proteins with no domains - then delete records
my $sth = $dbh->prepare('SELECT protein.id from protein left outer join domain on protein.id=domain.protein_id where domain.protein_id IS NULL');
$sth->execute();

my $proteinIds = returnArrayOfResultsOnCol($sth,0);

print "Found ".scalar(@$proteinIds)." protein ids to delete\n";
foreach my $proteinId (@$proteinIds){
    my $error = deleteProteinRecordsByProteinId($dbh, $proteinId);
    if($error){
	print "ERROR deleting $proteinId\n";
    }

}

$dbh->commit();

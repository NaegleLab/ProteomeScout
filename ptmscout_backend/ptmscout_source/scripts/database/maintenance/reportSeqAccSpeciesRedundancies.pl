#! /usr/bin/perl
use strict;
use warnings;
use DBTools::dbIO;
use DBTools::queryTools;

my $dbh = returnProductionDBHNOCommit();

# get sequences that are redundant
my $sth = $dbh->prepare('SELECT sequence, species, count(*) as NUM from protein group by species, sequence having NUM > 1');
$sth->execute();
my $sth2=$dbh->prepare('SELECT * from protein where sequence=? and species=?');
my $count = 0;
my $sth3 = $dbh->prepare('SELECT * from protein where sequence=? and species=? and acc_gene=?');
my %redundant;
while(my @row = $sth->fetchrow_array()){
    $count += 1;
    my $sequence = $row[0];
    my $species = $row[1];
    #now get protein info for that sequence
    $sth2->execute($sequence, $species);
    my $results = returnMultipleResultsForExSTH($sth2);
    foreach my $result (@$results){
	my $gene; 
	$gene = $result->{'acc_gene'};
	# look up to see if there's redundancy in sequence, species, and acc_gene
	$sth3->execute($sequence, $species, $gene);
	my $results3 = returnMultipleResultsForExSTH($sth3);
	if(scalar(@$results3) > 1){
	    my $key = $gene."_".$species;
	    
	   # print "$gene\t$species\n";
	    if(not defined $redundant{$key}){
		$redundant{$key} = [];
		foreach my $r (@$results3){
		    push @{$redundant{$key}}, $r->{'id'};
#		    print "\t $r->{'id'}\n";
		}
	    }
	    
	}

    }
}

foreach my $k (keys %redundant){

    print "$k\n";
    foreach my $id (@{$redundant{$k}}){
	print "\t$id\n";
    }
}

$dbh->rollback();


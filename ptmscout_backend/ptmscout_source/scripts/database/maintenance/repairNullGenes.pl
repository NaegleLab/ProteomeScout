#! /usr/bin/perl
use strict;
use warnings;
use DBTools::dbIO;
use DBTools::dbEntry;
use DBTools::resultParsing;
use entrezTools;

sub handleRichSeqHashObj($$);
sub repeatInBatches($$$);

my $dbh=returnProductionDBHNOCommit();
my $log = "repairNullGenes.log";
open(LOG, ">$log") || die "Can't open $log for writing\n";
print LOG "Protein Id\tAcc\tWarning\n";

# get all protein ids with NULL gene names
my $sth = $dbh->prepare('SELECT id from protein where acc_gene is NULL or acc_gene=?');
$sth->execute('NULL');
my $proteinIds = returnArrayOfResultsOnCol($sth, 0);

print "Found ",scalar(@$proteinIds)," proteins with null gene names\n";

my @swiss;
my @gi;
foreach my $proteinId (@$proteinIds){
    my $accArr = returnAccValuesByProteinId($dbh, $proteinId);
    my $accHash = returnAccHashByTypeFromArr($accArr);
    if(defined($accHash->{'swissprot'})){
	my $swissArr = $accHash->{'swissprot'};
	push @swiss, $swissArr->[0];
    }
    elsif(defined($accHash->{'gi'})){
	my $giArr = $accHash->{'gi'};
	push @gi, $giArr->[0];
    }
}

# now - go through and get gene and gene synonyms 
#for swiss: - - go through to change this to handle batches of size 100 (like in other)



#my $richSeqHash = returnGenPeptQueryByBatch(\@test);
#handleRichSeqHashObj($richSeqHash, \@test);
repeatInBatches($dbh, 1, \@swiss);
repeatInBatches($dbh, 1, \@gi);

    
print "Check $log for warnings\n";
close LOG;
$dbh->commit();


# break all accessions into batches of 100 and handle
sub repeatInBatches($$$){
    my ($dbh, $COMMIT, $missAcc) = @_;

    my $start = 0;
    my $stop;
    my $groupSize = 100;

    print "Retreiving ".scalar(@$missAcc)." in groups of $groupSize\n";    
    while(scalar(@$missAcc)){
	if(scalar(@$missAcc) < $groupSize){
	    $stop = scalar(@$missAcc) - 1;
	}
	else{
	    $stop = $groupSize - 1;
	   
	}
	my @miss = @$missAcc[$start..$stop];
	splice(@$missAcc, $start, $stop-$start+1);
	my $richSeqHash = returnGenPeptQueryByBatch(\@miss);
	handleRichSeqHashObj($richSeqHash, \@miss);
	if($COMMIT){
	    $dbh->commit();
	}
    }
    if($COMMIT){
	$dbh->commit();
    }
}



# This is the sub that goes through each accession and the richSeqHash to handle inserting gene accessions etc. 
sub handleRichSeqHashObj($$){
    my ($richSeqHash, $accArr) = @_;
    
    foreach my $acc (@$accArr){

	my $richSeq = $richSeqHash->{$acc};
	#get proteinId for acc
	my $proteinId = returnProteinIdByAcc($dbh, $acc);
	my ($gene, $gene_syn) = returnGeneGeneSynNames($richSeq);
	print "For ProteinId: $proteinId  accession:$acc has gene: $gene and gene_synonyms: @$gene_syn\n";
	my ($WARN, $warnTxt) = checkForWarnings($richSeq);
	if($WARN){
	    print LOG "$proteinId\t$acc\t $warnTxt\n";
	}
	my $insertedAccIds = handleGeneSynonymsForProteinId($dbh, $gene_syn, $proteinId);
	if(@$insertedAccIds){
	    print "Inserted @$insertedAccIds\n";
	}
	if($gene){
	    #print "Want to put $gene in for $proteinId\n";
	    updateAccGeneFieldForProteinId($dbh, $proteinId, $gene);
	}
	
    }
}

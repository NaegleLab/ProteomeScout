#! /usr/bin/perl
use strict;
use warnings;
use DBTools::dbIO;
use entrezTools;
use DBTools::resultParsing; 
use DBTools::insertFunctions;
#use SBMLTools;


# \@probesetIds = findProbeSetIdForProteinId($dbh, $protein_id)
# Returns the probeset ids for a protein id, looks for swissprot and refseq accessions and then finally acc_gene if none of those link to the annotation table
# Inputs: $dbh - database handle
#         $proteinId - id to protein table
# Outputs: \@probesetIds - ref. to array of probeset ids
# Kristen Naegle
# August 11, 2008
# Updated January 15, 2009 to account for species!
sub findProbesetIdForProteinId($$){
    my ($dbh, $proteinId) = @_;
 
    my @probesetIds;
    #Get species info
    my $proteinHash = returnProteinDescForProteinId($dbh, $proteinId);
    
    my $species = $proteinHash->{'species'};
    if($species ne 'homo sapiens' and $species ne 'mus musculus'){
	my @x;
	push @x, -1;
	return \@x;

    }
    #if mouse, only look in the gene accession field and aliases field  (do it for both, then skip remaining steps if mouse);  

    # ADD GENE Synononym search and don't look if acc_gene is '' 
    #my $acc_gene = $proteinHash->{'acc_gene'};
    my $genes = returnGeneListForProteinId($dbh, $proteinId);
    
    foreach my $gene (@$genes){

	my $ALIAS_SEARCH = 0;
	my $ids = returnProbesetIdByGene($dbh, $gene, $species, $ALIAS_SEARCH);
	#  print "Testing $acc_gene for $species\n";
	if($ids->[0] ne '-1'){
	    push @probesetIds, @$ids;		
	    #print "found @probesetIds\n";
	}
    }
    if($species eq 'homo sapiens'){
	my $accArr = returnAccValuesByProteinId($dbh, $proteinId);
	
	foreach my $acc (@$accArr){
	    my $type = returnAccType($acc);
	    if($type eq 'swissprot'){
# 		#  print "Found swissprot type acc\n";
		my $ids = returnProbesetIdByUniprot($dbh, $acc);
 		if($ids->[0] ne '-1'){
 		    push @probesetIds, @$ids;		
# 		    #print "found @probesetIds\n";
 		}
	    }
	    elsif($type eq 'entrez_protein' or $type eq 'refseq'){
		my $ids = returnProbesetIdByRefSeq($dbh, $acc);
		if($ids->[0] ne '-1'){
		    push @probesetIds, @$ids;		
# 		    #print "found @probesetIds\n";
		}
	    }
	    
	}
    }
#handle negative case 
    if(scalar(@probesetIds) == 0){
	push @probesetIds, -1;
    }
    

    # uniquify the probeset;
    my $hashRef = createUniqHash(\@probesetIds);
    my @ids = keys %$hashRef;
    return \@ids;
}

# ($numLinks, $numInserted) = handleExpressionForProteinId($dbh, $proteinId)
# handles input into protein_expression table for a particular protein id
# Inputs: $dbh - database handle
#         $proteinId - id to protein table
# Outptus: $numLinks - total number of links to annotation table found for protein
#          $numInserted - number of those that were inserted (the others already existed)
# Kristen Naegle
# August 12, 2008
sub handleExpressionForProteinId($$){
    my ($dbh, $proteinId) = @_;
    my $probesetIds = findProbesetIdForProteinId($dbh, $proteinId);
    my $numFound = 0;
    my $numInserted =0;
    if($probesetIds->[0] ne '-1'){
	foreach my $probeId (@$probesetIds){
	    $numFound+=1;
	    my %args;
	    $args{'protein_id'} = $proteinId;
	    $args{'probeset_id'} = "'".$probeId."'";
	    my $exist = checkForExistence($dbh, '*', 'protein_expression', \%args);
	    if(!$exist){
		$numInserted += 1;
		my $id = insertProteinExpression($dbh, $proteinId, $probeId);
		# print "Inserted: id=$id\n";
	    }
	    
	}
    }
    return($numFound, $numInserted);
}

# ($numProteinsHandled, $numFailed) = maintainExpressionLinkages($dbh);
# maintain linkages between proteins and expression
# find all protein ids not in the protein_expression table and handle them
# Inputs: $dbh - database handle
# Outputs: $numProteinsHandled - number of proteins that were attempted to be linked
#          $numFailed - number that could not be
# Kristen Naegle
# August 12, 2008
sub maintainExpressionLinkages($){
    my ($dbh) = @_;
    my $proteinIds = returnProteinsWithoutExpression($dbh); #only mouse and human
    my $numFailed = 0;
    my $numProteinsHandled = scalar(@$proteinIds);
    foreach my $proteinId (@$proteinIds){
	my ($numFound, $numInserted) = handleExpressionForProteinId($dbh, $proteinId);
	if($numFound==0){
	    $numFailed += 1;
	}
    }
    return ($numProteinsHandled, $numFailed);

}



# \@probeset_ids = returnProbesetIdByUniprot($dbh, $uniprot)
# Does a regexp search against the uniprot field of expression_ann to find the corresponding probeset_ids
# Inputs: $dbh - database handle
#         $uniprot - uniprot id
# Outputs: \@probeset_ids - ref to array of probeset ids, first will be -1 if not found
# Kristen Naegle
# August 11, 2008
# Updated January 16, 2009 to include word boundaries
sub returnProbesetIdByUniprot($$){
    my ($dbh, $uniprot) = @_;

    my $sth=$dbh->prepare('SELECT probeset_id from expression_ann where uniprot REGEXP ?');
    my $reg = '[[:<:]]'.$uniprot.'[[:>:]]'; #add word boundaries..these fields are words
    $sth->execute($reg);
    my $probeset_ids = returnArrayOfResultsOnCol($sth, 0);
    return $probeset_ids;
}

# \@probeset_ids = returnProbesetIdByRefSeq($dbh, $refseq)
# Does a regexp search against the refseq field of expression_ann to find the corresponding probeset_ids
# Inputs: $dbh - database handle
#         $refseq - refseq id
# Outputs: \@probeset_ids - ref to array of probeset ids, first will be -1 if not found
# Kristen Naegle
# August 11, 2008
# Updated January 16, 2009 to include word boundaries
sub returnProbesetIdByRefSeq($$){
    my ($dbh, $refseq) = @_;

    my $sth=$dbh->prepare('SELECT probeset_id from expression_ann where refseq REGEXP ?');
    my $reg = '[[:<:]]'.$refseq.'[[:>:]]'; #add word boundaries..these fields are words
    $sth->execute($reg);
    my $probeset_ids = returnArrayOfResultsOnCol($sth, 0);
    return $probeset_ids;
}

# \@probeset_ids = returnProbesetIdByGene($dbh, $gene_name, $species, $ALIAS_SEARCH)
# Does two searches: exact search on symbol field and regexp search on the aliase field (if ALIAS_SEARCH is true)
# Inputs: $dbh - database handle
#         $gene_name - acc_gene in protein field and name in expression tables
#         $species - species
#         $ALIAS_SEARCH - if 1, look in the alias field as well (only if acc_gene -> symbol field returned nothing), else ignore alias. 
# Outputs: \@probeset_ids - ref to array of probeset ids, first will be -1 if not found
# Kristen Naegle
# August 11, 2008
sub returnProbesetIdByGene($$$$){
    my ($dbh, $gene_name, $species, $ALIAS_SEARCH) = @_;
    my @probeset_ids;
    my $sth=$dbh->prepare('SELECT probeset_id from expression_ann where symbol=? and species=?');
    $sth->execute($gene_name, $species);
    my $symbol_ids = returnArrayOfResultsOnCol($sth, 0);
    if($symbol_ids->[0] ne '-1'){
	push @probeset_ids, @$symbol_ids;

    }
    if($ALIAS_SEARCH and scalar(@probeset_ids)==0){
	my $sth2 = $dbh->prepare('SELECT probeset_id from expression_ann where species=? and aliases REGEXP ?');
	$sth2->execute($species, $gene_name);
	my $alias_ids = returnArrayOfResultsOnCol($sth2, 0); #hm, if we found something by alias..should we extract the gene and go back and look for gene with ALIAS_SEARCH off?
	if($alias_ids->[0] ne '-1'){
	    push @probeset_ids, @$alias_ids;
	}
    }
    if(scalar(@probeset_ids)==0){
	push @probeset_ids, -1;
    }

    return \@probeset_ids;
}

1;

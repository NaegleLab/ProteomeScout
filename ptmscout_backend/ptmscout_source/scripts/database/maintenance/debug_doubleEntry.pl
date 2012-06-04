#! /usr/bin/perl 
use strict;
use warnings;
use DBTools::dbIO;
use DBTools::dbEntry;
use entrezTools;


my $dbh = returnProductionDBHNOCommit();

my $acc = "gi|27370424";

my $richSeq = returnRichSeqFromAcc($acc);
my $date = '01-2008';



my ($errorCode, $sequence, $species, $gene, $name, $primaryAcc) = getProteinFromRichSeq($richSeq);

my $proteinId_fromProtein = returnProteinIdBySeqGene($dbh, $sequence, $gene, $species);

print "Protein Id from lookup: $proteinId_fromProtein\n";
print "Parmas: acc_gene=$gene, species=$species\n";


#my ($proteinId, $iProteinFlag, $iAccFlag) = handleProteinByRichSeq($dbh, $acc, $richSeq, $date);
#print "protein id: $proteinId\tiProteinFlag:$iProteinFlag\tiAccFlat:$iAccFlag\n";


$dbh->rollback();

use strict;
use warnings;
use LWP::Simple;
use Bio::SeqIO;
use Bio::DB::GenPept;
use Bio::DB::GenBank;
use pELMTools; 
use fileTools;


# performEntrezQuery($query, $outputFile);
# Retrieve protein query from Entrez and write to output file so that it can be parsed and used
# Output file is in asn.1 format 
#
# Kristen Naegle 5/22/07 
# Based on code written by Paul Stothard, Canadian Bioinformatics Help Desk
# based on a script by Oleg Khovayko http://olegh.spedia.net
#
# This script uses NCBI's Entrez Programming Utilities to perform
# batch requests to NCBI Entrez.
#
# Values in the 'edit these' section should be changed to meet your requirements
#
# See 'Entrez Programming Utilities' for more info at
# http://eutils.ncbi.nlm.nih.gov/entrez/query/static/eutils_help.html
sub performEntrezQuery($$){
    my $query = shift;
    my $outputFile = shift;

    my $db = "protein";
    my $report = "fasta";
    
    my $url = "http://www.ncbi.nlm.nih.gov/entrez/eutils";
    my $esearch = "$url/esearch.fcgi?" . "db=$db&retmax=1&usehistory=y&term=";

    my $esearch_result = get($esearch . $query);

    $esearch_result =~  m/<Count>(\d+)<\/Count>.*<QueryKey>(\d+)<\/QueryKey>.*<WebEnv>(\S+)<\/WebEnv>/s;

    my $count = $1;
    my $query_key = $2;
    my $web_env = $3;

    my $retstart;
    my $retmax = 5;

    open (OUTFILE, ">" . $outputFile) or die ("Error: Cannot open $outputFile : $!");

    print "$count entries to retrieve\n";

    for ($retstart = 0; $retstart < $count; $retstart = $retstart + $retmax) {
	print "Requesting entries $retstart to " . ($retstart + $retmax) . "\n";
	my $efetch = "$url/efetch.fcgi?" . "rettype=$report&retmode=text&retstart=$retstart&retmax=$retmax&" . "db=$db&query_key=$query_key&WebEnv=$web_env";   

	my $efetch_result = get($efetch);

	print (OUTFILE $efetch_result);
	sleep(3);
    }

    close (OUTFILE) or die( "Error: Cannot close $outputFile file: $!");
}

# $hashRichSeq = returnGenPeptQueryByBatch($accArr)
# Returns a hash of rich seq objects by batch query to entrez, key is accession
# Inputs: $accArr - reference to array of accessions (ASSUMES All accessions are of the same type, has to handle swissprot and gi differenly) PLEASE call this with blocks of less than 500 per batch
# Outputs: $hashRichSeq - ref to hash of rich seq objects, key is accession if bad accession, returns 0 instead of a richSeq object so check for if($hash->{$acc})
# Kristen Naegle
# March 24, 2008
sub returnGenPeptQueryByBatch($){
    my ($accArr) = @_;

    # uniquify to make sure we aren't retrieving something more than once.
    #$accArr = uniq($accArr);
    print "Size of Accession ".scalar(@$accArr)."\n";

    my %hash;
    foreach my $a (@$accArr){
	$hash{$a} = -1;
    }

    my $mode = 'batch';
    my $gb = Bio::DB::GenPept->new();
    my %params = $gb->get_params($mode);
    my $seqio = $gb->get_Stream_by_id($accArr);
    my $index = 0;
    my $type = returnAccType($accArr->[0]);

    my $seq;
    while(defined($seq=$seqio->next_seq)){
	my $acc;
	if($type eq 'gi'){
	    $acc = $seq->primary_id();
	    $acc = "gi|".$acc;
	}
	else{
	    $acc = $seq->display_id();
	}
	
	if($hash{$acc}){
	    if($seq){
		$hash{$acc} = $seq;		
	    }
	    else{
		$hash{$acc} = 0;		
	    }
	}
	else{
	    handleError('returnGenPeptQueryByBatch', "Accession $acc not of orignial type, can't match stream return to accession", \@_);
	}
    }
    my @failedArr;
    # Do single batch queries for any accession that is still set as -1;
    foreach my $acc (keys %hash){
	if($hash{$acc} == -1){
	    push @failedArr, $acc;
	    delete $hash{$acc};
	    handleError('returnGenPeptQueryByBatch',"Could not assing $acc a seqio object in batch query",\@_);
# # my $seq = returnRichSeqFromAcc($acc);
# # 	    if($seq){
# # 		$hash{$acc} = $seq;
# # 	    }
# # 	    else{
# # 		$hash{$seq} = 0;
# # 	    }
	}

    }

    return \%hash;
   

}



# getStuffAboutGI($gi)
# Given a GI, return stuff -- this is really more to record the key points of sequence objects returned from GenPept so I can copy and paste in the future. 
sub getStuffAboutGI($){
    my $query = shift;
    my $gb = new Bio::DB::GenPept();
    my $richSeq = $gb->get_Seq_by_id($query);
    
#if record doesn't exist then richSeq is empty
    if(!$richSeq){
	print("ERROR no id\n");
    }
# print accession number 
    my $primAcc = $richSeq->accession();
    print("Primary accession: $primAcc;\n");
    my $seq = $richSeq->seq();
# #my ($gi, $acc, $locus);
# #(undef, $gi, undef, $acc, $locus) = split(/\|/, $richSeq->display_id);
    print("Sequence: $seq\n");
    #print("GI: $gi\t ACC: $acc\t LOCUS: $locus\n");
    my $displayID = $richSeq->display_id;
    print ("Display ID: $displayID\n");
    

# Key words
    my @keywords = $richSeq->get_keywords();
    print("Keywords: @keywords\n");

#accessions
    my @acc = $richSeq->get_secondary_accessions();
    print("Accessions: @acc\n");

# # annotation types: 'comment', 'segment' 'origin' 'reference' 'date_changed' 'keyword' 'secondary_accession' 'dblink'
# # Get annotations of Seq_object 
    my $ann = $richSeq->annotation(); #annotation object
    print("$ann\n");
    
    foreach my $ref ($ann->get_Annotations('reference')){
	print("Reference: ", $ref->title,"\n");
    }

# Name
    my $name = $richSeq->display_name();
    print("Name: $name\n");
    my $namespace = $richSeq->namespace();
    print("Namespace: $namespace\n");

# Species - creates a species object
    my $species = $richSeq->species();
    my $spec = $species->species();
    my $organelle = $species->organelle();
    my $commonName = $species->common_name();
    print("Species Common: $commonName\n");
    print("Species Scientific: $spec\n");
    print("Species Organelle: $organelle\n");

#RNA/DNA or protein sequence?    
    my $alphabet = $richSeq->alphabet;
    print("Sequence alphabet is: $alphabet\n");

# #features
# #   my $numFeats = $richSeq->feature_count();
# #    my @feat = $richSeq->get_all_SeqFeatures(); #an array of hashes
# #    my %hash = $feat[0];
# #    print("Features: @feat\n");
# #    print("Feature count: $numFeats\n");
# #    foreach my $key (keys %hash){
# #	print "$key\n";
# #    }

# xrefs?
    my @features = $richSeq->get_SeqFeatures;
    my $location_type;
    my $feature_type;
    my $feature_tag; 
    my %seen_tag = ();
    my @tags = ();
    foreach my $feature(@features){
	$location_type = $feature->location->location_type;
	$feature_type = $feature->primary_tag;
	# #$feature_tag = $feature->tag;
# #	print "Feature type: $feature_type\n";
	foreach my $tag (@tags = $feature->get_all_tags){
	    $seen_tag{$tag}++;
	    #print("Tag: $tag\n");
	    if($tag eq "db_xref"){
		print("DB_XREF: $tag", $feature->get_tag_values('db_xref'),"\n"); 
	    }
	    #can I get the gene name?
	    if($tag eq "gene"){
		print("GENE: ", $feature->get_tag_values('gene'), "\n");
	    }
	}
	#Ah-hah so this code didn't work because db_xref was showing up as a tag vs. primary tag, but still these tags are in the feature area...we want DBSOURCE tags;
	if($feature_type eq "source"){
	    my @db_xref = $feature->get_tag_values('db_xref') if exists $seen_tag{'map'};
	    print "db_xref: @db_xref\n";
	}
	
     

    }

#xrefs a different way?
    my $ac = $richSeq->annotation;
    my @values = $ac->get_Annotations();
    my %tagname_type = map{$_->as_text, ($_->tagname . " " . ref($_))} @values;
    print("Values: @values\n");
   # for my $value (@values){
	#print("DBSOURCE ",$value->tagname, $value->as_text,"\n");
    #}
foreach my $key (keys %tagname_type){
   # print("Key: $key -> ",$tagname_type{$key}, "\n");

}
    #my $dbsource = $richSeq->get_dbsource();

#hmm..yet another wa
    my $acc = $richSeq->annotation();
    foreach my $key ($acc -> get_all_annotation_keys() ){
	@values = $acc->get_Annotations($key);
	foreach my $value (@values){

	    print "Annotation ", $key, " stringified value ", $value->as_text, "\n";
	}
    }



}

# $seq = getSeqFromGI($gi);
# Return sequence from Entrez (GenPept) given a GI number
# Inputs:  $gi - string representing gi number (can have gi| or just the number)
# outputs: $seq - string sequence from Entrez
# Kristen Naegle, May 23, 2007
sub getSeqFromGI($){
    my $query = shift;
    my $gb = new Bio::DB::GenPept();
    my $richSeq = $gb->get_Seq_by_id($query);


    if($richSeq){
	my $seq = $richSeq->seq();
	my $alphabet = $richSeq->alphabet();
	my $desired = 'protein';
	if ($alphabet ne $desired){
	    print("ERROR: Sequence is $alphabet instead of $desired\n");
	    return "NA";
	}
	else {
	    print("$seq\n");
	    return $seq;
	}
    }
    else{
	return "DNE";
    }
    
}

# $species = returnSpecies($richSeq)
# returns the species (common_name) form a richSeq object
# Inputs: $richSeq - richseq object (see returnRichSeqFromAcc
# Outputs: $species - the common name of species in richseq object
# Kristen Naegle
# March 7, 2008
# Modified April 2, 2008 - now returns binomial name
sub returnSpecies($){
    my $richSeq = shift;

    my $species = $richSeq->species();
# #    my $spec = $species->species();
# #    my $organelle = $species->organelle();
# #    my $commonName = $species->common_name();
    if(!$species){
	print "ERROR: species does not exist for rich seq object\n";
	return 'NA';
    }
    else{
	my $bi = $species->binomial();

	return lc($bi);
    }

}

# $primAcc = returnPrimaryAccession($richSeq)
# returns the primary Accession from a richSeq object
# Inputs: $richSeq - richseq object (see returnRichSeqFromAcc
# Outputs: $primAcc - primary accession
# Kristen Naegle
# March 27, 2008
sub returnPrimaryAccession($){
    my ($richSeq) = @_;
    my $primAcc = $richSeq->accession();
    return $primAcc;


}

# $richSeq = returnRichSeqFromAcc($acc)
# returns the richSeq object of a GenPept accession
# Inputs: $acc - accession (e.g. gi number)
# Outputs: $richSeq - richseq objec
# Kristen Naegle
# March 7, 2008
sub returnRichSeqFromAcc($){
    my $query = shift;
    my $gb = new Bio::DB::GenPept();
    my $richSeq;
    eval{
	### try block
	$richSeq = $gb->get_Seq_by_id($query);
    };
    if($@){
	### catch block
	$richSeq = -1;
	
    }
    return $richSeq;

}

# ($swissprot, $gi) = convertGeneToSwissProt($gene, $species)
# returns the first swiss prot accession for a gene and it's gi. If no swissprot exists returns the last gi record found and -1 for swissprot. If there was an error in looking up the gene, returns -1 for both
# Inputs: $gene - gene name
#         $species - species name (two words, e.g. Homo sapiens) by which to limit search
# Outputs: $swissprot - swiss prot protein accession
#          $gi - corresponding gi accession number
# Kristen Naegle
# May 13, 2008
sub convertGeneToSwissProt($$){
    my ($gene, $species) = @_;;
    my $gb = new Bio::DB::GenBank();
    my $richSeq;
    my $query = Bio::DB::Query::GenBank->new
	(-query => "$species\[Organism] AND $gene",
	 -db => 'protein');
    my $swissprot = -1;
    my $gi = -1;

    eval{
	
	### try block
	my $seqio = $gb->get_Stream_by_query($query);
	my $giNum;
	while (my $seq = $seqio->next_seq){
	    my $acc = $seq->accession;
	    $giNum = $seq->primary_id;
	    my $type = returnAccType($acc); 
	    if($seq->alphabet eq 'protein' && $type eq 'swissprot'){
		$swissprot = $acc;
		$gi = "gi|".$giNum;
		last;
	   }
	}
	if($swissprot eq '-1'){
	    $gi = "gi|".$giNum; #it will get the LAST GI..is this ok?
	}
    };
    if($@){
	### catch block
	$swissprot = -1;
	$gi = -1;
	
    }
    return ($swissprot, $gi);

}


# $proteinAcc = convertGeneToGI($gene, $species)
# returns the first gi accession for a gene. Assumes you didn't find a swissprot using convertGeneToSwissprot
# Inputs: $gene - gene name
#         $species - species name (two words, e.g. Homo sapiens) by which to limit search
# Outputs: $proteinAcc - gi protein accession
# Kristen Naegle
# May 27, 2008
sub convertGeneToGI($$){
    my ($gene, $species) = @_;;
    my $gb = new Bio::DB::GenBank();
    my $richSeq;
    my $query = Bio::DB::Query::GenBank->new
	(-query => "$species\[Organism] AND $gene",
	 -db => 'protein');
    my $protAcc = -1;
    eval{
	### try block
	my $seqio = $gb->get_Stream_by_query($query);
	
	while (my $seq = $seqio->next_seq){
 	    my $acc = $seq->accession;
 	    my $gi = $seq->primary_id();
 	    my $type = returnAccType($acc); 
 	    if($seq->alphabet eq 'protein' && $gi){
 		$protAcc = $gi;
 		last;
 	   }
	}
    };
    if($@){
	### catch block
	$protAcc = -2;
	
    }
    return $protAcc;

}

# ($errorCode, $seq, $species, $gene, $name, $locus) = getProteinFromAcc($acc);
# Return sequence from Entrez (GenPept) given an accession number
# Inputs:  $acc - accession
# outputs: $errorCode - returns 0 if retrieval worked, -1 if accession could not be found in Entrez and -2 if accession sequence was not protein, -3 if there was a gene name problem
#          $seq - string sequence from Entrez
#          $species - species
#          $gene - gene id
#          $name - name
# Kristen Naegle, May 23, 2007
sub getProteinFromAcc($){
    my ($query) = @_;

    my $richSeq = returnRichSeqFromAcc($query);
    
    my ($errorCode, $sequence, $species, $gene, $geneSynonyms, $name, $locus) = getProteinFromRichSeq($richSeq);
    return($errorCode, $sequence, $species, $gene, $geneSynonyms, $name, $locus);
    
}

# ($errorCode, $seq, $species, $gene, $geneSynonyms, $name, $primaryAcc, $newGI) = getProteinFromRichSeq($richSeq);
# Return sequence from Entrez (GenPept) given an accession number
# Inputs:  $richSeq - richSeq object
# outputs: $errorCode - returns 0 if retrieval worked, 1 if accession could not be found in Entrez and 2 if accession sequence was not protein, 3 if there was a gene name problem, 4 if there was a WARNING
#          $seq - string sequence from Entrez
#          $species - species
#          $gene - gene id
#          $name - name
#          $primaryAcc - primary accession
#          $newGI - this field is only returned if there was a change in the GI - to grandfather code in
# Kristen Naegle, May 23, 2007
# updated Feb. 2009 to include geneSynonyms
sub getProteinFromRichSeq($){
    my ($richSeq) = @_;

    my $errorCode = 0;
    my ($species, $sequence, $gene, $geneSynonyms, $WARN, $warnTxt, $name, $primaryAcc);
    
    if($richSeq){
	$species = returnSpecies($richSeq);
	$sequence = returnSequence($richSeq);
	if($sequence eq "NA"){
	    $errorCode = 2;
	    handleError('getProteinFromAcc', 'Accession does not lead to protein sequence', \@_);
	}
	if($species eq "NA"){
	    $errorCode = 2;
	    handleError('getProteinFromRichSeq', 'Accession does not have a species', \@_)
	}
#	my $name = $richSeq->accession();
	#$locus = $richSeq->display_name();
	$primaryAcc = returnPrimaryAccession($richSeq);
	$name = $richSeq->description();
	($gene, $geneSynonyms) = returnGeneGeneSynNames($richSeq);
	($WARN, $warnTxt) = checkForWarnings($richSeq);
	if($gene eq "ERROR"){
	  #  $errorCode = 3; don't throw error, allow null gene name
	    handleError('getProteinFromAcc', 'Could not get a gene', \@_);
	    $gene = 'NULL';
	}
    }
    else{
	$errorCode = 1;
	handleError('getProteinFromAcc', 'Could not get richSeq object, bad accession?', \@_);
	
    }
    my($CONV, $newGI)=0;
    if($WARN){
	$errorCode = 4;
	handleError('getProteinFromAcc', "WARNING Thrown: $warnTxt", \@_);
	# get major issue from warning text here
	($CONV, $newGI) = translateWarningFromRichSeq($warnTxt);
	print "DEBUG: Found Warning\t$warnTxt\n NEW GI: $newGI\n";
	return($errorCode, $sequence, $species, $gene, $geneSynonyms, $name, $primaryAcc, $newGI);

    }
    return($errorCode, $sequence, $species, $gene, $geneSynonyms, $name, $primaryAcc);
    
}

# ($errorCode, $seq, $species, $gene, $geneSynonyms, $name, $primaryAcc, $newGI) = getProteinFromRichSeq($richSeq);
# same as getProteinFromRichSeq but forces a return of newGI (0 if record was not replaced)
# Return sequence from Entrez (GenPept) given an accession number
# Inputs:  $richSeq - richSeq object
# outputs: $errorCode - returns 0 if retrieval worked, 1 if accession could not be found in Entrez and 2 if accession sequence was not protein, 3 if there was a gene name problem, 4 if there was a WARNING
#          $seq - string sequence from Entrez
#          $species - species
#          $gene - gene id
#          $name - name
#          $primaryAcc - primary accession
#          $newGI - 0 if no new record, otherwise returns gi
# Kristen Naegle, May 14, 2009
sub getProteinFromRichSeq_requireAll($){
    my ($richSeq) = @_;
    my ($errorCode, $sequence, $species, $gene, $geneSynonyms, $name, $primaryAcc, $newGI);
    my @values;
    @values = getProteinFromRichSeq($richSeq);
    if(scalar(@values) ==8){
	($errorCode, $sequence, $species, $gene, $geneSynonyms, $name, $primaryAcc, $newGI) = @values;
    } 
    else{
	($errorCode, $sequence, $species, $gene, $geneSynonyms, $name, $primaryAcc) = @values;
	$newGI = 0;
    }
    return ($errorCode, $sequence, $species, $gene, $geneSynonyms, $name, $primaryAcc, $newGI);
}

# ($CONV, $value) = translateWarningFromRichSeq($warnTxt);
# Intended to parse gi from warning and to eventually handle other tpes of warnings if they become an issue later on
# Inputs: $warnTxt - the [WARNING] text from a richSeq record
# Outputs: $CONV - a fixed conversion e.g. if a replacement it's REPLACE, else OTHER for now
#          $value - the value of that record, takes the first gi number
# Kristen Naegle
# May 15, 2009
sub translateWarningFromRichSeq($){
    my ($warnTxt) = @_;
    my($CONV, $value);
    $value =0;
    # could have multiple replacements - check date to take most recent. 
    if($warnTxt =~ m/replaced/){
	$CONV = "REPLACED";
	# parse gi number - assume first hit is most up to date
	$warnTxt =~ m/gi:([0-9]+)/g;
	$value = "gi|".$1;
	if($2) {
	    print "In translateWarningFromRichSeq found multiple gi updates, taking first: $value\n";
	}
    }
    else{
	$CONV = "OTHER";
    }
    return ($CONV, $value);
    
}

# $t = translateAccGetErrorCode($errorCode);
# Translates error codes in accession accesses to a readable string
# Inptus: $errorCode - errorCode (0-4)
# Outputs: $t - string describing code
# Kristen Naegle
# March 2008
sub translateAccGetErrorCode($){
    my $error = shift;
    my $string; 
    if($error == 0){
	$string = 'NO Error';
    }
    elsif($error == 1){
	$string = 'Bad Accession';

    }
    elsif($error == 2){
        $string = 'Not a protein entry (RNA or DNA)';
  
    }
    elsif($error == 3){

	$string = 'Could not retrieve Gene name';
    }
    elsif($error == 4){

	$string = 'Record had WARNING';
    }
    else {
	$string = 'Bad error code (should be 0, 1, 2, 3, or 4)';
    }
    return $string;

}

# $geneName = returnGeneName($richSeq)
# returns the gene name a richSeq object
# Inputs: $richSeq - richseq object (see returnRichSeqFromAcc
# Outputs: $geneName - the gene name of species in richseq object
# Kristen Naegle
# March 7, 2008
sub returnGeneName($){
    my $richSeq = shift;
    my @features = $richSeq->get_SeqFeatures;
    my $feature_tag; 
    my %seen_tag = ();
    my @tags = ();
    my @genes;
    foreach my $feature(@features){
	foreach my $tag (@tags = $feature->get_all_tags){
	    $seen_tag{$tag}++;
	    if($tag eq "gene"){
		#print("GENE: ", $feature->get_tag_values('gene'), "\n");
		push @genes, $feature->get_tag_values('gene');
	    }
	    

	}
    }
    my $old = $genes[0];
    my $finalGene;
    foreach my $gene (@genes){
	if ($old ne $gene){
	    $finalGene = "ERROR";
	    return $finalGene;
	}
	$old = $gene;
		
    }
    $finalGene = $old;
	
    #print("Genes: @genes\n");
    return $finalGene;
    
}

# ($geneName, $gene_synonyms) = returnGeneGeneSynName($richSeq)
# returns the gene name a richSeq object
# Inputs: $richSeq - richseq object (see returnRichSeqFromAcc
# Outputs: $geneName - the gene name of species in richseq object (the one with the most record votes becomes the pirmary gene name 
#          $geneSyn - an array of other gene names. 
# Kristen Naegle
# Feb. 9, 2009 - fixed Feb. 20, 2009 to parse gene_syn field for semicolon
sub returnGeneGeneSynNames($){
    my $richSeq = shift;
    my @features = $richSeq->get_SeqFeatures;
    my $feature_tag; 
    my %seen_tag = ();
    my @tags = ();
    my @genes;
    my @gene_syn;

    # First pass, gene is in primary tag
    foreach my $feature(@features){
	if($feature->primary_tag eq 'gene'){
	    for my $tag ($feature->get_all_tags){
		if($tag eq 'gene'){
		    push @genes, $feature->get_tag_values($tag);
		    
		}
		elsif($tag eq 'gene_synonym'){
		    push @gene_syn, $feature->get_tag_values($tag);
		}
	    }
	}

    }
    

    # second pass, gene is in tag
    foreach my $feature(@features){
	foreach my $tag (@tags = $feature->get_all_tags){
	    $seen_tag{$tag}++;
	    if($tag eq "gene"){
		#print("GENE: ", $feature->get_tag_values('gene'), "\n");
		push @genes, $feature->get_tag_values('gene');
	    }
	    elsif($tag eq "gene_synonym"){
		push @gene_syn, $feature->get_tag_values('gene_synonym');
	    }

	}
    }
    
    if(scalar(@genes) == 0){
	return ('', \@gene_syn);
    }

    # split @gene_syn for ; - agnostic with regards to where @gene_syn came from
    my @gene_synSplit;
    foreach my $val (@gene_syn){
	if($val =~ /;/){
	    my @genes = split(';', $val);
	    push @gene_synSplit, @genes;
	}
	else{
	    push @gene_synSplit, $val;
	}

    }

    # convert to hashes then 
    my $geneHash = arrToHash(\@genes);
    my $synHash = arrToHash(\@gene_synSplit);
    
#vote for gene name, remaining keys become synonyms
    my $i = 0;
    my @uniqueGenes = keys %$geneHash;
    my $votedGene = $uniqueGenes[$i];
    my $count = $geneHash->{$votedGene};
    for(my $i=1; $i<=$#uniqueGenes; $i++){
	my $nextGene = $uniqueGenes[$i];
	if($geneHash->{$nextGene} > $count){
	    $votedGene = $nextGene;
	    $count = $geneHash->{$nextGene};
	}

    }
    # at the end $votedGene is final gene, push all others onto synonym hash and return keys of that for gene_syn! Remove $gene from both hashes. 
    delete($geneHash->{$votedGene});
    delete($synHash->{$votedGene});
    foreach my $gene (keys %$geneHash){
	if(not defined($synHash->{$gene})){
	    $synHash->{$gene} = 0;
	}
	$synHash->{$gene} += 1;
    }
    
    # Make sure there is no space at the beginning - don't know why this is a bug but it is
    $votedGene =~ s/^ //;
    my @voted_syn; # = keys %$synHash;
    foreach my $k (keys %$synHash){
	$k =~ s/^ //; #sometimes get spaces at the beginning of synonym names
       # add the substituion here for ' in gene syn names
#	$k =~ s/\'/\\\'/g;
	push @voted_syn, $k;
    }
    return ($votedGene, \@voted_syn);
    
}

# ($WARN, $warningTxt) = checkForWarnings($richSeq);
# Given a richSeq object, check the comments field in annotations for WARNING, return the warning string and a BOOL value
# Inputs: $richSeq - richSeq object
# Outputs: $WARN - boolean value, 1 if there was a warning, 0 otherwise
# Kristen Naegle 
# Feb 9, 2009
sub checkForWarnings($){
    my ($richSeq) = @_;
    my $WARN = 0;
    my $warnStr;
    my $anno_collection = $richSeq->annotation;
    my @annotations = $anno_collection->get_Annotations('comment');
    foreach my $ann (@annotations){
	if ($ann =~ m/WARNING/){
	    $WARN = 1;
	    $warnStr = $ann;
	}
    }
    return ($WARN, $warnStr);
}

# $sequence = returnSequence($richSeq)
# returns the PROTEIN sequence form a richSeq object
# Inputs: $richSeq - richseq object (see returnRichSeqFromAcc
# Outputs: $sequence - the sequence, if the alphabet is not protein, returns "NA"
# Kristen Naegle
# March 7, 2008
sub returnSequence($){
    my $richSeq = shift;
    
    my $seq = $richSeq->seq();
    my $alphabet = $richSeq->alphabet();
    my $desired = 'protein';
    if ($alphabet ne $desired){
	#print("ERROR: Sequence is $alphabet instead of $desired\n");
	return "NA";
    }
    else {
	#print("$seq\n");
	return $seq;
    }
    

}

# (\@positions,\@codes, $pepCap) = returnPhosphoPos($peptide)
#  Returns the arrays of positions and codes of the phosphorylated residues in a peptide sequence as indicated by pS, pT, and pY
# Also returns a capitalized form without the p in front of the peptide
# inputs: $peptide (from MS with pY, pS, or pT representing phosphorylated versions;
sub returnPhosphoPos($){
    my $pep = shift;
    #first replace pT, pS, pY with t,s,y 
    $pep =~ s/pS/s/g;
    $pep =~ s/pY/y/g;
    $pep =~ s/pT/t/g;
    my $pepCap = uc($pep);
    #print("$pepCap\n");
    my @pos;
    my @codes;
    my @pep = split(//, $pep);
    #now find instances of LC s, y, and t
    my $count = 0;
    foreach my $aa (@pep){
	#if the $aa is a lowercase then it's phosphorylated
	if($aa eq lc($aa)){
	    push(@pos, $count);
	    push(@codes, uc($aa));
	}
	$count+=1;
    }
    if(@pos){
	#print "Position @pos\n";
	return (\@pos, \@codes, $pepCap);
    }

    else{
#	print("Could not find phosphorylation in peptide: $pep\n");
	return("NA", "NP", $pepCap);
    }
}

# (\@pos, \@codes) = returnAlignedandCodes($gi, $pep, $numAA)
# Given a gi number return the position and codes and the aligned sequences blown out to +/-numAA
#If peptide is no phosphorylated returns "NA"
# inputs: $gi - gi accession number
#         $pep - MS peptide fragment (with pS, pY or pT marking alignment)
#         $numAA - number of aa on each side to return
# outputs: \@seqPosCodes - reference to an array of codes and their position in the full protein (e.g. Y110)
#          \@aligned - corresponding aligned sequences where pS is s pY is y and pT is t
# If peptide could not be found in sequence returns "NA"
# If peptide not found in sequence b/c record is RNA or DNA then returns "NA" in code and "RNA" in alignment too, if sequence did not exist in protein returns "NA" and "DNE"
# Kristen Naegle
# May 23, 2007
sub returnAlignedandCodes($$$){
    my $gi = shift;
    my $pep = shift;
    my $numAA = shift;
    my $seq; 
    $seq = getSeqFromGI($gi);
    #should be able to get length here
    
    #need to strip spaces
    $pep =~ s/ //g;

    #if sequence is RNA or DNA alphabet
    if($seq eq "NA"){
	return("NA", "RNA");
    }
    elsif($seq eq "DNE"){
	return("NA", "DNE");
    }
    else{

# #now find the position of the phosphorylations 
# #given a peptide in MS format, return capitalized peptide and array of positions
    my $codes;
    my $pos;
    my $pepCap;
    ($pos, $codes, $pepCap) = returnPhosphoPos($pep);
# #   print("$pep\n");
# #   print("$pepCap\n");

# #handle case where there is no phosphorylation position 
    if($pos eq "NA"){
	return($pos, $codes); #these contain error messages
    }
    else{

#now get the position of that residue in the full sequence
    my $index = index($seq, $pepCap);
    print "$pepCap";
    print "$index\n";

    if($index == -1){ #peptide was not found in sequence
#	if($VERBOSE){
	    print "ERROR: Cannot find position of phosphorylation of $pep in $seq\n";
	#}
	return("NA", "PNS", $pepCap);
	#print $index;
    }
    else{
	print("Index of $pepCap is $index\n");
	my $i;
	my @pos = @$pos;
	my @aligned;
	my @seqPos; 
	my @seqPosCode;
	for($i=0; $i<=$#pos; $i++){
	    #print("$pos->[$i] $codes->[$i]\n");
	    push(@seqPos, $pos->[$i]+$index+1); #have to add one for alignedSequences
	    push(@aligned, returnAlignedSequence($numAA,$seq, $codes->[$i], $seqPos[$i]));
	    push(@seqPosCode, $codes->[$i] . $seqPos[$i]);
	    
	} #end for
	print("Aligned sequences: @aligned\n");
#in the end return the T163 and the aligned sequence
	return(\@seqPosCode, \@aligned);
    } #end else
}
}

}


# appendAlignedSeqToFileGivenGI($inputFile, $outputFile, $giCol, $pepCol, $numAA)
# Take an inputfile with at least a gi accession number and peptide detected from MS and append at the beginning a peptide identifier and at the end the aligned sequence
# Modified 8/8/07 to handle header line in file
# inputs: $inputFile - txt file of data (pS, pY, and pT)
#         $outputFile - will be overwritten with original data and appended sequence/code and aligned peptide
#         $giCol - 1's based column number of gi acc number
#         $pepCol - 1's based column number of peptide fragment
#         $numAA - number of amino acids to align on each side 
# updated 10/25/07 - to shift for fileTools zero-based, and to put multiple sequence alignments in a single file.   (copied old below)
sub appendAlignedSeqToFileGivenGI($$$$$){
    my $inputFile = shift;
    my $outputFile = shift;
    my $giCol = shift;
    my $pepCol = shift;
    my $numAA = shift;

    open(IN, $inputFile) || die "Can't open input file: $inputFile\n";

    if(!-e $outputFile) { system("touch $outputFile");}
    open(OUT, ">$outputFile") || die "Can't open output file for writing: $outputFile\n";

    #write the header line to the output
    my $header = returnHeader($inputFile);
    chomp $header;
    my $headerNew = "pep:uniq\t".$header."\tsite:aligned\tpep:aligned\n";
    print OUT "$headerNew";
    my $line = <IN>; #skip header row
    my $peptideCounter = 0;
    while(defined($line=<IN>)){
	chomp $line;
	my @lineArr = split("\t", $line);
	my $gi = $lineArr[$giCol];
	my $pep = $lineArr[$pepCol];
	my $codePos;
	my $aligned; 
	($codePos, $aligned) = returnAlignedandCodes($gi, $pep, 7);

#If y code exists then double check it
	my $i;

	if ($aligned ne "RNA" & $codePos ne "NA"){
	    #my @a = @$aligned;
	    print(OUT "$peptideCounter\t");
	    print OUT $line;
	    my $codeT = $codePos->[0];
	    my $alignedT = $aligned->[0];
	    
	    if(scalar(@$aligned) > 1){
		for($i = 1; $i< scalar(@$aligned); $i++){
		    $codeT .= ";$codePos->[$i]";
		    $alignedT .= ";$aligned->[$i]";

		}
	    }
	    print OUT "\t$codeT\t$alignedT\n";
	    
# # 	    for($i=0; $i<=$#a; $i++){
# # 		#strip \n
# # #
# # 	#	print(OUT $line);
		
# # 		$appendate .= "\t$codePos"
# # 	    }
	}
	else{
	    print(OUT "$peptideCounter\t$line");
# #if($aligned eq "RNA"){
# #	print(OUT "\tRNA\n");
# #
# #    }
# #   else{
# #print(OUT "\tNA\n");
# #}
	    print(OUT "\t$codePos\t$aligned\n"); #Should say RNA or DNE if there was an error
	}

	chomp $line;
	$peptideCounter+=1;

    }

    close(IN);
    close(OUT);
}

# $gene = getGeneNameGivenGI($gi)
# returns the gene name from the feature fields of a genpept record 
# Checks to make sure all gene names of features agree
# inputs: $gi - accession number
# outputs: $gene - returns "ERROR" if the array of features do not agree
# July 3, 2007
# Kristen Naegle
sub getGeneNameGivenGI($){
    my $gi = shift;
    my $gb = new Bio::DB::GenPept();
    my $richSeq = $gb->get_Seq_by_id($gi);
    my @features = $richSeq->get_SeqFeatures;
    my $feature_tag; 
    my %seen_tag = ();
    my @tags = ();
    my @genes;
    foreach my $feature(@features){
	foreach my $tag (@tags = $feature->get_all_tags){
	    $seen_tag{$tag}++;
	    if($tag eq "gene"){
		print("GENE: ", $feature->get_tag_values('gene'), "\n");
		push @genes, $feature->get_tag_values('gene');
	    }
	    

	}
    }
    my $old = $genes[0];
    my $finalGene;
    foreach my $gene (@genes){
	if ($old ne $gene){
	    print "ERROR: Gene names in record $gi do not agree, $old and $gene\n";
	    $finalGene = "ERROR";
	    return $finalGene;
	}
	$old = $gene;
		
    }
    $finalGene = $old;
	
    print("Genes: @genes\n");
    return $finalGene;
    
}

# appendGeneNameToFileGiveGI($inputFile, $outputFile, $giCol)
# Append Gene name to file given accession number
# inputs: $inputFile - tab dileneated file 
#         $outputFile - output file that will have a last column appended gene name
# UPDATED 8/10/07 to handle header
# July 3, 2007, Updated 10/26/07 to only look up genes once 
# Kristen Naegle
sub appendGeneNameToFileGivenGI($$){
    my $inputFile = shift;
    my $outputFile = shift;
# #   my $giCol = shift;

# # First get unique Accession numbers in the file
    my $accCol = returnAccCol($inputFile);
    my $giHash = returnUniqColumnHash($inputFile, $accCol);
    foreach my $item (keys %$giHash){
	print "GI: $item\n";
    }


# Next get gene names for each 
    my %giGeneHash; 
    foreach my $gi (keys %$giHash){
	if(not $giGeneHash{$gi}){
	    $giGeneHash{$gi} = getGeneNameGivenGI($gi);
	    print "$gi -> $giGeneHash{$gi}\n";
	}
    }

# Finally write this to the output 
    my $giCol = returnAccCol($inputFile);
    open(FH_3, $inputFile) || die "Can't open input file: $inputFile\n";

    if(!-e $outputFile) { system("touch $outputFile");}
    open(OUT, ">$outputFile") || die "Can't open output file for writing: $outputFile\n";
    #take care of header line
    my $line = <FH_3>; 
    my $header = returnHeader($inputFile);
    chomp $header;
    $header = $header."\tacc:gene";
    print OUT "$header\n";
    
    my $peptideCounter = 0;
    while(defined($line=<FH_3>)){
	chomp $line;
	my @lineArr = split("\t", $line);

	my $gi = returnField($line, $giCol);
	my $gene = $giGeneHash{$gi};
	if(defined($gene)){
	    print(OUT "$line\t$gene\n");
	}
	else{ print(OUT "~~~\n"); }

    }
    close(OUT);
    close(FH_3);
    

}

# \%richSeqHash = returnRichSeqHash($accHash)
# Given a hash of accessions, returns a hash of GenPept rich seq objects
# Inputs: $accHash - ref. to hash of accessions
# Outputs: \%richSeqHash - ref to hash of rich seq objects with key accession
# Kristen Naegle
# ? September, 2007 ?
sub returnRichSeqHash($){
    my $accHash = shift;
    my %richSeqHash; 
    my $gb = new Bio::DB::GenPept();
    $gb->verbose(1);
    foreach my $acc (keys %$accHash){

	$acc =~ s/ //g;
# #$acc =~ s/gi\|//g;
# #print "Hello\n";
	print "Checking $acc\n";
	#print "World\n";
	$richSeqHash{$acc} = $gb->get_Seq_by_id($acc);
	#print "BAM\n";
    }

    return \%richSeqHash;

}

# $gene = returnGeneFromRichSeq($richSeq)
# returns the first gene name of features in richSeq object.
# Kristen Naegle
# ? September 2007 ? obsolete
sub returnGeneFromRichSeq($){

    my $richSeqFeature = shift;
        #print Gene name
    my @features = $richSeqFeature->get_SeqFeatures;
    my $feature_tag; 
    my %seen_tag = ();
    my @tags = ();
    my @genes;

    foreach my $feature(@features){
	foreach my $tag (@tags = $feature->get_all_tags){
	    $seen_tag{$tag}++;
	    if($tag eq "gene"){
		#print("GENE: ", $feature->get_tag_values('gene'), "\n");
		push @genes, $feature->get_tag_values('gene');
		last;
	    }
	    

	}
    }
    my $gene = $genes[0];

    return $gene;
}


# appendGeneNameToFileFromHash($inputFile, $outputFile, $richSeqHash)
# Append Gene name to file given a richSeq Hash from unique acc numbers
# inputs: $inputFile - tab dileneated file 
#         $outputFile - output file that will have a last column appended gene name
#         $richSeqHash - hash of richSeq (DB gets)
# November 5, 2007
# Kristen Naegle
sub appendGeneNameToFileFromHash($$$){
    my $inputFile = shift;
    my $outputFile = shift;
    my $richSeqHash = shift;

    my $giCol = returnAccCol($inputFile);
    open(FH_3, $inputFile) || die "Can't open input file: $inputFile\n";

    if(!-e $outputFile) { system("touch $outputFile");}
    open(OUT, ">$outputFile") || die "Can't open output file for writing: $outputFile\n";
    #take care of header line
    my $line = <FH_3>; 
    my $header = returnHeader($inputFile);
    chomp $header;
    $header = $header."\tacc:gene";
    print OUT "$header\n";
    
    while(defined($line=<FH_3>)){
	chomp $line;
	my $gi = returnField($line, $giCol);
	#my $gene = $giGeneHash{$gi};
	my $gene;
	if($richSeqHash->{$gi}){
	    $gene = returnGeneFromRichSeq($richSeqHash->{$gi});
	}
	if(defined($gene)){
	    print(OUT "$line\t$gene\n");
	}
	else{ print(OUT "~~~\n"); }

    }
    close(OUT);
    close(FH_3);
}


# appendAlignedSeqToFileFromHash($inputFile, $outputFile, $richSeqHash, $numAA)
# Take an inputfile with at least a gi accession number and peptide detected from MS and append at the beginning a peptide identifier and at the end the aligned sequence
# inputs: $inputFile - txt file of data (pS, pY, and pT)
#         $outputFile - will be overwritten with original data and appended sequence/code and aligned peptide
#         $richSeqHash - hash of rich seq features
#         $numAA - number of amino acids to align on each side 
# November 5, 2007
# Kristen Naegle
sub appendAlignedSeqToFileFromHash($$$$){
    my $inputFile = shift;
    my $outputFile = shift;
    my $richSeqHash = shift;
    my $numAA = shift;

    open(IN, $inputFile) || die "Can't open input file: $inputFile\n";

    if(!-e $outputFile) { system("touch $outputFile");}
    open(OUT, ">$outputFile") || die "Can't open output file for writing: $outputFile\n";

    #write the header line to the output
    my $header = returnHeader($inputFile);
    chomp $header;
    my $headerNew = "pep:uniq\t".$header."\tsite:aligned\tpep:aligned\n";
    print OUT "$headerNew";
    my $line = <IN>; #skip header row
    my $peptideCounter = 0;
    my $giCol = returnAccCol($inputFile);
    my $pepCol = returnColumnNumber($inputFile, "pep:tryps");
    while(defined($line=<IN>)){
	chomp $line;
	my @lineArr = split("\t", $line);
	my $gi = returnField($line, $giCol);
	my $pep = returnField($line, $pepCol);
	my $codePos;
	my $aligned; 
	($codePos, $aligned) = returnAlignedandCodesFromHash($gi, $pep, $richSeqHash, $numAA);
#If y code exists then double check it
	my $i;

	if ($aligned ne "RNA" & $codePos ne "NA"){
	    #my @a = @$aligned;
	    print(OUT "$peptideCounter\t");
	    print OUT $line;
	    my $codeT = $codePos->[0];
	    my $alignedT = $aligned->[0];
	    
	    if(scalar(@$aligned) > 1){
		for($i = 1; $i< scalar(@$aligned); $i++){
		    $codeT .= ";$codePos->[$i]";
		    $alignedT .= ";$aligned->[$i]";
		    
		}
	    }
	    print OUT "\t$codeT\t$alignedT\n";
	    
	}
	else{
	    print(OUT "$peptideCounter\t$line");
	    print(OUT "\t$codePos\t$aligned\n"); #Should say RNA or DNE if there was an error
	}

	chomp $line;
	$peptideCounter+=1;

    }

    close(IN);
    close(OUT);
}

# (\@pos, \@codes) = returnAlignedandCodes($gi, $pep, $numAA)
# Given a gi number return the position and codes and the aligned sequences blown out to +/-numAA
#If peptide is no phosphorylated returns "NA"
# inputs: $gi - gi accession number
#         $pep - MS peptide fragment (with pS, pY or pT marking alignment)
#         $numAA - number of aa on each side to return
# outputs: \@seqPosCodes - reference to an array of codes and their position in the full protein (e.g. Y110)
#          \@aligned - corresponding aligned sequences where pS is s pY is y and pT is t
# If peptide could not be found in sequence returns "NA"
# If peptide not found in sequence b/c record is RNA or DNA then returns "NA" in code and "RNA" in alignment too, if sequence did not exist in protein returns "NA" and "DNE"
# Kristen Naegle
# May 23, 2007
sub returnAlignedandCodesFromHash($$$$){
    my $gi = shift;
    my $pep = shift;
    my $richSeqHash = shift;
    my $numAA = shift;
    my $seq; 
  
    if($richSeqHash->{$gi}){
	$seq = getSeqFromHash($richSeqHash->{$gi});
#	print "Sequence: $seq\n";
    }
    else{
	$seq = "DNE";
	print "ERROR: $gi does not exist as key in rich sequence hash in pepalign attempt\n";
    }
    #should be able to get length here
    
    #need to strip spaces
    $pep =~ s/ //g;

    #if sequence is RNA or DNA alphabet
    if($seq eq "NA"){
	return("NA", "RNA");
    }
    elsif($seq eq "DNE"){
	return("NA", "DNE");
    }
    else{

# #now find the position of the phosphorylations 
# #given a peptide in MS format, return capitalized peptide and array of positions
    my $codes;
    my $pos;
    my $pepCap;
    ($pos, $codes, $pepCap) = returnPhosphoPos($pep);

#handle case where there is no phosphorylation position 
    if($pos eq "NA"){
	return($pos, $codes); #these contain error messages
    }
    else{

#now get the position of that residue in the full sequence
	my $index = index($seq, $pepCap);
# #	print "$pepCap";
# #	print "$index\n";

	if($index == -1){ #peptide was not found in sequence
	    return("NA", "PNS", $pepCap);
	    #print $index;
	}
	else{
	    #print("Index of $pepCap is $index\n");
	    my $i;
	    my @pos = @$pos;
	    my @aligned;
	    my @seqPos; 
	    my @seqPosCode;
	    for($i=0; $i<=$#pos; $i++){
		#print("$pos->[$i] $codes->[$i]\n");
		push(@seqPos, $pos->[$i]+$index+1); #have to add one for alignedSequences
		push(@aligned, returnAlignedSequence($numAA,$seq, $codes->[$i], $seqPos[$i]));
		push(@seqPosCode, $codes->[$i] . $seqPos[$i]);
		
	    } #end for
# #	    print("Aligned sequences: @aligned\n");
# #in the end return the T163 and the aligned sequence
	    return(\@seqPosCode, \@aligned);
	} #end else
    }
}
    
}

# $seq = getSeqFromHash($richSeq);
# Return sequence from a richSeq object
# Inputs:  $richSeq - richSeq object from entrez query 
# outputs: $seq - string sequence from Entrez
# Kristen Naegle, May 23, 2007
sub getSeqFromHash($){
    my $richSeq = shift;

    #if($richSeq){
    my $seq = $richSeq->seq();
    my $alphabet = $richSeq->alphabet();
    my $desired = 'protein';
    if ($alphabet ne $desired){
	print("ERROR: Sequence is $alphabet instead of $desired\n");
	return "NA";
    }
    else {
	#print("$seq\n");
	return $seq;
    }
    
}

# $type = returnAccType($acc)
# Returns the accession type enum: 'gi', 'entrez_protein', 'swissprot', and 'undefined'
# Inputs: $acc - string: accession
# Outputs: $type - type of known databases 
# Kristen Naegle
# March 6, 2008
sub returnAccType($){
    my $acc = shift;
    my $type;
    if($acc =~ /gi/){
	$type = 'gi';

    }
    elsif($acc =~ /^[NXZ]P_\d+/){
	$type = 'entrez_protein';


    }
# #   elsif($acc =~ /^([A\-N]|[R\-Z])\d([A\-Z])([A\-Z]|\d)([A\-Z]|\d)\d/){
# #   elsif($acc =~ /^([O|P|Q])\d([A\-Z]|\d)([A\-Z]|\d)([A\-Z]|\d)\d$/){
    elsif($acc =~ /^[O|P|Q]\d...\d$/){
	$type = 'swissprot';
    }
    elsif($acc =~ /^[A-N|R-Z]\d[A-Z]..\d/){
	$type = 'swissprot';
    }
    elsif($acc =~ /^[A-Z]{3}\d{5}$/){
	$type = 'genbank';
    }
    else{
	$type = "undefined";
	print "DEBUG: Type undefined for acc: $acc\n";

    }
    return $type;    
}

# $newAcc = removeAccIsoforms($acc)
# Given an accession remove the isoform tails -number or .number (used in refseq and entrez)
# Inputs: $acc - original accession
# Outputs: $newAcc - original accession stripped of isoforms
# Kristen Naegle
# June 28, 2009
sub removeAccIsoforms($){
    my ($acc) = @_;
    
    $acc =~ s/-\d+//;
    $acc =~ s/\.\d+//;
    return $acc;

}

# ($sites, $codes) = parseSiteCodesFromRP($rp)
# Given an RP node line from a swissprot richseq object, parse out the sites of phosphorylation
# Inputs: $rp - the rp text line of: PHOSPHORYLATION [LARGE SCALE ANALYSIS] AT"
# Outputs: $sites - ref. to array of sites
#          $codes - ref. to corresponding array of codes (S, Y or T)
# Kristen Naegle
# Feb. 3, 2009
sub parseSiteCodesFromRP($){
    my ($rp) = @_;
    my @sites;
    my @codes;

    my %convHash;
    $convHash{'SER'} = 'S';
    $convHash{'TYR'} = 'Y';
    $convHash{'THR'} = 'T';

    if($rp =~ m/PHOSPHORYLATION/ and $rp =~ m/LARGE SCALE ANALYSIS/){
	# match TYR-Number, SER-number, THR- number
	
	while($rp =~ m/(SER|TYR|THR)-(\d+)/){
	    push @sites, $2;
	    push @codes, $convHash{$1};
	    $rp =~ s/$2//;
	    #print "test: $1\n";
	}
	
    }
    return (\@sites, \@codes);

}


1;

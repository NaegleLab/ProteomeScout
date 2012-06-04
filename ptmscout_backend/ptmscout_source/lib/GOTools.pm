use strict;
use warnings;
use GO::AnnotationProvider::AnnotationParser;
use DBTools::queryTools;
use DBTools::insertFunctions;
use DBTools::dbEntry;
use errorHandling;
use GO::Parser;
use commonTools;
use globalVars;
#use mirTools;

# $standardName = returnStandardNameForAlias($annotationFile, $alias)
# Given an annotation file, find the alias and return StandardName
# When more than one standardName exists for the correct taxonomy, chooses the one with more entries.
# Inputs: $annotationFile - annotation file from geneOntology.org
#         $alias - the name of the gene for which you believe appears in alias
# Outputs: $standardName - standard name for gene
# Kristen Naegle
# April 15, 2008
sub returnStandardNameForAlias($$$){
    my ($annotationFile, $alias, $taxonomy) = @_;
    my $standardName=-1;
    my $id = -1;
    my @lines = `grep $alias $annotationFile`;
    my %idHash;
    my %nameHash;
    foreach my $line (@lines){
	my($id, $name,$taxon) = parseGOALine($line);
	if(checkAliases($line, $alias)){
#	print "Found $id\t$name\t$taxon\n";
	    if($taxon =~ /$taxonomy/){
		if(not defined $idHash{$id}){
		    my %mini;
		    $mini{'count'} = 1;
		    $mini{'name'} = $name;
		    $idHash{$id} = \%mini;
#		$idHash{$id}{'name'} = $name;
		}
		else{
		    $idHash{$id}->{'count'} += 1;
		    #	$nameHash{$name} +=1;
		}
	    }
	} # end check that alias was correctly in line
    }
    my @keys = keys %idHash;
    my $i = 0;
    my $k = $keys[$i];
    if(scalar(@keys) > 1){

	$i+=1;
	while($i <= $#keys){
	    if($idHash{$keys[$i]}->{'count'} > $idHash{$k}->{'count'}){
		$k = $keys[$i];
	    }
	    $i+=1;
	}

    }

    if($k){
	$standardName = $idHash{$k}->{'name'};
	$id = $k;
    }
    return ($id, $standardName);


}

# $check = checkAliases($aliasesFromLine, $alias)
# Given a line, make sure that the alias matches exactly (i.e. word boundary match, more than just a grep
# Inputs: $aliasesFromLine - this could be a whole GOA line, or just the alias portion
#         $alias - the alias you are searching for
# Outputs: $check - 1 if found, 0 otherwise
# Kristen Naegle
# April 15, 2008
sub checkAliases($$){
    my ($aliasesFromLine, $alias) = @_;
    my $check = 0;
    # value should fit full length of word in line, either based on space character or comma, or colon, i.e. should not be followed by another alphanumeric character
   # print "Finding $alias in $aliasesFromLine\n";
    if($aliasesFromLine =~ /(^|\s)$alias(:|,|;|\s)/){
	$check = 1;
    }
    return $check;
}


# ($dbId, $standardName, $aliases, $taxonomy) = parseGOALine($line);
# Given a GOA line return portions of it
# Inputs: $line
# Outputs: $dbId- the first GO useable id 
#          $standardName - GO useable standard name
#          $taxonomy - taxonomy, example: Human is taxon:9606
# Kristen Naegle
# April 15, 2008
sub parseGOALine($){

    my ($line) = @_;
    chomp $line;
    my @line = split("\t", $line);
    my ($db, $dbId, $standardName, $GO, $REF, $foo, $source, $type, $aliases,$foo2,$foo3,$foo4, $taxonomy, $bar, $source2) = @line;
    return($dbId, $standardName, $taxonomy);
}

# ($error, $associations) = returnGOAssociations($annotationParser, $names, $databaseIds, $aspect);
# Returns array of GO annotations for the databaseIds and the names given
# Inputs: $annotationParser - annatation parser object from GO::AnnotationProvider
# REMOVED - new version of GO now looks up by alias        $annotationFile - file parser was created from, need this because parser lookup does not seem to be working with aliases..so built in house alias look up
# CHANGED to pass in arrays of names and databaseIds instead      $gene - the gene of interst, swissprot/uniprot preferabl
#         $names - array of names that have annotations for teh gene of interest#         $databaseIds - array of databaseIds for gene of interest
#         $aspect - F, P, or C for function, process or cellular location
# REMOVED - now assumes annotationParser is for correct species         $taxonomy - taxonomy of interest (also for alias look up) e.g. human is 9606
# Outputs: $error - error is raised if annotation for gene could not be found
#          $associations - the associations GO
# Kristen Naegle
# April 15, 2008
# Mod May 21st, 2009 - now gets all go ids from all possible names no longer use aliases.
sub returnGOAssociations($$$$){
    my($annotationParser, $names, $databaseIds, $aspect) = @_;
    my $error = 0;
    
    my @associations; 
    
    foreach my $name (@$names){
	my $associations = $annotationParser->goIdsByName(name=>$name, aspect=>$aspect);
	if($associations){
	    push @associations, @$associations;
	}
	
    }
    foreach my $dbId (@$databaseIds){
	my $associations = $annotationParser->goIdsByDatabaseId(databaseId=>$dbId, aspect=>$aspect);
	if($associations){
	    push @associations, @$associations;
	}
    }
    # create a unique array of go terms
    
    my $unique = returnUniqueArray(\@associations);
    
    return ($error, $unique);
    
    
}

# handleGOTermsForProteinId($dbh, $annotationParser, $ontology, $proteinId)
# Retrieve associations and insert into GO tables for a proteinId
# Inputs: $dbh - databasehandle
#         $annotationParser - correct annotation parse of a species specific annotation file
#         $ontology - ontology parser
#         $proteinId - protein id of interest for insertion
# Outputs: \%GOHash - ref. to GO hash, keys are aspect terms and value is an array of the GO terms that are loaded 
# Kristen Naegle
# May 22, 2009
sub handleGOTermsForProteinId($$$$$$){
    my ($dbh, $annotationParser, $ontology, $proteinId, $OBO_VER, $ANN_VER) = @_;
    print "DEBUG: Ontology Version: $OBO_VER\tAnnotation Version: $ANN_VER\n"; 

    my @aspects = ('C', 'F', 'P');
    my %GOTerms;
    
    # first get names, then get GO for each aspect and insert
    my $AMB = 0; #eventually want to pass this in as argument?
    my($databaseIds, $names) = returnNamesForProteinId($dbh, $proteinId, $annotationParser, $AMB);
    
    if(!$databaseIds and !$names){
	handleError('handleGOTermsForProteinId', "Protein $proteinId has no annotated terms based on protein accessions", \@_);
	return \%GOTerms;
    }
    
    foreach my $aspect (@aspects){
	my ($error, $associations) = returnGOAssociations($annotationParser, $names, $databaseIds, $aspect);
	my ($hashRef) = insertGOTerms($dbh, $proteinId, $associations, $ontology, $aspect, $OBO_VER, $ANN_VER);
	$GOTerms{$aspect} = [];
	foreach my $GO (@$associations){
	    push @{$GOTerms{$aspect}}, $ontology->get_term($GO);
	}
    }
    return \%GOTerms;
}

# handleGOTermsForProteinIds($dbh, $ontologyFile, $proteinIds);
# Inserts F, C, and P associations for a list of protein ids - reports error if no accession can be found for a protienId for whcich to look up GO
# Inputs: $dbh = database handle
#         $proteinIds - array of protein.id 
# Outputs: $GOProteinHash - ref to hash with top key protein Id points to GOHash. GOHash has keys equal to aspect and values equal to arrays of GO hash terms. 
# Kristen Naegle
# April 16, 2008
# May 25, 2009 - updated to find taxonomy with each protein
# June 23, 2009 - updated to look up ontology file inside here 
sub handleGOTermsForProteinIds($$){
    my ($dbh, $proteinIds) = @_;
    
    my %GOProteinHash;

    my ($rJunk, $rJunk2, $ontologyFile) = returnGlobalGOFileInfo();
    my $parser = new GO::Parser({handler=>'obj'}); 
    $parser->parse($ontologyFile);
    my $ontology = $parser->handler->graph;

    my ($ontologyVer, $fileJ) = returnLinkRelease($ontologyFile);

    
    # for a list of protein Ids, get the species so we can get the annotationparserhash
    my $proteinStr = makeSelectINString($proteinIds, 0);
    my $sth = $dbh->prepare("SELECT species from protein where id IN $proteinStr group by species");
    $sth->execute();
    my $speciesArr = returnArrayOfResultsOnCol($sth, 0);

    my ($annotationHash, $verHash) = returnAnnotationParserHash($speciesArr);
    my @badSpecies;

    foreach my $proteinId (@$proteinIds){
	my $proteinHash = returnProteinDescForProteinId($dbh, $proteinId);
	my $annotationParser = $annotationHash->{$proteinHash->{'species'}};
	my $ANN_VER = $verHash->{$proteinHash->{'species'}};
	if(!$annotationParser){
	    push @badSpecies, $proteinHash->{'species'};
	    next;
	}
	my $GOHash = handleGOTermsForProteinId($dbh, $annotationParser, $ontology, $proteinId, $ontologyVer, $ANN_VER);
	print "DEBUG: Success on parser for protein Id: $proteinId\n";
	#printHash($GOHash);
	$GOProteinHash{$proteinId} = $GOHash;
    }
    
    #write errors to log file
    if(@badSpecies){
	my $uniqueBad = returnUniqueArray(\@badSpecies);
	handleError('handleGOTermsForProteinIds', "Could not get GO terms for these species @$uniqueBad\n", \@_);
	print "DEBUG: These species have no GO terms @$uniqueBad\n";
    }
    return (\%GOProteinHash);
}


# handleGoTermsForProteinIds_viaExpId($dbh, $expId);
# Only works on human taxonomy right now 
# DOES NOT COMMIT
# Inputs: $dbh - database handle
#         $expId - id from experiment table
# Kristen Naegle
# May 25th, 2009 - totally updated
sub handleGOTermsForProteinIds_viaExpId($$){
    my($dbh, $expId) = @_;
    #my $SLIM = 1;
   # my $taxonomy = '9606';
    my $proteinIds = returnProteinIdsForExpId($dbh, $expId);
    my $GOHash = handleGOTermsForProteinIds($dbh, $proteinIds);
    
    #$dbh->commit();
    return $GOHash;
    
}



# $taxonomy = convertBinomialSpeciesToTaxonomy($species)
#given the two word sceintific species name, return the taxonomy
# Inputs: $species - two name species, fom protein table
#         $taxonHash - ref. to taxonomy hash see returnSpeciesTaxonHash();
# Outputs: $taxonomy - numeric taxonomy term for GO
# Kristen Naegle
# May 19, 2009
sub convertBinomialSpeciesToTaxonomy($$){
    my ($species, $taxonHash) = @_;
    if($taxonHash->{$species}){
	return $taxonHash->{$species};
       
    }
    else{
	return -1;
    }

}

# $taxonHash = returnSpeciesTaxonHash();
# Returns a hash with key equal to the two species name and value equal to the taxon:id found in NCBI.  http://www.ncbi.nlm.nih.gov/Taxonomy/
# Outputs: $taxonHash - described above
# Kristen Naegle
# May 19, 2009
sub returnSpeciesTaxonHash(){
    my %hash;
    $hash{'arabidopsis thaliana'} = 3702;
    $hash{'bos taurus'} = 9913;
    $hash{'caenorhabditis elegans'} = 6239;
    $hash{'chlamydomonas reinhardtii'} = 3055;
    $hash{'danio rerio'} = 7955;
    $hash{'dictyostelium discoideum'} = 44689;
    $hash{'drosophila melanogaster'} = 7227;
    $hash{'escherichia coli'} = 562;
    $hash{'homo sapiens'} = 9606;
    $hash{'mus musculus'} = 10090;
    $hash{'mycoplasma pneumoniae'} = 2104;
    $hash{'oryza sativa'} = 4530;
    $hash{'plasmodium falciparum'} = 5833;
    $hash{'pneumocystis carinii'} = 4754;
    $hash{'rattus norvegicus'} = 10166;
    $hash{'saccharomyces cerevisiae'} = 4932;
    $hash{'schizosaccharomyces pombe'} = 4896;
    $hash{'takifugu rubripes'} = 31033;
    $hash{'xenopus laevis'} = 8355;
    $hash{'zea mays'} = 4577;
    $hash{'gallus gallus'} = 9031;
    $hash{'bacillus subtilis'} = 1423;
    $hash{'oryctolagus cuniculus'} = 9986;
    $hash{'sus scrofa'} = 9823;
    $hash{'canis lupus'} = 9612;

    return(\%hash);
}

# (\%annotationHash, \%versionHash) = returnAnnotationParserHash($speciesArr)
# Given an array of species names (binomial convention), return a hash of annotation providors with keys based on species name
# Inputs: \@species - ref. to array of species of inerest
# Outputs: $hash - hash of annotation parsers, keys are species 
#     example to get homo sapiens parser $parser = $hash->{'homo sapiens'};
#          $verHash - ref. to hash, same species keys but now value is the version of that annotation file
# Kristen Naegle
# May 21, 2009 
# Implemented this way because parser name lookup is taxonomy independent which is a problem if we were to use the large file.
# June 23, 2009 - added version capability
sub returnAnnotationParserHash($){
    my ($speciesArr) = @_;
    my %hash;
    my %verHash;
    my ($rootDir, $rootFile, $oboFile) = returnGlobalGOFileInfo();
    
    foreach my $species (@$speciesArr){
	my $speciesName = $species;
	$speciesName =~ s/ /_/g;
	my $file = $rootDir.$rootFile.$speciesName;
	# check for existence - if file doesn't exist report and go to next species
	if(!-e $file){
	    print "ERROR $file does not exist, therefore $species is not supported\n";
	    next;
	}
	my $annotationParser = GO::AnnotationProvider::AnnotationParser->new(annotationFile=>$file);
	$hash{$species} = $annotationParser;
	my ($release, $f) = returnLinkRelease($file); 
	$verHash{$species} = $release;
    }
    return (\%hash, \%verHash);
}

# ($name, $type) =  returnNameForProteinId($dbh, $proteinId, $annotationParser, $AMBIGUOUS)
# Given a parser and a protein Id, return the non-ambiguous name. Returns the first name that exists and is non-ambigous that is tested for the protein
# Inputs: $dbh - database handle
#         $proteinId - id of protein table
#         $annotationParser - GO parser object of annotation file - assums annotation parser is for correct species DOES NOT CHECK HERE
#         $AMBIGUOUS - 1 if you want to map ambiguous names to all non-ambiguous identifiers, 0 if ignore an ambiguous name
# Outputs: $dbIds - ref. to array of database ids for that protein
#          $names - ref to array of names for protein
# Kristen Naegle
# May 21, 2009
sub returnNamesForProteinId($$$$){
    my ($dbh, $proteinId, $annotationParser, $AMB) = @_;


    my (@dbIds, @names);
    

    # check acc_gene first 
    my  $sth = $dbh->prepare('SELECT acc_gene from protein where id=?');
    $sth->execute($proteinId);
    
    my $gene = returnSingleResultOnCol($sth, 0);

    if($gene){
	print "DEBUG: Looking by acc_gene $gene\n";
	my $goodName = checkName($gene, $annotationParser);
	if($goodName==1){
	    push @names, $gene;
	}
	elsif($goodName==2 and $AMB){
	    my @databaseIds = $annotationParser->databaseIdsForAmbiguousName($gene);
	    push @dbIds, @databaseIds;
	}
    }

    # get accessions for a proteinId
    $sth = $dbh->prepare('SELECT * from acc where protein_id=?');
    $sth->execute($proteinId);
    my $results = returnMultipleResultsForExSTH($sth);
    foreach my $result (@$results){
	my $rType = $result->{'type'};
	if($rType eq 'swissprot'){
	    
	    my $value = $result->{'value'};
	    #print "DEBUG: Looking by name $value\n";
	    my $goodName = checkName($value, $annotationParser);
	    if($goodName==1){
		push @dbIds, $value;
		}
	    }
	elsif($rType eq 'gene_synonym'){
	     my $value = $result->{'value'};
	    #print "DEBUG: Looking by name $value\n";
	    my $goodName = checkName($value, $annotationParser);
	    if($goodName==1){
		push @names, $value;
	    }
	     elsif($goodName==2 && $AMB){
		 my @databaseIds = $annotationParser->databaseIdsForAmbiguousName($value);
		 push @dbIds, @databaseIds;
	     }
	 }
    }
    
    

    return(\@dbIds,\@names); # handle last case
}


# $goodName = checkName($name, $annotationParser);
# checks a name in the annotationParser to see if it exists AND is not ambigous (1) otherwise returns 0.
# Inputs: $name - name to check in database id field, or alias field
#         $annotationParser - annotation parser object of GO
# Outputs: $goodName - 0 if does not exist, 1 exists and is not ambiguous, 2 if exists but is ambiguous
# Kristen Naegle
# May 21, 2009
sub checkName($$){
    my($name, $annotationParser) = @_;
    my $goodName = 0;

    if($annotationParser->nameIsAnnotated(name=>$name)){
	if(!$annotationParser->nameIsAmbiguous($name)){
	    print "$name is not ambigous!\n";
	    $goodName = 1;
	}
	else{
	    my @databaseIds = $annotationParser->databaseIdsForAmbiguousName($name);
	    $goodName = 2;
	    #print "$name was ambigous\n Maps to @databaseIds\n"; } 
	    handleError("GO:checkName", "$name is ambiguous, but is annotated", \@_); 
	}
    }
    return $goodName;
}

# ($rootDir, $rootFile, $oboFile) = returnGlobalGOFileInfo()
# Outputs: $rootDir - location of root directory for GO_DATA files (annotations and ontology files
#          $rootFile - this is the start of teh species specific file name
#          $oboFile - this is the name of the obo file ( a symbolic link to a specific version)
# Kristen Naegle
# June 23, 2009
sub returnGlobalGOFileInfo(){
    my ($rootDir, $rootFile, $oboFile);

    #$rootDir = "/data/knaegle/data/GO/";
    $rootDir = $globalVars::GO_PATH;
    $rootFile = 'gene_association.';
    $oboFile = $rootDir."gene_ontology.obo";
    return ($rootDir, $rootFile, $oboFile);
}

# removeIEATermsFromAnnotationFile($inputFile, $outputFile)
# remove all IEA terms from an annotation file
# Inputs: $inputFile - gene association file from which to remove IEA terms
#         $outputFile - output of same gene association file without IEA terms
# Kristen Naegle
# August 13, 2009
sub removeIEATermsForAnnotationFile($$){
    my ($inputFile, $outputFile) = @_;

    open(GO, $inputFile) || die "Can't open $inputFile for reading\n";
    open(GO_OUT, ">$outputFile") || die "Can't open $outputFile for writing\n";
    my $linesRemoved = 0;
    while(defined(my $line = <GO>)){
	if($line =~ m/^!/){
	    print GO_OUT $line;
	}
	else{
	    my @line = split('\t', $line);
	    my $evidenceCode = $line[6];
	    if($evidenceCode ne 'IEA'){
		print GO_OUT $line;
	    }
	    else{
		$linesRemoved +=1;
	    }
	}
    }
    print "Removed $linesRemoved lines\n";
    close(GO);
    close(GO_OUT);

}

# $fileName = retrieveAnnotationFile($species, $linkLocation)
#Given a binomial species name and the gene ontology consortium location, retreive the file and save it to the data directory.
# Inputs: $species - binomial species name with an underscore
#         $linkLocation - location of species file on gene ontology see GO file gene_association_fileLocations.txt
# Outputs: $fileName - destination of final file, unzipped
# Kristen Naegle
# August 14, 2009
sub retrieveAnnotationFile($$){
    my ($species, $linkLoc) = @_;
    
    my @loc = split('/', $linkLoc);
    my $fileName = $loc[$#loc];
    my $newName = $fileName;
    $newName =~ s/\?rev=HEAD//;
    `wget $linkLoc`;
    
    `mv $fileName $newName`;
    `gunzip $newName`;
    $fileName = $newName;
    $fileName =~ s/\.gz//;
    

    my ($rootDir, $rootFile, $oboFile) = returnGlobalGOFileInfo();
    my ($sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst) = localtime();
    my $relDate = sprintf('rel%02d-%02d-%04d',$mon+1,$mday,$year+1900);
    #$newName = $rootDir.$rootFile.$species."_relDATE";
    my $link = $rootDir.$rootFile.$species;
    $newName = $link."_$relDate";
    `mv $fileName $newName`;
    return ($newName, $link);
}

# $hash = returnSpeciesLinkLocations(); 
# Open and read the link locations file into a hash with keys species and values equal to the wget location on gene ontology
# Inputs - none, contains global link to moving file
# Outputs: $hash - keys species and values link locations
# Kristen Naegle
# August 14, 2009
sub returnSpeciesLinkLocations(){
    my ($rootDir, $rootFile, $oboFile) = returnGlobalGOFileInfo();
    my $fileName = $rootDir.'gene_association_fileLocations.txt';
    open(FILE, $fileName) || die "Can't open $fileName for reading\n";
    my %hash;
    while(defined(my $line = <FILE>)){
	chomp $line;
	my ($species, $loc) = split('\t', $line);
	$hash{$species} = $loc;
    }

    close(FILE);
    return \%hash;
}

# retrieveAllGOAnnotationFiles($REMOVE_IEA)
# Retrieve all annotation files, remove IEA terms if specified and then update symbolic link in the root directory
# Inputs: REMOVE_IEA - bool value specifies whether to remove IEA values from annotation file before linking.  Both files will exist if true, but link will be to non-IEA file
# Kristen Naegle
#August 14, 2009
sub retrieveAllGOAnnotationFiles($){
    my ($REMOVE_IEA) = @_;
    #process
#1. Get a file, remove IEA terms (if specified), symbolic link the file.
    my $hash = returnSpeciesLinkLocations();
    foreach my $species (keys %$hash){
	print "Working on Species: $species with link location: $hash->{$species}\n";
	my ($file, $linkName) = retrieveAnnotationFile($species, $hash->{$species});
	my $newFile = $file;
	if($REMOVE_IEA){
	    #store release number to tack back on 
	    $file =~ m/_rel(.+)/;
	    my $VER = $1;
	    $newFile =~ s/_rel(.+)/noIEA_rel$VER/;
	    removeIEATermsForAnnotationFile($file, $newFile);
	}
	# link symoblic link
	`ln -fs $newFile $linkName`;
    }
}

1;

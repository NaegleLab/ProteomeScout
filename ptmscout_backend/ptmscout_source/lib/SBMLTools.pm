use strict;
use warnings;
use fileTools;


# $arrayRef = returnProtIDArr($proteinFile, $species)
# Takes a tab separated file of proteins (generated first by saving csv of proteins in cell designer and then converting to txt using excel. 
# Returns protein ids for all proteins which contain the name in $species
# Inputs: $proteinFile - tab separated file of proteins in cell designer
#         $species - species name like "EGFR" or "SPRY2"
# Outputs: $arrRef - reference to an array of protein ids that positively matched a species identifier
# Kristen Naegle
# July 31, 2007
sub returnIDArr($$){
    my $proteinFile = shift;
    my $SPECIES = shift;
    
    $SPECIES = lc($SPECIES);
    my $NAME_COL = returnColumnNumber($proteinFile, "name");
    my $PROTID_COL = returnColumnNumber($proteinFile, "id");    
    open(PROTEIN, $proteinFile) || die "Can't open protein file $proteinFile\n";
    my $line = <PROTEIN>; #remove header
    my @IDS;
    my $count = 0;
    while($line = <PROTEIN>){
	chomp $line;
	$count +=1; 
	my @line = split("\t", $line);
	# print "$line[$NAME_COL]\n";
	my $speciesName = lc($line[$NAME_COL]);
	my @species = split("@", $speciesName);

	if(@species){
	    if ($species[0] =~ /$SPECIES/){
		#print "$line[$NAME_COL] \t\t--> \t\t $line[$PROTID_COL]\n";
		push @IDS, $line[$PROTID_COL];
	    } 
	}
	else{
	    if ($speciesName =~ /$SPECIES/){	    
		#print "$line[$NAME_COL] \t\t--> \t\t $line[$PROTID_COL]\n";
		push @IDS, $line[$PROTID_COL];
	    } 
	}
	
    }
    
    #print "There were $count proteins in file\n";
    #print "".scalar(@IDS)." match $SPECIES\n";
    close(PROTEIN);
    return \@IDS;
}

# $reactionHash = createReactionArray($reactionFile, $speciesHashRef)
# Given a boolean species hash, create a boolean reaction hash from a reactionFile 
# Inputs: $reactionFile - the reaction file created by cell designer in csv and converted to txt (could use excel for this)
#         $speciesHashRef - hash with key of species and value being a 1 or 0 as to whether this species includes a protein of interest for the new submodel
# Outputs: $reactionHashRef - reference to hash of reactions contained in reaction file with boolean value as to whether these have a protein of interest in either reactants, products or modifiers
# Kristen Naegle
# July 31, 2007 
sub createReactionArray($$){
    my $reactionFile = shift;
    my $speciesHashRef = shift;
    #check reactants, products and modifiers for species that are allowed 
    my $ID_COL = returnColumnNumber($reactionFile, "id");

    my @COLS;
    push @COLS,returnColumnNumber($reactionFile, "reactants");
    push @COLS, returnColumnNumber($reactionFile, "products");
    #push @COLS, returnColumnNumber($reactionFile, "modifiers");

    my %reactionHash;

    open(REACTIONS, $reactionFile);
    my $line = <REACTIONS>; 
    while($line = <REACTIONS>){
	chomp $line;
	$line =~ s/\"//g;
	my @line = split("\t", $line);
	my $reaction = $line[$ID_COL];
	my @species;
	foreach my $col (@COLS){
	    if($line[$col]){
		my @temp = split(",",$line[$col]);
		push @species, @temp;
		#print "@species\n";
	    }
	    
	}
	foreach my $species (@species){
	    if($speciesHashRef->{$species}){
		$reactionHash{$reaction} = 1;
		print "REACTION $reaction allowed because of $species\n";
		last;
	    }
	    else{
		$reactionHash{$reaction} = 0;
		#print "Reaction $reaction NOT ALLOWED\n";
	    }
	}



    }
    close(REACTIONS);
    return \%reactionHash;

}

# removeReactionsFromSBML($SBMLfile, $outputFile, $reactionHashRef)
# Takes a reaction hash reference and an SBML file and rewrites the SBML file without those reactions that have a 0 in the hash 
# Inputs: $SBMLfile - SBML file
#         $outputFile - output SBML file
#         $reactionHashRef - reference to a reaction hash, where keys are reactions and values are 0 or 1 depending on whether they should be included
# Outputs: Outputs a modified file to the output file
# August 3, 2007
# Kristen Naegle
sub removeReactionsFromSBML($$$){
    my $SBMLIn = shift;
    my $SBMLOut = shift;
    my $reactionHashRef = shift;

    open(SBMLIN, $SBMLIn)||die "Can't open input SBML file\n";
    open(SBMLOUT, ">$SBMLOut") || die "Can't open output SBML file $SBMLOut for writing\n";

    my $line;
    while($line = <SBMLIN>){
	my $FLAG_RE = 0;
	#for now just remove <reaction id="x"> ... </reaction> sections of XML
	#Full SBML section of reaction id equivalents
		## Handle Catalyzed Reaction lines (remove the single line)
	if($line =~ /reaction=/){
	    my $re = getReactionNum($line);
	    print "FOUND CATALYSIS: $re\n";
	    if($reactionHashRef->{$re}){
		print "ALLOWING Catylsis reaction $re\n";
	    }
	    else{
		print "REMOVING Catalysis reaction $re\n";
		$line = <SBMLIN>;
	    }
	}
	if($line =~ /<reaction id=/){
	    my $re = getReactionNum($line);

	    #print "$re\n";
	    if($reactionHashRef->{$re}){
		print "Allowing $re\n";
		
	    }
	    else{
		$FLAG_RE = 1;
		print "Removing $re\n";
		#print "NOT PRINTING $line";
		while(not $line =~ /reaction>/){
		    if($re eq "re138"){
			print "NOT PRINTING IN WHILE  $line";
		    }
		    $line = <SBMLIN>;		    
		    
		} #end while 
		
		if($re eq "re138"){
		    print "NOT PRINTING AFTER $line\n";
		}

	    } #end check for reaction in hash

	       
	}
#skip printing line </reaction> but don't move onto the next b/c then the new $line = <SBMLIN> at while will skip something
	if(!$FLAG_RE){
	    print SBMLOUT $line;
	}
	else{
	    #skip line and reset flag
	}

    } #end file 

    close(SBMLIN);
    close(SBMLOUT);

}

sub getReactionNum($){
    my $line = shift;
    my @rxn = split(/\"/, $line);
    return $rxn[1];
    
}


# colorSpeciesInSBML($SBMLfile, $outputFile, $speciesHashRef)
# Takes a species hash reference and an SBML file and rewrites the SBML file without those reactions that have a 0 in the hash 
# Inputs: $SBMLfile - SBML file
#         $outputFile - output SBML file
#         $reactionHashRef - reference to a reaction hash, where keys are reactions and values are 0 or 1 depending on whether they should be included
# Outputs: Outputs a modified file to the output file
# October 8, 2007
# Kristen Naegle
sub colorSpeciesInSBML($$$$){
    my $SBMLIn = shift;
    my $SBMLOut = shift;
    my $speciesHashRef = shift;
    my $color = shift;
    my @SPECIES;
    open(SBMLIN, $SBMLIn)||die "Can't open input SBML file\n";
    open(SBMLOUT, ">$SBMLOut") || die "Can't open output SBML file $SBMLOut for writing\n";

    my $line;
    while($line = <SBMLIN>){
	if(($line =~ m/speciesAlias id=/gi) or ($line =~ m/speciesAlias compartmentAlias/gi)){
	    #print "FOUND LINE: $line";
	    #print SBMLOUT $line;
	    $line =~ m/species=(.+?)>/gi;

	    my $species = $1;
	    $species =~ s/\"//g;
	    push @SPECIES, $species;

	    if($speciesHashRef->{$species}){
		
		print "Species found: $species\n";
#		while(not($line=~ m/paint color=/gi)){
		while(not($line=~ m/celldesigner:speciesAlias>/gi)){
		    if($line =~ m/paint color=/gi){
			my $colorLine = '<celldesigner:paint color="'.$color.'" scheme="Color"/>'."\n";
			print SBMLOUT $colorLine;
		    }
		    else{print SBMLOUT $line;}
		    #print $line;
		    $line = <SBMLIN>;
		}
		print SBMLOUT $line;
		#print "JUST outside while LINE: $line";
		
#		my $colorLine = '<celldesigner:paint color="ffffff66" scheme="Color"/>'."\n";
		#print SBMLOUT $colorLine;
		#print $colorLine;
		
	    }
	    else { #FOUnd a speciesAlias line but it's not a species in the hash
		print "FOUND $species but not allowed\n";
		print SBMLOUT $line;

	    }
	    

	} # end found start of a definition for shape 
	else { print SBMLOUT $line;}
    
    
    } #end file 

    foreach my $species (@SPECIES){
	#print "Found $species in <complexSpeciesAliasId>\n";
    }
    close(SBMLIN);
    close(SBMLOUT);

}

sub returnSpeciesHash($$){
    my $proteinFile = shift;
    my $SPECIESRef = shift;
    my @IDS_all;
    foreach my $species (@$SPECIESRef){
	my $idRef = returnIDArr($proteinFile, $species);
	my @IDS = @$idRef;
	push @IDS_all, @IDS;
#print "There were $count proteins in file\n";
	#print "".scalar(@IDS)." match $species\n";
    }
    #print "TOTAL PROTEINS: ".scalar(@IDS_all)."\n";
#print @IDS_all;
#there are potentially redundant entries since two could be in complex, etc.
    my $IDS_uniq = createUniqHash(\@IDS_all);
    return $IDS_uniq;
}

1;

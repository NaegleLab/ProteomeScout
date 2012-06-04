use strict;
use warnings;
use fileTools;

# $domainHash = returnBindHash($file, $ligandHash)
# Return the domainHash, has two dimensional structure, the hash has keys that correspond to the PDZ domain, and the value is an array for all the peptides tested against it with values of the binding affinity.  0 means an array negative, 1 means a FP negative and every other binding affinity values are in nM. $domainHash->{'peptideArray'}->[$index] gives the name of the peptide tested at that index
# Inputs: $file - the binding file (2dimensional, domain x peptide)
#         $ligandHash - a hash of peptide names
# Outputs: $domainHash - returns the reference to a domain hash
# Kristen Naegle
# December 18, 2007
sub returnBindHash($$){
    my $file = shift;
    my $ligandHash = shift;
    
#first line is the protein so create this hash..must 
#check that these are all in $ligandHash
    #have to convert names to peptides 
    
    open(FH_BIND, $file) || die "Can't open $file for reading in $0\n";
    my $line = <FH_BIND>;
    chomp $line;
    my @ligands = split("\t", $line);
    my @errArr;
    my %pepHash;
    my %domainHash;
    push @{$domainHash{'peptideArray'}}, @ligands;
    my @ligandPep;
    foreach my $i (1..$#ligands){
	my $ligandName = $ligands[$i];
	$ligandName =~ s/ //g;
	my $ligandPep = $ligandHash->{$ligandName}->[0];
	if(not $ligandHash->{$ligandName}){
	    push @errArr, $ligandName;
	    print "ERROR: $ligandName in $file not found in ligand hash\n";
	}
	else{
	    push @ligandPep, $ligandPep;
	}
	
    }

    
    #now contine with the domain lines and create an array that corresponds to the pepArray
    while(defined($line = <FH_BIND>)){
	chomp $line;
	if(defined($line)){
	    my @line = split("\t", $line);
	    my $key = $line[0];
	    $key =~ s/ //g;
	    #print "Domain Hash Key: $key\n ";#$line\n";
	    if(not $domainHash{$key}){
		$domainHash{$key} = [];
	    }

	    my $last = $#line;
	    my @array;
	    push @{$domainHash{$key}}, @line[1..$last];
	    #$domainHash{$key} = $line[1..$last];
	}
    }
    close(FH_BIND);
    return \%domainHash;
   
}

# $pepSequence = convertNameToPeptide($transHash, $nameArray)
# Return the aa sequence for the given peptide names
# Inputs: $transHash - hash with peptide names as keys and sequences as values
#         $nameArray - ref to array of names desired for translation
# Outputs: $pepSequence - corresponding array of sequences that match name array
# Kristen Naegle
# December 18, 2007
sub convertNameToPeptide($$){
    my $transHash = shift;
    my $nameArray = shift;
    my @pepArray; 
    my @errArray;

    foreach my $name (@$nameArray){
	$name =~ s/ //g;
	if($transHash->{$name}){
	    push @pepArray, $transHash->{$name}->[0];
	}
	else{
#	    print "ERROR in $0: Could not find $name in translateHash\n";
	    push @errArray, $name;
	}
    }
    return (\@pepArray, \@errArray);
}

# $domainKeys = returnAllDomainKeys($domainHash);
# Return all the domain keys in a domain hash (excluding the peptideArray)
# Inputs: $domainHash - reference to domain hash
# Outputs: $domainKeys - reference to array of domains contained in hash
# Kristen Naegle
# December 18, 2007
sub returnAllDomainKeys($){
    my $domainHash = shift;
    
    my @domainKeys;
    foreach my $k (keys %$domainHash){
	print "Key: $k\n";
	if(!($k =~ m/array/gi)){
	    push @domainKeys, $k;
	}
	#else{
	   # print "FOUND PEPTIDE ARRAY key, and skipped it\n";
	#}
    }
    
    return \@domainKeys;
}

# $pepName = returnPeptideNameAtIndex($domainHash, $index)
# Given an index return the peptide name that refers to it (uses domainHash{peptideArray}
# Inputs: $domainHash - reference to domainHash 
#         $index - zero based index
# Outputs: $pepName - string that describes peptide (to get actual peptide sequence must translate through a key, convertNameToPeptide)
# Kristen Naegle 
# December 18, 2007
sub returnPeptideNameAtIndex($$){
    my $domainHash = shift;
    my $index = shift;
    my $key = "peptideArray";
    my $peptide = $domainHash->{$key}->[$index];
    return $peptide;
}

# $pepArray = returnPepArrayBelowCutoff($domainHash, $domainKeys, $threshold)
# Given a Hash of domains, an array of domain keys to look up and an affinity threshold, return all the peptides that bind to those domains with greater than the given affinity
# Inputs: $domainHash - hash of domains, where $domainHash->{'peptideArray'} gives the array reference to the peptides and the remainder are affinity values that line up with the peptideArray
#         $domainKeys - array of keys for which to return binders
#         $threshold - affinity threshold that matches value in domainHash (i.e. for PDZ it's in nM
# Outputs: $pepArray - array of peptides that bound to desired domains ($domainKeys) below the threshold
# Kristen Naegle
# December 18, 2007
sub returnPepArrayBelowCutoff($$$){
    my $domainHash = shift;
    my $domainKeys = shift;
    my $threshold = shift;
    my @pepArray;
    
     if($domainKeys->[0] =~ m/all/){
 	$domainKeys = returnAllDomainKeys($domainHash);
	print "Retrieving all domains in hash\n";
     }

    foreach my $key (@$domainKeys){
	if(not $domainHash->{$key}){
	    print "ERROR: Your domain $key was not found in domain hash in $0\n";
	}
	my $l = scalar(@{$domainHash->{$key}}) - 1;
	#print "The array length is $l\n";
	foreach my $i (0..$l){
	    #print $i."\n";
	    my $affinity = $domainHash->{$key}->[$i];
	    #print "value: ".$affinity."\n";
 	    if(($affinity > 0) && ($affinity < $threshold)){
 		my $peptideName = returnPeptideNameAtIndex($domainHash, $i);
		#print $domainHash->{$key}->[$i];
		push @pepArray, $peptideName;
 	    } 
	}
    }
    
    return \@pepArray;
}


# $hashRef = createUniqHash($arrRef)
# Takes a non-unique array and creates a unique hash
# Inputs: $arrRef - reference to an array
# Outputs: $hashRef - reference to a uniq hash based on array items
# July 31, 2007
# Kristen Naegle
sub createUniqHash($){
    my $arr = shift;
    my %hash;
    my @uniq;
    foreach my $item (@$arr){
	$hash{$item} = 1;
    }
    return \%hash;
}



# $array = hashToArr($hashRef);
# convert hash to array of keys
# Inputs: $hashRef - reference to hash
# Outputs: $array - referenc to array of keys
# Kristen Naegle
# December 18, 2007
sub hashToArr($){
    my $hash = shift;
    my @array;
    foreach my $key (keys %$hash){
	push @array, $key;
    }
    return \@array;
}

1;

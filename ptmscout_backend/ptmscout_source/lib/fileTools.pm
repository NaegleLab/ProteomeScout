#! /usr/bin/perl
use warnings;
use strict;
use proteinArrayTools;
use errorHandling;

use constant ACC => "acc:";
use constant NAME => "name";
use constant PEPTIDE => "pep";
use constant PEPTIDE_ALIGNED => "aligned";
use constant DATA_FIELD => "data:"; 

# \@colNums = returnColumnNumberArr($inputFile, $field)
# Returns an array of column numbers (zero-based) from header in a tab dilineated data file
# e.g. If you want all data labels then call returnColumnNumber($file, "data:");
#else if you want a specific label then call returnColumnNumber($file, "data:t0");
#inputs: $inputFile - tab seperated file where the first row is the header info
#        $field - label descriptor that you are looking for
#outputs: $col - reference to an array of column numbers that matched that descriptor
# Kristen Naegle
# July 09, 2007
sub returnColumnNumberArr($$){
    my $inputFile = shift;
    my $field = shift;
    $field = lc($field);
    open(FH_3, $inputFile) || die "Can't open input file $inputFile in returnColumnNumber:fileTools.pm";

    my $line = <FH_3>; #header is always first line
    close(FH_3);
    chomp $line;
    my @line = split("\t", $line);
    my $i;
    my @col;
    #print "$field\n";
    $field =~ s/'"'//g; #have to strip quotes...can't regex compare
    for($i=0; $i <= $#line; $i++){
	#print "$line[$i]\n";
	#if(lc($line[$i]) =~ /lc($field)/){
	if(lc($line[$i]) =~ /$field/ or $field eq "all"){
	    push @col, $i;
	    #print "You are getting column $i with label $line[$i]\n";
	}
    }
    if(not @col){
	push @col, -1; #indicates error
	#print "ERROR: Could not find column..perhaps your file does not have a header? See fileTools.pm for header types\n";
    }
    close(FH_3);
    return \@col;
    
}

# $header = printHeader($file);
# Print the header text in a file to see what the labels are
# Inputs: $inputFile - tab seperated file where the first row is header text
# Outputs: $header - string of header row (can be split on tabs)
#          PRINTS Header to stdout
#Kristen Naegle
# July 09, 2007
sub printHeader($){
    my $inputFile = shift;
    open(FH_2, $inputFile) || die "Can't open input file $inputFile in returnColumnNumber:fileTools.pm";

    my $line = <FH_2>; #header is always first line
    close(FH_2);
    print "HEADER ROW in $inputFile: $line\n";
    return $line;
}
# $header = returnHeader($file);
# Returns the header text in a file to see what the labels are
# Inputs: $inputFile - tab seperated file where the first row is header text
# Outputs: $header - string of header row (can be split on tabs)
#Kristen Naegle
# July 09, 2007
sub returnHeader($){
    my $inputFile = shift;
    open(FH_2, $inputFile) || die "Can't open input file $inputFile in returnColumnNumber:fileTools.pm";

    my $line = <FH_2>; #header is always first line
    close(FH_2);
    return $line;
}

# $headerArrRef = returnHeader($file);
# Return the header array for when you are looking at various versions of header types (e.g. acc but want to know if there is acc:gi and acc:sp)
# Inputs: $inputFile - tab seperated file where the first row is header text
# Outputs: $headerArrRef - array of header items 
#Kristen Naegle
# September 13, 2007
sub returnHeaderArr($){
    my $inputFile = shift;
    open(FH_2, $inputFile) || die "Can't open input file $inputFile in returnColumnNumber:fileTools.pm";

    my $line = <FH_2>; #header is always first line
    close(FH_2);
    my @line = split("\t", $line);
    return \@line;
}

# $missingHeadersArr = headerCheck($inputFile, $headerColsArr)
# Returns an array of items from $headerColsArr that were not in the header row
# Inputs: $inputFile - the tab separated input file with a header
#         $headerColsArr - reference to an array of items that you are checking for in header
# Outputs: $missingHeadersArr - reference to arry of items in input that didn't match the header row
# Kristen Naegle
# October 24, 2007
sub headerCheck($$){
    my $inputFile = shift;
    my $headerArrRef = shift;
    my @missingHeaderCol;
    foreach my $headerCol (@$headerArrRef){
	my $colArr = returnColumnNumberArr($inputFile, $headerCol);
	if($colArr->[0] == -1){
# 	    if(!(defined($missingHeaderCol[0]))){
# 		@missingHeaderCol = []; 
# 	    }
	    push @missingHeaderCol, $headerCol;
	    print "$headerCol\n";
	}

    }

    return \@missingHeaderCol; 

} 

# $column = createSingleColumn($arrayRef);
# From an array of column numbers return the only element in a single array or tell the user that they have too many choices, prints the header options
# Inptus: $arrayRef - reference to array of potential column numbers
# Outputs: $column - column number when there was only one element in array, else 0
# Kristen Naegle
# July 12, 2007
sub createSingleColumn($){
    my $arrayRef = shift;
    if(scalar(@$arrayRef) > 1){
	#print "You're column array has multiple entries. Choose a more specific one to convert ot a single column\n";
	#print "You're headers:\n";
	foreach my $item (@$arrayRef){
	    #print "$item\n";
	    return -2;
	}

    }
    else{
	return $arrayRef->[0];
    }
}

# $colNumber = returnColumnNumber($file, $field)
# Returns a single column number for a field in $file with a header row
# Essentially gets an array first and then checks for one or more elements
# Errors: -2 means there was more than one header that matched (USE a more specific field name)
#         -1 means field did not exist in header
# Output: $colNumber - zero based column number for $field
# Kristen Naegle
# 7/19/07
sub returnColumnNumber($$){
    my $inputFile = shift;
    my $field = shift;
    my $arrRef = returnColumnNumberArr($inputFile, $field);
    my $colNumber = createSingleColumn($arrRef);
    return $colNumber;
}

# $throw = catchColumnError($colNameArrRef, $colNumArrRef, $scriptName);
# Returns a throw (1) if the any of the columns returned a negative number and prints the exact error 
# Inputs: $colNameArrRef - ref to array of names of the columsn
#         $colNumArrRef - array of col numbers (-1 and -2 are error codes)
#         $scriptName - name of script this sub was called from for stack catch use $0
# Outputs: $throw - returns 1 if error and 0 if none
# Kristen Naegle
# 12/5/07
sub catchColumnError($$$){
    my $colNameArr = shift;
    my $colNumArr = shift;
    my $scriptGen = shift;
     
    my $throw = 0;
    my $i;
    for ($i=0; $i < scalar(@$colNumArr); $i++){
    if($colNumArr->[$i] < 0){
	my $colName = $colNameArr->[$i];
	print "ERROR in $scriptGen: COLUMN HEADER $colName";
	if($colNumArr->[$i] == -1){
	    print " was not found\n";
	    $throw = 1;
	    
	}
	elsif($colNumArr->[$i] == -2) {
	    print " was not specified enough - multiple matches\n";
	    $throw = 1;
	    
	}
    }
 
    
}
return $throw;
}


# # # $accCol = returnAccCol($dataFile);
# # # Return the accession column - chooses in order of preference GI (gi), SwissProt (sp), RefSeq (ref), Gene Symbol (gene)
# # # Inputs: $inputFile - datafile with header
# # # Outputs: $accCol - accession column (zero-based) of accession number with highest priority, if non exist it returns -1;
# # # Kristen Naegle
# # # September 13, 2007
# # sub returnAccCol($){
# #     my $inputFile = shift;
# #     my @accArrChoice = ("gi", "sp", "ref", "gene");
# #     my $i = 0;
# #     my $accArrRef = returnColumnNumber($inputFile, "acc:".$accArrChoice[$i]);
# #     while ($accArrRef < 0){
# # 	if($i > $#accArrChoice){
# # 	    print "ERROR: Could not find an accession number in your file. Please ma
# # ke sure your header includes an acc:gi, acc:sp, or acc:ref\n";
# # 	    return -1;
# # 	    exit;
# # 	}
# # 	$i += 1;
# # 	$accArrRef = returnColumnNumber($inputFile, "acc:".$accArrChoice[$i]);
# #     }
    
# #     return $accArrRef;
# # }

# $accCol = returnAccCol($dataFile);
# Return the accession column - column is annotated with acc
# Inputs: $inputFile - datafile with header
# Outputs: $accCol - accession column (zero-based) of accession number, returns -1 if accession col doesn't exist and -2 if more than one column is labeled with keyword acc.
# Kristen Naegle
# July 7, 2009 - rewrote to handle just acc, don't care about all of them having the same type
sub returnAccCol($){
    my $inputFile = shift;
   
    my $accCol = returnColumnNumber($inputFile, "acc");
    if($accCol == -1){
	print "ERROR: Could not find an accession number in your file. Please ma
ke sure your header includes an acc:gi, acc:sp, or acc:ref\n";
	return $accCol;
    }
    
    
    return $accCol;
}

# $pepCol = returnPeptideColumn($dataFile)
# Given a data file, return the column number (zeros based) of the peptide
# Inputs: $dataFile - tab separated data file
# Outputs: $pepCol - zeros based column number where peptide resides (indicated by pep:tryps, pep, pep:aligned.  If more than one of these fields exists than it will return pep:tryps then pep:aligned in order of priority
# Kristen Naegle
# July 7, 2009
sub returnPeptideCol($){
    my($dataFile) = @_;
    
    my $pepCol = returnColumnNumber($dataFile, 'pep');
    if($pepCol == -2){
	my $pepTryps = returnColumnNumber($dataFile, 'pep:tryps');
	if($pepTryps == -1){
	    $pepCol = returnColumnNumber($dataFile, 'pep:aligned');
	    
	}
	else{
	    $pepCol = $pepTryps;
	}
    }
    if($pepCol == -1){
	print "Error: No peptide column exists\n";
    }
    
    return $pepCol;

}

# $acc = returnAcc($line, $colNumber);
# returns a stripped accession number (number only)
# Inputs: $line from input file, tab dileneated
#         $colNumber - zeros based column number corresponding to desired accession
# Outputs: $acc - the accession number stripped of descriptors (gi|, NP_) and spaces
# Kristen Naegle
# October 24, 2007
sub returnAcc($$){
    my $line = shift;
    my $accCol = shift;
    chomp $line;
    my @line = split("\t", $line);
    my $acc = $line[$accCol];
  #  $acc =~ s/gi\|//;
  #  $acc =~ s/NP_//;
    $acc =~ s/ //g;
    return $acc;
}

# Given a full accession description, return a stripped version with only the number
sub returnAccNumber($){
    my $acc = shift;
    $acc =~ s/[A-Z]//gi;
    $acc =~ s/ //g;
    $acc =~ s/\|//g;
    $acc =~ s/_//g;
    return $acc;

    
}



# \%hash = returnUniqColumnHash($dataFile, $colNumber);
# Returns all the uniq entries in a file in a particular column, the value of the hash is the number of times that entry appeared. This includes taking multiple fields from a column
# Inputs: $inputFile - tab separated file from which you want to gather a uniq hash
#         $colNumber - column number from which you want to get hash
# Outputs: $hashRef - a hash with uniq keys, which are taken from ; and : separated fields, value is number of times that appeared in the column
# October 15, 2007 
# Kristen Naegle
sub returnUniqColumnHash($$){
    my $inputFile = shift;
    my $colNumber = shift; # zeros-based
    my %hash;
    open(D3, $inputFile) || die "Can't open $inputFile for reading in returnColumnHash (fileTools.pm)\n";
    my $line = <D3>; # eat the header
    while(defined($line = <D3>)){
	#chomp $line; 
	#my @line = split("\t", $line);
	#my $key = $line[$colNumber];
	my $keys = returnAllFields($line, $colNumber);
	foreach my $key (@$keys){
##Eventually need to handle multiple keys that come from ; separations
	    if(not $hash{$key}){
		$hash{$key} = 0;
	    }
	    $hash{$key} +=1; 
	}
    }
    close(D3);
    return \%hash;


}

# \%hash = returnUniqColumnHashWithVal($dataFile, $keyColNumber, $valColNumber);
# Returns all the uniq entries in a file in a particular column, the value of the hash is the number of times that entry appeared. This includes taking multiple fields from a column
# Inputs: $inputFile - tab separated file from which you want to gather a uniq hash
#         $colNumber - column number from which you want to get hash
#         $valColNumber - column number from which you want to get value
# Outputs: $hashRef - a hash with uniq keys, which are taken from ; and : separated fields, value is from the field passed in by $valColNumber
# December 17, 2007 
# Kristen Naegle
sub returnUniqColumnHashWithVal($$$){
    my $inputFile = shift;
    my $colNumber = shift; # zeros-based
    my $valCol = shift;
    my %hash;
    open(D3, $inputFile) || die "Can't open $inputFile for reading in returnColumnHash (fileTools.pm)\n";
    my $line = <D3>; # eat the header
    while(defined($line = <D3>)){
	#chomp $line; 
	#my @line = split("\t", $line);
	#my $key = $line[$colNumber];
	my $keys = returnFields($line, $colNumber);
	my $vals = returnFields($line, $valCol);
	foreach my $key (@$keys){
##Eventually need to handle multiple keys that come from ; separations
	    if(not $hash{$key}){
		$hash{$key} = 0;
	    }
	    $hash{$key} = $vals; 
	}
    }
    close(D3);
    return \%hash;


}

# \%hash = returnUniqColumnHashArrWithVal($dataFile, $keyColNumber, $valColNumber);
# Returns all the uniq entries in a file in a particular column, the value of the hash is the number of times that entry appeared. This includes taking multiple fields from a column. Pushes into array when multiple keys are found
# Inputs: $inputFile - tab separated file from which you want to gather a uniq hash
#         $colNumber - column number from which you want to get hash
#         $valColNumber - column number from which you want to get value
# Outputs: $hashRef - a hash with uniq keys, which are taken from ; and : separated fields, value is from the field passed in by $valColNumber
# May 27, 2008 
# Kristen Naegle
sub returnUniqColumnHashArrWithVal($$$){
 my $inputFile = shift;
    my $colNumber = shift; # zeros-based
    my $valCol = shift;
    my %hash;
    open(D3, $inputFile) || die "Can't open $inputFile for reading in returnColumnHash (fileTools.pm)\n";
    my $line = <D3>; # eat the header
    while(defined($line = <D3>)){
	#chomp $line; 
	#my @line = split("\t", $line);
	#my $key = $line[$colNumber];
	my $keys = returnFields($line, $colNumber);
	my $vals = returnFields($line, $valCol);
	foreach my $key (@$keys){
##Eventually need to handle multiple keys that come from ; separations
	    if(not $hash{$key}){
		$hash{$key} = [];
	    }
	    push @{$hash{$key}}, $vals; 
	}
    }
    close(D3);
    return \%hash;


}

# \@fields = returnFields($line, $col)
# Returns reference to an array of fields contained in a column of the tab separated line
# splits fields on ;  and then takes the desriptor name as the first in a : separated descriptor
# Inputs: $line - tab separated line
#         $col - column number
# Outputs: \@fields - array of fields
# October 26, 2007
# Kristen Naegle
sub returnFields($$){
    my $line = shift;
    chomp $line;
    my $col = shift;
    my $field = returnField($line, $col);
    my @fields;
    if($field =~ /;/){
	my @fieldT = split(";", $field);
	foreach my $item (@fieldT){
	    push @fields, returnFieldDescriptor($item);
	}
    }
    else{
	push @fields, returnFieldDescriptor($field);
    }
    return \@fields;


}

# \@fields = returnAllFields($line, $col)
# Extends returnFields to split on comma separated as well
# Inputs: $line - the line of a tab separated data file
#         $col - column number of that file of interest 
# Outputs: \@fields - array of all fields, given by field descriptor, separated by ; or ,
# Kristen Naegle
# January 8, 2008
sub returnAllFields{
    my $line;
    my $col;
    my $uniqueOnly;
    if(scalar(@_) == 3){
	($line, $col, $uniqueOnly) = @_;
    }
    elsif(scalar(@_) ==2){
	($line, $col) = @_;
	$uniqueOnly = 0;
    }
    else{
	print "ERROR: Incorrect number of arguments in returnAllFields in fileTools.pm, only 2 or optional third argument allowed\n";
	exit;
    }
    my $field = returnFields($line, $col);
    my @fields;
    $field =~ s/\"//g;
       
    foreach my $f (@$field){
	$f =~ s/\"//g;
	my @f;
	if($f =~ /,/){
	    @f = split(',', $f);
	    push @fields, @f;
	}
	else{
	    push @f, $f;
	    push @fields, $f;	    
	}
    }
    my $fref;
    if($uniqueOnly){
	$fref = uniq(\@fields);
    }
    else{
	$fref = \@fields;
    }
    return $fref;
}

# $descriptor = returnFieldDescriptor($field);
# Returns the first item of a : separated $field
# Inputs: $field - as in what comes from returnField that has been split on ;
# Outputs: $descriptor - the first field in the : seperated field
# October 26, 2007
# Kristen Naegle
sub returnFieldDescriptor($){
    my $field = shift;
    my $descriptor;
    if($field =~ /:/){
	my @d = split(":", $field);
	$descriptor = $d[0];
    }
    else{
	$descriptor = $field;
    }
    return $descriptor;
}

# $field = returnField($line, $col)
# returns the field in the column of tab-spearated $line
# Inputs: $line - tab separated line
#         $col - column in line requested
# Outputs: $field - space removed field
# October 26, 2007
# Kristen Naegle
sub returnField($$){
    my $line = shift;
    my $col = shift;
    chomp $line;
    my @line = split("\t", $line);
    my $field = $line[$col];
    $field =~ s/ //g;
    return $field;
}

# [$dirSource, $fileName, $ext] = returnFileNames($)
# Splits up a file name into a root a name and an extension
# Inputs: $fileName
# Outputs: 
sub returnFileNames($){
    my $file = shift;
    my @root = split('/', $file);
    my $fileName = $root[$#root];
    my @ext = split('\.', $file); #need to handle no extensions and no root directory
    my $ext = '.'.$ext[$#ext];
    $fileName =~ s/$ext//;
    my $i;
    my $root;
    if(scalar(@root) < 1){
	$root = "./"; #actually, this should be a get current path..probably this is a special perl var..but actually, the rest assumes dependency as well, so ./ is probably fine
    }
    else{
	foreach $i (0..$#root-1){
	    $root .= $root[$i]."/";
	}
    }
    return($root, $fileName, $ext);
    
}

# expressHeaderTypes;
# prints the possible header options accepted in data files 
# Kristen Naegle
# July 12, 2007
sub expressHeaderTypes{
    print "Here are the standard header types. If you use the most general case returnColumnNumber will return an array with all types\n\n";
    print "ACCESSIONS:\t\t\t\tacc:\t\n";
    print "\t\t\tGI\t\tacc:gi\n";
    print "\t\t\tSwissprot\t\tacc:sp\n";
    print "\t\t\tRefSeq\t\tacc:ref\n";
    print "\t\tGene Symbol\t\tacc:gene\n";
    print "\tUnique Identifier\t\tacc:peptide\n";
    
    print "\n\n";
    print "DATA:\t\t\t\tdata:\t\n";
    print "\t\tTimes\t\tdata:time=LABEL\n";
    print "\t\tHeader\t\tdata:LABEL\n";
    
    print "\n\n";
    print "PEPTIDES: \t\t\t\tpep\n";
    print "\t\t\tMS Fragment:\tpep\n";
    print "\t\t\tAligned: \tpep:aligned\n";
    
    print "\n\n";
    print "NAMES: \t\t\t\t\tname\n";
    print "\t\t\tLong Name\tname:long\n";
    print "\t\t\tShort Name\tname:short\n";
    
    print "\n\n";
    print "SITES: \t\t\t\t\tsite:\n";
    print "\t\t\tSite Number\t\tsite:number\n";
    print "\t\t\tPhosphorylation Type\tsite:type\n";
    
    print "\n\n";
    print "Notes: \t\t\t\tnote:LABEL\n";

    print "\n\n";
    print "Empty column: \t\t\t\tempty\n";


}

# $array = returnLineArrayOfMatch($file,$argsHashArray);
# returns an array of file lines that match the array of criteria
# Inputs: $file - tab separated file with column header
#         $argsHashArray - array of hash objects, hash has a criteria, regexp, and columnHeader field (regexp value is 1 or 0);
# Outputs: $array - reference to array of lines that matched the intersection of all the criteria
# Kristen Naegle
# January 31, 2008
sub returnLineArrayOfMatch($$){
    my $file = shift;
    my $args = shift;
    
    my $hashRef = $args->[0];
    my $criteria = $hashRef->{'criteria'};
    my $column = $hashRef->{'column'};
    my $regexp = $hashRef->{'regexp'};
    my $columnNumber = returnColumnNumber($file, $column);
    my $throw = catchColumnError([$column], [$columnNumber], "returnLineArrayOfMatch");
    if($throw){exit;}
    my $fileHash = returnLineHash($file, $columnNumber);
    my $lineHash = returnSingleCriteriaHash($fileHash, $criteria, $regexp); #this is the starting to a hash of lines 

    my $numberArgs = scalar(@$args);
    # Get the first criteria hash (hash of lines)
    foreach my $i (1..$numberArgs-1){
	my @array;
	$hashRef = $args->[$i];
	$criteria = $hashRef->{'criteria'};
	$column = $hashRef->{'column'};
	$regexp = $hashRef->{'regexp'};
	$columnNumber = returnColumnNumber($file, $column);
	$throw = catchColumnError([$column], [$columnNumber], "returnLineArrayOfMatch");
    if($throw){exit;}
	my $fileHash = returnLineHash($file, $columnNumber);
	my $lineHashTemp = returnSingleCriteriaHash($fileHash, $criteria, $regexp);
	foreach my $line (keys %$lineHash){
	    if(not $lineHashTemp->{$line}){
		delete $lineHash->{$line};
	    }
	}
	
    }
    #my $arrayUniq = uniq(\@array);
    #convert hash keys of lines to an array
   # my @array = keys %$lineHash;
    return $lineHash;
}

# $hashRef = returnSingleCriteriaHash($fileHash, $criteria, $regexp)
# Given a hash of lines with keys that are the column values, return the subset of that hash that matches a criteria (either via regexp or exact match
# Inputs: $fileHash - hash with keys equal to column value and value equal to array of lines of that column value
#        $criteria - string value for criteria match
#        $regexp - boolean, 1 if you want to check criteria based on regexp 0 if must be equal to
# Outputs: $hashRef - keys now are lines that had field and criteria match and value is number of times that line appeared in $fileHash
# Kristen Naegle
# January 31, 2008
sub returnSingleCriteriaHash($$$){
    my $fileHash = shift;
    my $criteria = shift;
    my $regexp = shift;
    my %lineHash;

    if($regexp){
	foreach my $key (keys %$fileHash){
	    if($key =~ /$criteria/){
		my $values = $fileHash->{$key};
		foreach my $value (@$values){
		    if(not $lineHash{$value}){
			$lineHash{$value} = [];
		    }
		    push @{$lineHash{$value}}, $key;
		    #print "putting $key in hash\n";
		    my @temp = @{$lineHash{$value}};
		    #print "value of hash: ".$temp[$#temp]."\n";
		}
	    }
	}
    }
    else{
	my $values = $fileHash->{$criteria};
	foreach my $value (@$values){
	    if(not $lineHash{$value}){
	    $lineHash{$value} = [];
	    }
	    push @{$lineHash{$value}}, $criteria;
	}
    }
    return \%lineHash;
    

}

# $hashRef = returnCriteriaHash($column, $criteria, $regexp);
# Converts a column header, criteria condition, and boolean regexp to a hash object, for use in args development
# Inputs: $column - string of column header 
#         $criteria - string for criteria comparison
#         $regexp - boolean, 1 if you want to check criteria based on regexp 0 if must be equal to
# Outputs: $hashRef KEYS: 'column', 'criteria', 'regexp'
# Kristen Naegle
# January 31, 2008
sub returnCriteriaHash($$$){
    my $column = shift;
    my $criteria = shift;
    my $regexp = shift;
    
    my %hash; 
    $hash{'criteria'} = $criteria;
    $hash{'column'} = $column;
    $hash{'regexp'} = $regexp;

    return \%hash;

}

# $hash = returnLineHash($file, $keyCol)
# Returns a hash with key according to $keyCol and value is line of file
# Inputs: $file - tab separated pELM file
#         $keyCol - column number to use for key (i.e. pep:tryps)
# Outputs: $hash - hash reference, keys have replaced pY with y, etc.
# Kristen Naegle
# January 18, 2008
sub returnLineHash($$){
    my $file = shift;
    my $keyCol = shift;

    my %hash;

    open(FX_98, $file) || die "Can't open $file file for reading \n";
    my $line = <FX_98>;
    while(defined($line = <FX_98>)){
	chomp $line;
	my @line = split("\t", $line);
	#my $key = returnTrypsPepKey($line[$keyCol]);
	my $key = $line[$keyCol];
	if(not $hash{$key}){
		$hash{$key} = [];
	    }
	    push @{$hash{$key}}, $line;
	}


    return \%hash;
    close(FX_98);

}

# $hash = returnAllLineHash($file, $keyCol)
# Returns a hash with key according to $keyCol and value is line of file
# Inputs: $file - tab separated pELM file
#         $keyCol - column number to use for key (i.e. pep:tryps)
# Outputs: $hash - hash reference, keys have replaced pY with y, etc.
# Kristen Naegle
# January 18, 2008
sub returnAllLineHash($$){
    my $file = shift;
    my $keyCol = shift;

    my %hash;

    my $uniqueOnly = 0; 

    open(FX_98, $file) || die "Can't open $file file for reading \n";
    my $line = <FX_98>;
    while(defined($line = <FX_98>)){
	chomp $line;
	my @line = split("\t", $line);
	#my $key = returnTrypsPepKey($line[$keyCol]);
#	my $key = $line[$keyCol];
	my $keys = returnAllFields($line, $keyCol, $uniqueOnly);
	foreach my $key (@$keys){
	    if(not $hash{$key}){
		$hash{$key} = [];
	    }
	    push @{$hash{$key}}, $line;
	}
    }


    return \%hash;
    close(FX_98);

}

# $uniqArray = uniq($array)
# Return a array of unique elements contained in input array
# Inputs: $array - reference to array
# Outputs: $uniqArray - reference to unique array
# Kristen Naegle
# December 18, 2007
sub uniq($){
    my $array = shift;
    my $hash = createUniqHash($array);
    my $uniq = hashToArr($hash);
    return $uniq;
}

# $outDir = checkDir($dir)
# appends / to the end of a directory if it isn't already there
# makes the directory if it doesn't exist
# Kristen Naegle
# February 4, 2008
sub checkDir($){
    my $dir = shift;
    my $outDir;
    if($dir !~ /\/$/){
	$outDir = $dir."/";
    }
    else{
	$outDir = $dir;
    }
    if(!-e $dir){
	my @dir = split("/", $dir);
	my $i = 1;
	my $dirName = "/".$dir[$i];	    
	print "Checking DIR: $dirName\n";
	if(!-e $dirName){ `mkdir $dirName`};
	for($i=2; $i <= $#dir; $i++){
	    $dirName = $dirName."/".$dir[$i];
	    print "Checking DIR: $dirName\n";
	    if(!-e $dirName){ `mkdir $dirName`};
	} 
    }

    return $outDir;

}

# \%dataHash = returnDataHash($header);
# Given a header, return a hash of datahash types. Datahash types have values 'type', 'label', 'priority'. 
# Inputs: $header - tab separated header from a data file
# Outputs: \%dataHash - keys are column numbers where data types exist and values are data hashes (i.e. hash of a hash)
# Kristen Naegle
# March 16, 2008
sub returnDataHash($){
    my ($header) = @_;
    my %dataHash;
    
    chomp $header;
    my $dataColArr = returnColumnNumberArrFromLine($header, "data");
    my @header = split("\t", $header);
    my $count = 1;
    if($dataColArr->[0] != -1){
	foreach my $colNum (@$dataColArr){
	    my %hash; 
	    my $col = $header[$colNum];
	    $col =~ s/[dD]ata://;
	    my @v = split(":", $col);
	    if(scalar(@v) > 2){
		handleError('returnDataHash', 'WARNING: The data columns have more than just data:type:label, ignoring remaining poritons', \@_);
	    }
	    $hash{'type'} = $v[0];
	    $hash{'label'} = $v[1];
	    if(!$hash{'label'}){
		handleError('returnDataHash', 'WARNING: A data columns is missing a label, should have data:type:label', \@_);
	    }
	    $hash{'priority'} = $count;
	    $count += 1;
	    $dataHash{$colNum} = \%hash;
	}
    }
    return \%dataHash;

}

# ($type, $label, $priority) = returnDataHashValues($dataHashRef)
# Returns the type, lable, and priority of a datahash ref
# Inputs: $dataHashRef - reference to data hash
# Outputs: $type - data type
#          $label - data label
#          $priority - data priority 
# Kristen Naegle
# March 17, 2008
sub returnDataHashValues($){
    my ($hashRef) = @_;
    my ($type, $label, $priority);

    $type = $hashRef->{'type'};
    $label = $hashRef->{'label'};
    $priority = $hashRef->{'priority'};

    return($type, $label, $priority);
    
}

# $run = returnRun($line, $runCol)
# Returns the run from a line (tab separated). Returns 'average' if runCol is -1
# Inputs: $line - tab separated line
#         $runCol - column number of run
# Outputs: $run - the run field of line
# Kristen Naegle
# March 18, 2008
sub returnRun($$){
    my ($line, $runCol) = @_;
    my $run;
    if($runCol == -1){
	$run = 'average';
    }
    else{
	$run = returnField($line, $runCol);
    }
    return $run;
}

# \@colNums = returnColumnNumberArrFromLine($header, $field)
# Returns an array of column numbers (zero-based) from header in a tab dilineated data file
# e.g. If you want all data labels then call returnColumnNumber($file, "data:");
#else if you want a specific label then call returnColumnNumber($file, "data:t0");
#inputs: $header - tab seperated header line of a file
#        $field - label descriptor that you are looking for
#outputs: $col - reference to an array of column numbers that matched that descriptor
# Kristen Naegle
# March 16, 2008
sub returnColumnNumberArrFromLine($$){
    my ($line, $field) = @_;
    $field = lc($field);
    chomp $line;
    my @line = split("\t", $line);
    my $i;
    my @col;
    #print "$field\n";
    $field =~ s/'"'//g; #have to strip quotes...can't regex compare
    for($i=0; $i <= $#line; $i++){
	#print "$line[$i]\n";
	#if(lc($line[$i]) =~ /lc($field)/){
	if(lc($line[$i]) =~ /$field/ or $field eq "all"){
	    push @col, $i;
	    #print "You are getting column $i with label $line[$i]\n";
	}
    }
    if(not @col){
	push @col, -1; #indicates error
	#print "ERROR: Could not find column..perhaps your file does not have a header? See fileTools.pm for header types\n";
    }
    return \@col;
}

1;

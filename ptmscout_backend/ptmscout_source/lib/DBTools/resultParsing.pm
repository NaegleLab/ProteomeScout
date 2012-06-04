use strict;
use warnings;
use DBI;

# \@arr = returnArrayOfResultsOnCol($sth, $col)
# Returns array of results based on a column number
# Inputs: $sth - statement handler
#         $col - column number
# Outputs: \@arr - reference to array of results
# Kristen Naegle
# March 5, 2008
sub returnArrayOfResultsOnCol($$){
    my ($sth, $col) = @_;
    my @arr;
    my @row = $sth->fetchrow_array;
    
    if(defined($row[$col])){
	push @arr, $row[$col];
	while(my @row = $sth->fetchrow_array){
	    push @arr, $row[$col];
	}
    }
    else{
	#push @arr, -1;
	$arr[0] = -1;
    }
    return \@arr;

}

# ($errorCode, $result) = returSingleResultOnCol($sth, $col)
# Returns array of results based on a column number
# Inputs: $sth - statement handler
#         $col - column number
# Outputs:$errorCode - 1 if there was more than one result 
#         $result - -1 if it doesn't exist, else the first value in array of results
# Kristen Naegle
# March 5, 2008 
sub returnSingleResultOnCol($$){
    my ($sth, $col) = @_;
    my $errorCode = 0;
    my $result;
    my @row = $sth->fetchrow_array;
    if(@row){
	$result = $row[$col];
	if(defined($sth->fetchrow_array)){
	    $errorCode = 1; #More than one result -- handle this outside the program since sth and col are not useful
	}
    }
    else{
	$result = -1;
    }
    return ($errorCode, $result);


}

# $hashRef = returnHashResult($sth, $keyCol, $valColArr)
# Parses sth return into a hash  
# Assumes that the key col are unique (i.e. returning a count(*) as val for a label which is keyCol - if not, then just pushes onto array 
# Inputs: $sth - executed statement handle
#         $keyCol - index into @row = $sth->fetchrow_array
#         $valCol - index into @row of value to be pushed onto array
# Outputs: $hashRef - ref to hash of arrays, key is value in Keycol and array is entries of @valCol
# Kristen Naegle
# July 1, 2008
# Modified July 20, 2008 to accept NULL value and report as NULL;
sub returnHashResult($$$){
    my($sth, $keyCol, $valColArr) = @_;
    my %hash;
    my $count = 0;
    while(my @row = $sth->fetchrow_array){
	my $key;
	if(defined $row[$keyCol]){
	    $key = $row[$keyCol];
	}
	else{
	    $key = 'NULL';
	}
#	print "@row\n";
	if(defined($hash{$key})){
	    handleError('returnHashResult', 'Function meant to handle unique key columns only! I.e. Return of SELECT grouped by a term', \@_);
	}
	#$hash{$key} = [];
	foreach my $valCol (@$valColArr){
	    #print "KEY: $key \t COUNT $row[$valCol]\n";
	    push @{$hash{$key}}, $row[$valCol];
	}
	
    }
    return \%hash;
}

1;

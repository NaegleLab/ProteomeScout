#! /usr/bin/perl
use strict;
use warnings;
#use globalVars;

# $timeStr = returnTimeStr;
# Returns date string in year-month-day hours:minutes:seconds
# Inputs: none
# Outputs: $str - date string
# January 22, 2008
sub returnTimeStr{
    my ($sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdist) = localtime(time);
    my $str = sprintf("%4d-%02d-%02d %02d:%02d:%02d\n", $year+1900, $mon+1, $mday, $hour, $min, $sec);
    return $str;
    
}

# ($wday,$hour) = returnDayTime;
# Returns the day and the hour, Sunday=0 Saturday=6;
# Inputs: none
# Outputs: $wday - 0-6 sunday-saturday
#          $hour - 24 hour time
# January 22, 2008
sub returnDayTime{
    my ($sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdist) = localtime(time);
    return ($wday, $hour);
    
}

# $timeStr = returnMonthYear
# Returns date string in Month-year (eg. 02-2008)
# Inputs: none
# Outputs: $str - Month-Year
# January 22, 2008
sub returnMonthYear{
    my ($sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdist) = localtime(time);
    my $str;
    $str = sprintf("%02d-%4d", $mon+1, $year+1900);
    return $str;
    
}

# $CPU = returnCPU();
# Returns the number of CPUs to use.. weekday worktimes is 6, other is 8. Work day from 8am to 4pm
# Outputs: $CPU - number of cpus to use
# Kristen Naegle
# March 24, 2008
# Mod June 9, 2010, now have CPU in globalVars.pm
sub returnCPU{
    my $CPU;
    my ($day, $hour) = returnDayTime();
    if($hour >= 16){
	#$CPU = 8;
	$CPU = $globalVars::CPU_NIGHT;
	
    }
    elsif($day == 6 || $day == 0){
	$CPU = $globalVars::CPU_NIGHT;
	#$CPU = 8;
    }
    else{
	#$CPU = 6;
	$CPU = $globalVars::CPU_DAY;
    }

#override = temp
    #$CPU = 10;
    return $CPU;
}

# $intersectArr = returnIntersection($array1, $array2)
# Inputs $arr1 - ref to array 1
#        $arr2 - ref to array 2
# Outputs: $arr - interesection of array 1 and array 2
# Kristen Naegle
# ?
sub returnArrayIntersection($$){
    my ($arr1, $arr2) = @_;
    my %hash;
    my @intersect;
    foreach my $item (@$arr1){
	if(not defined($hash{$item})){
	    $hash{$item} = 1;
	}
	
    }
    foreach my $i (@$arr2){
	if($hash{$i}){
	    push @intersect, $i;
	}
    }
    return \@intersect;
}

# $unqiueArr = returnUniqueArray($arr)
# Given a reference to an array, return the unique array
# Inputs: $arr - ref to an array of objects
# Outputs: $unique - ref to an array of the uniquified objects
# Kristen Naegle
# January 28, 2009
sub returnUniqueArray($){
    my ($arr) = @_;
    my %hash;
    foreach my $item (@$arr){
	if(not defined($hash{$item})){
	    $hash{$item} = 0;
	}
	$hash{$item} += 1;
    }
	
    my @unique = keys(%hash);
    return \@unique;
}

# $str = joinArrSingleQuotes($arr)
# Add single quotes around array elements and return as string (with space between elements as well);
# Inputs: $arr - ref to an array of elements
# Outputs: $str - the string of joined array elements surrounded in quotes
# Kristen Naegle
# January 28, 2009
sub joinArrSingleQuotes($){
    my ($arr) = @_;
    my $str;
    foreach my $item (@$arr){
	$str .= "'".$item."' ";

    }
    return $str;
}

# $arrString = makeSelectINString($arr, $QUOTES)
# Create a select var string from an array i.e. if $arr = 1, 2,3 output is (1, 2, 3) if QUOTES then output is ('1', '2', '3');
# Inputs: $arr - ref. to array of elements to put into string
#         $QUOTES - if 0, no quotes, if 1 QUOTES
# OUtputs: $string - string suitable for using in SELECT .. IN
# Kristen Naegle
# May 12, 2009
sub makeSelectINString($$){
    my ($arr, $QUOTES) = @_;
    my $arrStr = '(';
    for(my $i=0; $i<scalar(@$arr) - 1; $i++){
	$arr->[$i] =~ s/\'//g;
	if($QUOTES){
	    $arrStr .= "'$arr->[$i]',";
	    
	}
	else{
	    $arrStr .= $arr->[$i].",";
	}

    }
    $arr->[scalar(@$arr)-1] =~ s/\'//g;
    if($QUOTES){
	$arrStr .= "'".$arr->[scalar(@$arr)-1]."')";
	
    }
    else{
	$arrStr .= $arr->[scalar(@$arr)-1].")";
    }
    return $arrStr;
}


# $hashRef = arrToHash($arrRef);
# takes an array and creates a unique hash with keys in array and value is the number of times it occurs in array
# Inputs: $arrRef - referenc to an array
# Outputs: $hashRef - reference to a hash with unique keys in array and value is number of occurances
# Kristen Naegle
# November 1, 2007
sub arrToHash($){

    my $arrRef = shift;
    my %hash;
    foreach my $item (@$arrRef){
	if(not $hash{$item}){
	    $hash{$item} = 0;

	}
	$hash{$item} += 1;

    }
    return \%hash;
}

# ($RELEASE, $file) = returnLinkRelease($link);
# Given a symbolic link, return the release version (for example swisspfam points to swisspfam_relX return X
# Inputs: $link - symbolic link
# Outputs:$file - file symbolic link points to 
#         $X - Version 
# Kristen Naegle
# June 18, 2009
sub returnLinkRelease($){
    my ($link) = @_;
    my $file = readlink($link);
    $file =~ m/_rel(.+)/;
    my $VER = $1;

    return($VER, $file);
    
}

# $hashRef = returnHashAllowedSpecies()
# Returns a hash with keys equal to the binomial names of species allowed in PTMScout
# Outputs: $hashRef - see above, ref to described hash, values are taxonomy
# Kristen Naegle
# Sept. 23, 2009
sub returnHashAllowedSpecies(){
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
    $hash{'rattus norvegicus'} = 10166;
    $hash{'saccharomyces cerevisiae'} = 4932;
    $hash{'schizosaccharomyces pombe'} = 4896;
    $hash{'xenopus laevis'} = 8355;
    $hash{'zea mays'} = 4577;
    $hash{'gallus gallus'} = 9031;
    $hash{'bacillus subtilis'} = 1423;
    $hash{'oryctolagus cuniculus'} = 9986;
    $hash{'sus scrofa'} = 9823;
    $hash{'canis lupus'} = 9612;

    return \%hash;

}

1;

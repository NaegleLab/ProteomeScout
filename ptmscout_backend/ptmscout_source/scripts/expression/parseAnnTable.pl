#! /usr/bin/perl
use strict;
use warnings;
use fileTools;

#problem - annotation table does not have a line entry for each entry.  Many redundancies..so need to parse entries and copy them. 

### NEVER MIND - This is not the problem!!

my $dataFile = "gnf1b-anntable.txt"; #downloaded Jan 15, 2009 
my $outputFile = "gnf1b_anntable_parsed.txt";  # from here still need to go on and delete columns to create master table.

my $colNums_probe = returnColumnNumberArr($dataFile, 'Probeset ID');
my $colNums_reporters = returnColumnNumberArr($dataFile, 'Reporters');

my $probeCol = $colNums_probe->[0];
my $reporterCol = $colNums_reporters->[0];

print "FOUND probe: $colNums_probe->[0]\n";
print "FOUND reporters: $colNums_reporters->[0]\n";

open(DATA, $dataFile) || die "Can't open $dataFile for reading\n";
open(OUT, ">$outputFile") || die "Can't open $outputFile for writing\n";

my $line = <DATA>;
print OUT $line;

while(defined($line = <DATA>)){ 
#foreach line if there are multiple reporters reprint line 
    my %hash; 
    chomp $line;
    my @line = split("\t", $line);
    my $probe = $line[$probeCol];
    my $reporters = $line[$reporterCol];
    my $rs = returnReporters($reporters);
    print OUT $line."\n"; #print line
    $hash{$probe} = 1; 

    #print "$probe has $reporters\n";
    foreach my $r (@$rs){
	#print "\t$r";
	if(not defined $hash{$r}){
	    $hash{$r} = 1;
	    $line[$probeCol] = $r;
	    foreach my $item (@line){
		print OUT $item."\t";
	    }
	    print OUT "\n";
	}
    }
#    print "\n";


    


}

close(DATA);
close(OUT);


sub returnReporters($){
    my ($reporters) = @_;
    my @reporters = split(';', $reporters);
    my @parsed;
    foreach my $r (@reporters){
	my @rs = split(" ", $r);
	my $rstrip = $rs[0];
	$rstrip =~ s/ //g;
	push @parsed, $rstrip;

    }
    return \@parsed;
}

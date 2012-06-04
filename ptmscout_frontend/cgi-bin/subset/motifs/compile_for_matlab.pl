#! /usr/bin/perl -w
use strict;

my $in = "sites_final.csv";
my $out = "for_matlab.csv";

open(IN, $in) or die "Can't read $in";
my $line = <IN>;

my %data;

while(defined(my $line = <IN>)){
    chomp $line;
    $line =~ s/\"//g;
    my @line = split('\t', $line);
    my $id = $line[0];
    if (! defined $data{$id}){
	$data{$id} = ();
	$data{$id}->[0] = "$line[2] $line[7]$line[6]";
	$data{$id}->[1] = $line[8];
	$data{$id}->[2] = $line[9];
	$data{$id}->[3] = $line[10];
	$data{$id}->[4] = $line[11];
    }
    else { 
	$data{$id}->[0] .= ",$line[7]$line[6]";
    }
}
open(OUT, ">$out") or die "Can't write $out";
foreach my $key (sort keys %data) {
    print OUT join("\t", @{$data{$key}}) . "\n";
}

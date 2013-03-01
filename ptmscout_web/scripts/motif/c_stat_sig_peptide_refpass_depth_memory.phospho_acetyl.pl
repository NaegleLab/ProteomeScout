#! /usr/bin/perl -w
use strict;
#use Math::Counting ':long';
use POSIX;
use Getopt::Long;
$| = 1;

my $foreground;
my $background;
my $pattern;
my $cutoff;
my $numpept;
my $binomial = 0;

GetOptions('foreground=s' => \$foreground,
	   'background=s' => \$background,
	   'pattern:s'    => \$pattern,
	   'cutoff:f'     => \$cutoff,
	   'numpept:i'    => \$numpept,
	   'binomial'     => \$binomial); 

if (! -e $foreground || ! -e $background) {
    print "Usage: ./c_stat_sig_protein_singles.pl --foreground FILE --background FILE [--pattern PATTERN] [--cutoff CUTOFF] [--numpept NUMPEPT] [--binomial]\n";
    print "       PATTERN = pos1.aa1.pos2.aa2...\n";
    print "       CUTOFF  = p-value required to forward to recursion.\n";
    print "       NUMPEPT = Number of peptides required to forward to recursion.\n";
    print "       If the binomial flag is used, a binomial distribution is used to calculate \n";
    print "        significance.  By default, hypergeometric is used.\n";
    print "       Work directory requires a file named conversion that includes amino acid list.\n";    print "        See ./sample_conversion\n";
    exit;
}

if (! defined $cutoff) {$cutoff = 1e-2;}
if (! defined $numpept) {$numpept = 0;}

# Read in the foreground and background.
my @foreground = &read_peptides($foreground);
my @background = &read_peptides($background);

my $froot = $foreground;
my $broot = $background;
my @froot = split('/', $froot);
my @broot = split('/', $broot);
$froot = $froot[$#froot];
$broot = $broot[$#broot];
my $log = "log.$froot.$broot";
open(LOG, ">$log") or die "Can't write $log";
my $oldfh = select(LOG);  $| = 1;  select($oldfh);

my %done;
my %fdone;
my @for_singles;
my @for_doubles;
my @for_next_doubles;
my %fg_by_pattern;
my %bg_by_pattern;
my $pos = 14;

# Read in the conversion file.
my @aa;  my %aa_multi;  my %multi_aa;  my %aa_pattern;
&read_conversion(\@aa, \%aa_multi, \%multi_aa, \%aa_pattern);

my $Nf = scalar(@foreground);
my $Nb = scalar(@background);

my $fail_done = 0;
my $fail_sub = 0;
my $fail_count = 0;
my $fail_sig = 0;
my $fail_bestsig = 0;
my $fail_ancestry = 0;
my $fail_checkpat = 0;
my $num_tests_for_kristen = 0;
my $succeed = 0;

# I need to recurse the fg and bg to eliminate failed matches.
my %checked_pat;
# KMN Changes 3/1/2013 -- Parsing the foreground and initiating the pattern based on the presence of any amino acid in the central position occuring more than twice
my %patHash;
foreach my $pep (@foreground){
    my $centralResidue = substr($pep, 7, 1);
    if(not defined $patHash{$centralResidue}){
	$patHash{$centralResidue} = 0;
    }
    $patHash{$centralResidue} += 1;

}
foreach my $res (keys %patHash){
    if($patHash{$res} >= 2){
	my $str = ".......".$res.".......";
	print "DEBUG: adding pattern in $str\n";
	push @for_singles, $str;
    push @for_doubles, $str;
    &init_refhash($str, \@foreground, \%fg_by_pattern);
    &init_refhash($str, \@background, \%bg_by_pattern);
	
    }
}


# if (! defined $pattern){ 
#     push @for_singles, ".......y.......";
#     &init_refhash(".......y.......", \@foreground, \%fg_by_pattern);
#     &init_refhash(".......y.......", \@background, \%bg_by_pattern);
#     push @for_singles, ".......s.......";
#     &init_refhash(".......s.......", \@foreground, \%fg_by_pattern);
#     &init_refhash(".......s.......", \@background, \%bg_by_pattern);
#     push @for_singles, ".......t.......";
#     &init_refhash(".......t.......", \@foreground, \%fg_by_pattern);
#     &init_refhash(".......t.......", \@background, \%bg_by_pattern);
#     push @for_singles, ".......x.......";
#     &init_refhash(".......x.......", \@foreground, \%fg_by_pattern);
#     &init_refhash(".......x.......", \@background, \%bg_by_pattern);
#     push @for_singles, ".......k.......";
#     &init_refhash(".......k.......", \@foreground, \%fg_by_pattern);
#     &init_refhash(".......k.......", \@background, \%bg_by_pattern);
#     push @for_doubles, ".......y.......";
#     push @for_doubles, ".......s.......";
#     push @for_doubles, ".......t.......";
#     push @for_doubles, ".......x.......";
#     push @for_doubles, ".......k.......";
# }
# else {
#     push @for_singles, $pattern;
#     push @for_doubles, $pattern;
#     &init_refhash($pattern, \@foreground, \%fg_by_pattern);
#     &init_refhash($pattern, \@background, \%bg_by_pattern);

# }

undef @foreground;  undef @background;

print join(" ", sort @for_singles) . "\n";

print "Initializing with " . scalar(@for_singles) . " patterns.\n";
# Run on the initialized singles and doubles.
foreach my $pat (sort @for_singles) {
    &singles($pat, $Nf, $Nb);
    &search_paironce($pat, $Nf, $Nb);
}


close LOG;
print "\n";
print "FAIL_DONE $fail_done\n";
print "FAIL_SUB $fail_sub\n";
print "FAIL_COUNT $fail_count\n";
print "FAIL_SIG $fail_sig\n";
print "FAIL_BESTSIG $fail_bestsig\n";
print "FAIL_ANCESTRY $fail_ancestry\n";
print "FAIL_CHECKPAT $fail_checkpat\n";
print "SUCCEED $succeed\n";
print "NUM_TESTS_FOR_KRISTEN $num_tests_for_kristen\n";
print "MEM: " . `ps -o rss,vsz $$`;
exit;

################################################################
# Given a source array of scalar peptides and a pattern, add references
# to the peptides that match the pattern to the destination hash.
sub init_refhash {
    my $pattern = shift;
    my $source = shift;
    my $dest = shift;
    my $regex = &make_pattern($pattern, \%aa_pattern);
    $dest->{$pattern} = ();
    
    foreach my $peptide (@{$source}){
	if ($peptide =~ /$regex/){
	    push @{$dest->{$pattern}}, \$peptide;
	}
    }
}

################################################################
# Dive depth-first through singles.  Clear doubles on the way back out.
sub singles {
    my $pattern = shift;
    my $Nf = shift;
    my $Nb = shift;
    my @pattern = split('', $pattern);

    my $parent = $pattern;

    # Iterate over positions and amino acids.
    foreach my $p (0..$pos){
	if ($pattern[$p] ne ".") {next;} 
	foreach my $aa (@aa){
	    $pattern[$p] = $aa;
	    $pattern = join("", @pattern);

	    print "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b";
	    print $pattern . " singles";

	    if (defined($checked_pat{$pattern})) {
		$pattern[$p] = ".";
		$pattern = join("", @pattern);
		$fail_checkpat++;
		next; 
	    }	
	    $checked_pat{$pattern} = 1;

	    # Skip this pattern if it's been done before successfully..
	    if (defined($done{$pattern})) { 
		$pattern[$p] = ".";
		$pattern = join("", @pattern);
		$fail_done++;
		next; 
	    }

	    my $m = $Nf;

	    # Count pattern in foreground, making an entry in the hash of refs.
	    my $c_xj = &count_occurences($pattern, $parent, \%fg_by_pattern);
	    
	    # Skip this pattern if too few peptides contain it in the foreground.
	    if ($c_xj < $numpept) {
		delete $fg_by_pattern{$pattern};
		$pattern[$p] = ".";
		$pattern = join("", @pattern);
		$fail_count++;
		next; 
	    }

	    # If any parent of this pattern other than the one we are here from 
	    # is successfully done, this pattern's been checked (and failed).  Skip it.
	    if (&check_ancestry($pattern, $parent)) {
		$pattern[$p] = ".";
		$pattern = join("", @pattern);
		$fail_ancestry++;
		next;
	    }

	    # Test a lower bound on score.  If it fails, skip it.
	    if (&hypergeometric2($Nf, $Nb, $c_xj, $c_xj) > $cutoff) {
		delete $fg_by_pattern{$pattern};
		$pattern[$p] = ".";
		$pattern = join("", @pattern);
		$fail_bestsig++;
		next;
	    }

	    # Test if this motif is a submotif of a previously done motif, with
	    # the same membership.
	    if (&is_a_submotif_with_same_size($pattern, \%done, $c_xj, \%aa_multi)) { 
		delete $fg_by_pattern{$pattern};
		$pattern[$p] = ".";
		$pattern = join("", @pattern);
		$fail_sub++;
		next; 
	    }
	    
	    # Only if all of the above checks out is it acceptable to count the background.
	    # THIS IS THE SLOW STEP!
	    my $D = &count_occurences($pattern, $parent, \%bg_by_pattern);
	    my $p_xj = $D / $Nb;
	    my $sig;
	    
	    if ($binomial == 1){
		$sig = &binomial($m, $c_xj, $p_xj);
	    }
	    else {
		$num_tests_for_kristen++;
		$sig = &hypergeometric2($Nf, $Nb, $c_xj, $D);
	    }
	    
	    # If the significance is there, do all singles by depth, and pairs on the way back up.
	    if ($sig <= $cutoff) {
		$done{$pattern} = $c_xj;
		#push @for_doubles, $pattern;
		printf LOG "| $pattern | %5d / %5d | %5d / %5d | %4.2E |\n", $c_xj, $Nf, $D, $Nb, $sig;
		&singles("$pattern", $Nf, $Nb);
		&search_paironce($pattern);
		delete $fg_by_pattern{$pattern};
		delete $bg_by_pattern{$pattern};
		$succeed++;
	    }
	    else {
		delete $fg_by_pattern{$pattern};
		delete $bg_by_pattern{$pattern};
		$fail_sig++;
	    }
	    $pattern[$p] = ".";
	    $pattern = join("", @pattern);
	}
    }
}

##############################################################
# Given a pattern, test all doubles off of it, as they are found by going
# through the foreground.
sub search_paironce {
    my $pattern = shift;
    my $code = &pattern_to_regex($pattern);
    
    foreach my $peptide (@{$fg_by_pattern{$pattern}}){
	if ($$peptide =~ /^$code$/) {
	    &pairs_informed($pattern, $$peptide);
	}
    }

    delete $fg_by_pattern{$pattern};
    delete $bg_by_pattern{$pattern};

}

##############################################################
# Given a pattern and a peptide, test all pair-descendant patterns in that peptide.
sub pairs_informed {
    my $pattern = shift;
    my $peptide = shift;
    my @pattern = split("", $pattern);    
    my @peptide = split("", $peptide);

    my $parent = $pattern;
    foreach my $p (0..$pos-1){
	if($pattern[$p] ne ".") {next;}
	foreach my $p2 ($p+1..$pos){
	    if($pattern[$p2] ne ".") {next;}
	    foreach my $aa (keys %{$multi_aa{$peptide[$p]}}){
		$pattern[$p] = $aa;
		$pattern = join("", @pattern);
		foreach my $aa2 (keys %{$multi_aa{$peptide[$p2]}}){
		    $pattern[$p2] = $aa2;
		    $pattern = join("", @pattern);
		    print "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b";
		    print $pattern . " doubles";


		    if (defined($checked_pat{$pattern})) {
			$pattern[$p2] = ".";
			$pattern = join("", @pattern);
			$fail_checkpat++;
			next; 
		    }	
		    $checked_pat{$pattern} = 1;

		    if (defined($done{"$pattern"})) { 
			$pattern[$p2] = ".";
			$pattern = join("", @pattern);
			$fail_done++;
			next; 
		    }
		    my $m = $Nf;

		    # Count pattern in foreground, making an entry in the hash of refs.
		    my $c_xj = &count_occurences($pattern, $parent, \%fg_by_pattern);
	    
		    # Skip this pattern if too few peptides contain it in the foreground.
	    	    if ($c_xj < $numpept) {
			delete $fg_by_pattern{$pattern};
			$pattern[$p2] = ".";
			$pattern = join("", @pattern);
			$fail_count++;
			next; 
		    }

		    # If any parent of this pattern other than the one we are here from 
		    # is successfully done, this pattern's been checked (and failed).  Skip it.
		    if (&check_ancestry($pattern, $parent)) {
			$pattern[$p2] = ".";
			$pattern = join("", @pattern);
			$fail_ancestry++;
			next;
		    }

	    	    # Test a lower bound on score.  If it fails, skip it.
		    if (&hypergeometric2($Nf, $Nb, $c_xj, $c_xj) > $cutoff) {
			delete $fg_by_pattern{$pattern};
			$pattern[$p2] = ".";
			$pattern = join("", @pattern);
			$fail_bestsig++;
			next;
		    }
	    
		    # Test if this motif is a submotif of a previously done motif, with
		    # the same membership.
		    if (&is_a_submotif_with_same_size($pattern, \%done, $c_xj, \%aa_multi)) { 
			delete $fg_by_pattern{$pattern};
			$pattern[$p2] = ".";
			$pattern = join("", @pattern);
			$fail_sub++;
			next; 
		    }

		    # Only if all of the above checks out is it acceptable to count the background.
		    # THIS IS THE SLOW STEP!
		    my $D = &count_occurences($pattern, $parent, \%bg_by_pattern);
		    my $p_xj = $D / $Nb;
		    my $sig;
	    
		    if ($binomial == 1){
			$sig = &binomial($m, $c_xj, $p_xj);
		    }
		    else {
			$num_tests_for_kristen++;
			$sig = &hypergeometric2($Nf, $Nb, $c_xj, $D);
		    }
		    
		    # If the significance is there, do all singles by depth, and pairs on the way back up.
		    if ($sig <= $cutoff) {
			$done{"$pattern"} = $c_xj;
			push @for_singles, $pattern;
			push @for_next_doubles, $pattern;
			&singles("$pattern", $Nf, $Nb);
			&search_paironce($pattern);
			delete $fg_by_pattern{$pattern};
			delete $bg_by_pattern{$pattern};
			printf LOG "| $pattern | %5d / %5d | %5d / %5d | %4.2E |\n", $c_xj, $Nf, $D, $Nb, $sig;
			$succeed++;
		    }
		    else {
			delete $fg_by_pattern{$pattern};
			delete $bg_by_pattern{$pattern};
			$fail_sig++;
		    }
		    $pattern[$p2] = ".";
		    $pattern = join("", @pattern);
		}
		$pattern[$p] = ".";
		$pattern = join("", @pattern);

	    }
	}
    }
}
##############################################################
# Tests whether the 'submotif' is an ancestor of any known motif, and the same size.
sub is_a_submotif_with_same_size(){
    my $submotif = shift;
    my $motifs = shift;
    my $size = shift;
    my $multis = shift;
    foreach my $motif (keys %{$motifs}){
	if ($size != $motifs->{$motif}){ next; }
	my $result = 1;
	my @submotif = split('', $submotif);
	my @motif = split('', $motif);
	foreach my $i (0..$#motif){
	    if ($submotif[$i] eq ".") { next;}
	    else {
		if ($motif[$i] eq $submotif[$i]) { next;}
		else {
		    if (defined($multis->{$submotif[$i]}->{$motif[$i]})){ next; }
		    else{ $result = 0; last; }
		}
	    }
	}
	if ($result) { return 1;}
    }
    return 0;
}
##############################################################
# Counts the occurences of a motif in a hash of references, just from among a parent motif.
sub count_occurences(){
    my $pattern = shift;
    my $parent = shift;
    my $hash = shift;
    
    my $regex = &make_pattern($pattern, \%aa_pattern);
    # Find and count $pattern in $hash.
    my $count = 0;

    $hash->{$pattern} = ();
    
    foreach my $peptide (@{$hash->{$parent}}) {
	if ($$peptide =~ /$regex/) {
	    $count++;
	    push @{$hash->{$pattern}}, $peptide;
	}
    }
    return $count;
}

##############################################################
# Turn a human pattern into a regular expression.
sub make_pattern(){
    my $scalar = shift;
    my $aa = shift;

    my @p_arr = split('', $scalar);
    my $pattern = "";
    foreach my $p (@p_arr){
	if (defined($aa->{$p})){
	    $pattern .= $aa->{$p};
	}
	else { $pattern .= "."; }
    }
    return $pattern;
}
##############################################################
# Approximate a statistical significance using the binomial distribution.
sub binomial(){
    # Arguments: (here $m and the denominator of $p_xj are constant)
    my $m = shift;     # Size of selected data.
    my $c_xj = shift;  # Number of hits in selected data.
    my $p_xj = shift;  # Probability of hit in unselected data;
    my $sum = 0;

    foreach my $i ($c_xj..$m){
        $sum += (combination($m, $i)) * ($p_xj**$i) * (1-$p_xj)**($m-$i);
    }
    return $sum;
}
##############################################################
# Calculate a statistical significance using the hypergeometric distribution.
sub hypergeometric(){
    my $n = shift;     # Size of selected data.
    my $N = shift;     # Size of unselected data.
    my $k = shift;     # Number of hits in selected data.
    my $D = shift;     # Number of hits in unselected data.
    my $sum = 0;
    
    if ( ($k == 0) ) { return 1; }
    if ( ($n == $k) && ($N == $D) ) { return 1; }
    my $max = $n;
    if ($D < $n) { $max = $D; }
    foreach my $i ($k..$max){
	my $log = 0;
	foreach my $index (0..$i-1){
	    $log += (log10($D-$index) - log10($i-$index));
	}
	foreach my $index (0..$n-$i-1) {
	    $log += (log10($N-$D-$index) - log10($n-$i-$index));
	}
	foreach my $index (0..$n-1){
	    $log -= (log10($N-$index) - log10($n-$index));
	}
	$sum += 10**$log;
    }

    return $sum;
}

##############################################################
# Calculate a statistical significance using a gamma-function
# approximation to the hypergeometric that's fast and accurate to like
# 8 digits
sub hypergeometric2(){
    my $n = shift;     # Size of selected data.
    my $N = shift;     # Size of unselected data.
    my $k = shift;     # Number of hits in selected data.
    my $D = shift;     # Number of hits in unselected data.
    my $sum = 0;
    
    if ( ($k == 0) ) { return 1; }
    if ( ($n == $k) && ($N == $D) ) { return 1; }
    my $max = $n;
    if ($D < $n) { $max = $D; }
    foreach my $i ($k..$max){
	my $log = 0;
	$log += &factln($D) - &factln($i) - &factln($D-$i);
	$log += &factln($N-$D) - &factln($n-$i) - &factln($N-$D-$n+$i);
	$log -= &factln($N) - &factln($n) - &factln($N-$n);
	
	$sum += exp($log);
    }

    return $sum;
}
##############################################################
# Gamma function-related math.
sub factln {
    my $x = (shift) + 1;
    my $tmp = $x + 5.5;
    $tmp -= ($x + .5) * log($tmp);
    my $ser = 1.000000000190015
	    + 76.18009172947146    / ++$x
	    - 86.50532032941677    / ++$x
	    + 24.01409824083091    / ++$x
	    -  1.231739572450155   / ++$x
	    +  0.12086509738661e-2 / ++$x
	    -  0.5395239384953e-5  / ++$x;
    return log(2.5066282746310005*$ser/($x-6)) - $tmp;
}
##############################################################
sub read_peptides {
    my $input = shift;
    my @array;
    open(IN, $input) or die "Can't read $input";
    while(defined(my $line = <IN>)){
	chomp $line;
	push @array, $line;
    }
    close IN;
    return @array;
}
##############################################################
sub read_conversion {
    my $aa = shift;
    my $aa_multi = shift;
    my $multi_aa = shift;
    my $aa_pattern = shift;

    open(IN, "conversion") or die "Can't read conversion";
    while(defined(my $line = <IN>)){
	my @line = split(' ', $line);
	push @{$aa}, $line[0];
	$aa_pattern->{$line[0]} = $line[1];
	if ($line[1] !~ /\[/){
	    if (!defined $multi_aa->{$line[1]}){
		$multi_aa->{$line[1]} = {};
	    }
	    $multi_aa->{$line[1]}->{$line[0]} = 1;
	    if (!defined $aa_multi->{$line[0]}){
		$aa_multi->{$line[0]} = {};
	    }
	    $aa_multi->{$line[0]}->{$line[1]} = 1;
	}
	else {
	    my @multi = split("", $line[1]);
	    foreach my $i (1..$#multi-1){
		if (!defined $multi_aa->{$multi[$i]}){
		    $multi_aa->{$multi[$i]} = {};
		}
		$multi_aa->{$multi[$i]}->{$line[0]} = 1;
		if (!defined $aa_multi->{$line[0]}){
		    $aa_multi->{$line[0]} = {};
		}
		$aa_multi->{$line[0]}->{$multi[$i]} = 1;
	    }
	}
    }
    close IN;
}
##############################################################
sub pattern_to_regex(){
    my $pattern = shift;
    my @pattern = split('', $pattern);
    foreach my $i (0..$#pattern){
	if ($pattern[$i] eq ".") {next;}
	$pattern[$i] = $aa_pattern{$pattern[$i]};
    }
    return (join("", @pattern));
}
##############################################################
sub check_ancestry(){
    my $pattern = shift;
    my $parent = shift;
#    my $grandparent = shift;
    my @pattern = split('', $pattern);
    my $good_parent = 0;
    foreach my $p1 (0..$#pattern) {
	if ($p1 == 7) {next;}
	if ($pattern[$p1] eq "."){ next; }
	my $p1_orig = $pattern[$p1];
	$pattern[$p1] = ".";
	my $newp = join("", @pattern);
	if (defined($done{$newp})) {
	    if (($newp ne $parent)) {
		$good_parent = 1; last;
	    }
	}

	$pattern[$p1] = $p1_orig;
	
    }
    return $good_parent;
    
}
##############################################################

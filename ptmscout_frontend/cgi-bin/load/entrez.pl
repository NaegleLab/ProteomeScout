use entrezTools;

my ($acc,) = @ARGV;
my $type=returnAccType($acc);
print $type;
closeLOG;

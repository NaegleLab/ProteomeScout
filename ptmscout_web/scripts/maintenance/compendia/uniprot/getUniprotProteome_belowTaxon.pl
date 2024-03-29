use strict;
use warnings;
use LWP::UserAgent;
use HTTP::Date;

# Taxonomy identifier of top node for query, e.g. 2 for Bacteria, 2157 for Archea, etc.
# # (see http://www.uniprot.org/taxonomy)
if(scalar(@ARGV)<1){
    print "USAGE: getUniprotProteome_belowTaxon.pl <TAXON_NODE>\n";
    exit();

}
my $top_node = $ARGV[0];
#
my $reference = 0; # Toggle this to 1 if you want reference instead of complete proteomes.
my $proteome = $reference ? 'reference:yes' : 'complete:yes';
my $keyword = $reference ? 'keyword:1185' : 'keyword:181';
#
my $contact = 'proteomescout@seas.wustl.edu'; # Please set your email address here to help us debug in case of problems.
my $agent = LWP::UserAgent->new(agent => "libwww-perl $contact");

# Get a list of all taxons below the top node with a complete/reference proteome.
my $query_list = "http://www.uniprot.org/taxonomy/?query=ancestor:$top_node+$proteome&format=list";
my $response_list = $agent->get($query_list);
die 'Failed, got ' . $response_list->status_line .
    ' for ' . $response_list->request->uri . "\n"
     unless $response_list->is_success;

# For each taxon, mirror its proteome in FASTA format.
for my $taxon (split(/\n/, $response_list->content)) {
    #my $file = $taxon . '.fasta';
    my $dir = './data/'.$top_node;
    `touch $dir`;
    my $file = $dir.'/'.$taxon. '.txt';
    print $file;
    my $text = 'mod%20res';
    #my $query_taxon = "http://www.uniprot.org/uniprot/?query=$text+&organism:$taxon+$keyword&format=text&include=yes";
    my $query_taxon = "http://www.uniprot.org/uniprot/?query=$text+organism:$taxon+$keyword&format=txt&include=yes";
    my $response_taxon = $agent->mirror($query_taxon, $file);

    if ($response_taxon->is_success) {
        my $results = $response_taxon->header('X-Total-Results');
        my $release = $response_taxon->header('X-UniProt-Release');
        my $date = sprintf("%4d-%02d-%02d", HTTP::Date::parse_date($response_taxon->header('Last-Modified')));
        print "File $file: downloaded $results entries of UniProt release $release ($date)\n";
    }
    elsif ($response_taxon->code == HTTP::Status::RC_NOT_MODIFIED) {
        print "File $file: up-to-date\n";
    }
    else {
    die 'Failed, got ' . $response_taxon->status_line .
        ' for ' . $response_taxon->request->uri . "\n";
     }
}

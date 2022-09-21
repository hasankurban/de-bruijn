## Read in the blast output file from blasting oligos against all transcripts
## Output a list of hits for each oligo


use strict;
use warnings;

use Bio::SearchIO;


my %hit_hash = ();
my %hit_details = ();


open(OP1, ">oligo_hyb_list.tab");
open(OP2, ">oligo_hyb_count.tab");


my $infile = $ARGV[0];

my $report = new Bio::SearchIO(
	         -file=>"$infile",
	      -format => "blast"); 
    
while(my $result = $report->next_result) {

    my $query_accession = $result->query_accession; 
    
    my $hit_count = 0;

    print OP1 "$query_accession";
    
    print OP2 "$query_accession";

    while(my $hit = $result->next_hit) { 
	my $hit_name = $hit->name;
	my $hit_matches = $hit->matches('id');
	if($hit_matches  >= 20){
	    print "$hit_matches\n";
	    print OP1 "\t$hit_name";
	    $hit_count++;
	}
    }
    
    print OP1 "\n";
    print OP2 "\t$hit_count\n";

}

close OP1;
close OP2;

exit;

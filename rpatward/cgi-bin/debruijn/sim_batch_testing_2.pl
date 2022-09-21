#!/usr/bin/perl -w 

# created from sim_batch_testing_2.pl on 02/11/2006
# adding the alignment step

## Testing my similarity idea

### Usage:

# perl program_name -fam <input file name> -sub <substitution matrix file> 
# The input file should have sequences in FASTA format

use Getopt::Long;
use strict;
use Benchmark;

$|=1;

my $start_time = new Benchmark;


my ($i,$j,$m,$k,$alpha,$temp,$max,$max2,$max_main,$start,$key,$key1,$key2,$seq_len,$sequence,$maxkey1,$maxkey2,$motif,$path_weight, $reg_ex, $add, $break_flag);
my ($max_inner,$max_inner_key);
my ($subs_matrix, $protein_file, $seq_file);
my (@config_data,@str1,@str2,@seq,@len,%edge_dup,%edge_rec);
my ($node_len,$overlap,$offset,$threshold,$max_threshold, $graph_threshold);
my ($scale_threshold,$scale_max_threshold, $scale_graph_threshold,$nsites);
my ($score_messy, $score_per, $score_neg, $messy_match, $score);
my ($look_ahead_flag, $jump_count, $reg_ex_all, $reg_ex_thres, $reg_ex_look_ahead, $motif_look_ahead,$jump);
my ($reg_ex_temp, $reg_ex_middle, $option_count,@temp_maxkey1,$exact,$strict);
my ($triplet);


my ($sim_threshold, $sim_scaling_factor);
my %sim = ();



GetOptions ("sub=s" => \$subs_matrix , "fam=s" => \$protein_file,  "seq=s" => \$seq_file,  "nlen=s" => \$node_len, "ov=s" =>\$overlap, "th=s" =>\$scale_threshold,  "mth=s" =>\$scale_max_threshold, "gth=s" =>\$scale_graph_threshold, "n=s" =>\$nsites, "jump=s" => \$jump,  "triplet=s" => \$triplet);
   
# "sub" needs the name of substitution matrix file as the argument
# "fam" needs the name of the input file containing the protein family sequences in FASTA format as the argument
# "seq" needs the name of the file containing the protein sequence to be tested


#### Set the default values

if (!(defined($node_len))) {$node_len = 4;}
if (!(defined($scale_threshold))) {$scale_threshold = 70;}
if (!(defined($scale_max_threshold))) {$scale_max_threshold = 75;}
if (!(defined($scale_graph_threshold))) {$scale_graph_threshold = 100;}
if (!(defined($nsites))) {$nsites = 100;}
if (!(defined($jump))){$jump = 10;}
if (!(defined($sim_threshold))){$sim_threshold = 0;}
if (!(defined($sim_scaling_factor))){$sim_scaling_factor = 0.5; }


$overlap = $node_len - 1;
$offset = $node_len - $overlap;


####
 
my @matrix_file_data  = &open_file($subs_matrix);
my @protein_file_data = &open_file($protein_file);
my %matrix            = &read_matrix_file(@matrix_file_data);
my @string            = &process_file_data(@protein_file_data);

####

my $total_len = 0;

for($i=0;$i<=$#string;$i++){
    $len[$i] = length $string[$i];
    $total_len += $len[$i]; 
}

### create similarity measures
foreach my $amino_acid_1 (keys %matrix){
    foreach my $amino_acid_2 (keys %{$matrix{$amino_acid_1}}){
	if($matrix{$amino_acid_1}{$amino_acid_2} > $sim_threshold){
	    $sim{$amino_acid_1}{$amino_acid_2}=$matrix{$amino_acid_1}{$amino_acid_2};
	}
    }
}



my $number_of_seq = $#string +1;
print "<br>\nThe number of sequences in the family is <b>$number_of_seq</b><br>\n";

### Create the de Bruijn graph 

my %edge = ();
my %edge_original = (); 
my %nodes_array =();
my %smaller_nodes_array = ();
my %node_cluster=();
my %num_node_cluster =();
my %cluster = ();


for ($i=0;$i<$number_of_seq;$i++){
    $start = 0;
    my $pre_node = "";
    while($start < $len[$i]- $node_len ){

	my $node1 = substr($string[$i],$start,$node_len);
	my $node2 = substr($string[$i],$start+$offset,$node_len);

	if ((!(exists($edge_rec{$node1}{$node2}{$i})) || ($edge_rec{$node1}{$node2}{$i} <= $nsites)) && ($pre_node ne $node2 || $pre_node ne $node1) ){

	    $edge{$node1}{$node2}++;
	    $edge_original{$node1}{$node2}++;
	    $edge_rec{$node1}{$node2}{$i}++;

	}
	 
	$start++;  
        $pre_node = $node1; 	
    }

}

print "Done creating graph\n";

### To make it fuzzy

foreach my $node1 (keys %edge){
    my @chars = split("", $node1);
    my $sum =0;
    my %partial_string =();
    for( my $i=0;$i<=$#chars;$i++){
	$sum += $matrix{$chars[$i]}{$chars[$i]};
    }
#    print "$node1, $sum\n";
    for( my $i=0;$i<=$#chars;$i++){
	my %partial_string2 =();
	foreach my $key (keys %{$sim{$chars[$i]}}){
	    if($i ==0){
		$partial_string{$key}=$matrix{$chars[$i]}{$key};
		
	    }
	    else{
		foreach my $dev_string (keys %partial_string){
		    my $temp_string = $dev_string.$key;
		    $partial_string2{$temp_string} = $partial_string{$dev_string} + $matrix{$chars[$i]}{$key};
		}		
	    }
	}
	if($i !=0){
	    %partial_string = %partial_string2;
	}
    }

    

    foreach my $node2 (keys %{$edge{$node1}}) {

	$edge{$node1}{$node2} = $edge_original{$node1}{$node2};		
	### Find if any other node is very similar to node 2
	my $char_1 = substr($node2,-1*$offset);  	
	
	foreach my $sim_node1 (keys %partial_string){
#	    print "$sim_node1\n";
	    foreach my $inner_key (keys %{$sim{$char_1}}){
		my $sim_node2=substr($sim_node1,$offset).$inner_key;	 
		if(exists($edge_original{$sim_node1}{$sim_node2}) && ($inner_key ne $char_1 || $sim_node1 ne $node1)){  
#		    if (($matrix{$char_1}{$char_2} > 0) && ($key2 ne $node2)){
		    $edge{$node1}{$node2} += $sim_scaling_factor * $edge_original{$sim_node1}{$sim_node2}*($partial_string{$sim_node1}+$matrix{$inner_key}{$char_1})/($sum+$matrix{$char_1}{$char_1});
#		    print "$node1, $node2, $edge_original{$node1}{$node2}, $sim_node1, $sim_node2, $edge{$node1}{$node2}\n"; 
		}	
	    }
	}
    }
}


print "Done making it fuzzy\n";

### Find the edge with maximum weight

$max=0;
foreach $key1 (keys %edge){
    foreach $key2 (keys %{$edge{$key1}}){
	if ($edge{$key1}{$key2}>$max) {$max=$edge{$key1}{$key2}; $maxkey1=$key1;$maxkey2=$key2;}
    }
}


$max_main = $max;

$max_threshold = $max_main *  $scale_max_threshold / 100;
#$threshold = $max * $scale_threshold / 100;
$threshold = $max_main *  $scale_max_threshold * $scale_threshold / 10000;

my @char_maxkey1 = split("", $maxkey1);

my $max_temp=0;

foreach $key (keys %{$edge{$maxkey1}}){
    if ($edge{$maxkey1}{$key}>=$threshold) {
	$reg_ex_thres .= substr($key,$node_len-$offset);
    }
    if($edge{$maxkey1}{$key}>$max_temp){$max_temp = $edge{$maxkey1}{$key};$maxkey2=$key;}
    $reg_ex_all .= substr($key,$node_len-$offset);
}  


if($reg_ex_thres eq ""){
    $look_ahead_flag  = 1;
}
else{
    $look_ahead_flag = 0;
}

#### Create text files as input to Graphviz

$graph_threshold = $max * $scale_graph_threshold / 100 ;
  
=head
open(FD,">graph_input.txt");
open(FD1, ">weights.txt");

foreach $key1 (keys %edge){
    foreach $key2 (keys %{$edge{$key1}})  {
	if ($edge{$key1}{$key2}>= $graph_threshold){
	    print FD "$key1\n";
	    print FD "$key2\n";
	    print FD1 "$edge{$key1}{$key2}\n";
	}
    }
}


close(FD);
close(FD1);

=cut
    
    open(TE, ">>test_into.txt");
print TE "$protein_file\t$number_of_seq\t$triplet";



open(DI, ">detailed_info.txt");    ## to write all the additional info like forward max, backward max, best greedy motif, etc
open(CON, ">consensus_seqs.fasta");

### Traverse the graph starting from the maximum weight edge till edge weight is above threshold

$key1=$maxkey1;
$key2=$maxkey2;
$motif=$key1;
$path_weight=0;
my $motif_count = 0;

foreach $key1 (keys %edge){ 
    foreach $key2 (keys %{$edge{$key1}}){ 
	$edge_dup{$key1}{$key2} = 0;
    }
}

while ($max >= $max_threshold){ 
#    print  "The forward max is $max at $maxkey1 -> $maxkey2 <br>\n";
    print DI "The forward max is $max at $maxkey1 -> $maxkey2 \n"; 
    $break_flag = 0;
    #$look_ahead_flag = 0;
    $jump_count =0;
    $reg_ex_look_ahead = "";
    $motif_look_ahead = "";
    while (($edge{$key1}{$key2}>= $threshold || $look_ahead_flag ==1) && ($break_flag == 0) ){

	if($look_ahead_flag ==1){

	    if($reg_ex_thres eq ""){
                $reg_ex_look_ahead .= "{$reg_ex_all}";
            }
            else {
                $reg_ex_look_ahead .= "{$reg_ex_thres}";
            }

	    $reg_ex_thres = "";
            $reg_ex_all = "";

	    $motif_look_ahead  = $motif_look_ahead.substr($key2,$node_len-$offset);
	#$path_weight += $edge{$key1}{$key2};
	#$edge_dup{$key1}{$key2}=1;
	    $max_inner=0;
	    $max_inner_key = "";
	    #$reg_ex_look_ahead .= "{";
	    foreach $key (keys %{$edge{$key2}}){
		if ($edge{$key2}{$key}>= $threshold) {
		$reg_ex_thres .= substr($key,$node_len-$offset);
	    }
		$reg_ex_all .= substr($key,$node_len-$offset);
		if (($edge{$key2}{$key}>$max_inner) && (!(($key eq $key2) && ($key2 eq $key1))) && ($edge_dup{$key2}{$key} == 0) ){
		    $max_inner = $edge{$key2}{$key}; 
		    $max_inner_key = $key;
		}  
	    }

	}
	else{
	    if($reg_ex_thres eq ""){
                $reg_ex .= "{$reg_ex_all}";
            }
            else {
                $reg_ex .= "{$reg_ex_thres}";
            }

	    $motif = $motif.substr($key2,$node_len-$offset);
	    $path_weight += $edge{$key1}{$key2};
	    $edge_dup{$key1}{$key2}=1;
	    $max_inner=0;
	    $max_inner_key = "";
	    $reg_ex_thres = "";
	    $reg_ex_all = "";
	    foreach $key (keys %{$edge{$key2}}){
		if ($edge{$key2}{$key}>= $threshold) {
		    $reg_ex_thres .= substr($key,$node_len-$offset);
		}
		$reg_ex_all .= substr($key,$node_len-$offset);
		if (($edge{$key2}{$key}>$max_inner) && (!(($key eq $key2) && ($key2 eq $key1))) && ($edge_dup{$key2}{$key} == 0) ){
		    $max_inner = $edge{$key2}{$key};
		    $max_inner_key = $key;
		}
	    }
  
	}
	
	if ($max_inner_key eq "") {
	    $break_flag = 1;
	}
	else {
	    $key1=$key2;
	    $key2=$max_inner_key;
	    if($edge{$key1}{$key2} < $threshold){
		if($jump_count < $jump){
		    $jump_count++;
		    $look_ahead_flag = 1;
		}
		else{
		    $reg_ex .= "($reg_ex_look_ahead)?";
		    $motif .= $motif_look_ahead;
		    $look_ahead_flag =0;
		    $reg_ex_look_ahead = "";
		    $motif_look_ahead = "";
		}
	    }
	}
	if($edge{$key1}{$key2} >= $threshold && $look_ahead_flag == 1){
	    $reg_ex .= $reg_ex_look_ahead;
	    $motif .= $motif_look_ahead;
	    $look_ahead_flag = 0;
	    $jump_count = 0;
	    $reg_ex_look_ahead = "";
	    $motif_look_ahead = "";
	}
    }

#### Looking backward

    $key1=$maxkey1;
    $max_inner_key="";
    $max2=0;
    $reg_ex_thres = "";
    $reg_ex_all = "";

    foreach $alpha (keys %matrix){
	$key=$alpha.substr($key1,0,$overlap); 
	if (exists($edge{$key}{$key1})){
	    if ($edge{$key}{$key1} >=  $threshold){
		$reg_ex_thres= $alpha.$reg_ex_thres; 
	    }
	    $reg_ex_all = $alpha.$reg_ex_all;
	    if  (($edge{$key}{$key1} > $max2)&& ($edge_dup{$key}{$key1} == 0)){
		$max2=$edge{$key}{$key1};
		$max_inner_key=$key;
	    }
	}
    }

    print DI "The backward max is $max2 at $max_inner_key -> $maxkey1\n";
    $key2= $max_inner_key;  

    $break_flag = 0;
    if( $edge{$key2}{$key1}>= $threshold ){
	$look_ahead_flag = 0;
    }
    else{
	$look_ahead_flag = 1;
    }
    $jump_count =0;
    $reg_ex_look_ahead = "";
    $motif_look_ahead = "";

    while (($edge{$key2}{$key1}>= $threshold || $look_ahead_flag ==1 ) && ($break_flag == 0)){

	if($look_ahead_flag ==1){
	    if($reg_ex_thres eq ""){
                $reg_ex_look_ahead = "{$reg_ex_all}".$reg_ex_look_ahead;
            }
            else{
                $reg_ex_look_ahead = "{$reg_ex_thres}".$reg_ex_look_ahead;
            }

	    $reg_ex_thres = "";
            $reg_ex_all = "";

	    $motif_look_ahead = substr($key2,0,$offset).$motif_look_ahead;
	    #$path_weight +=$edge{$key2}{$key1};
	    #$edge_dup{$key2}{$key1}=1;
	    $max_inner=0;
	    $max_inner_key = "";
	    foreach $alpha (keys %matrix){
		$key=$alpha.substr($key2,0,$overlap);
		if (exists($edge{$key}{$key2})){
		    if ($edge{$key}{$key2} >= $threshold){
		    $reg_ex_thres= $alpha.$reg_ex_thres;
		    }
		    $reg_ex_all = $alpha.$reg_ex_all;
		    if (($edge{$key}{$key2} > $max_inner) && ($edge_dup{$key}{$key2} == 0) && (!(($key eq $key2) && ($key2 eq $key1)))){ 
		        $max_inner=$edge{$key}{$key2};
			$max_inner_key=$key; 
		    }
		}
	    }
	    #$reg_ex_look_ahead = "{".$reg_ex_look_ahead;
	}
	else{
	    if($reg_ex_thres eq ""){
                $reg_ex = "{$reg_ex_all}".$reg_ex;
            }
            else{
                $reg_ex = "{$reg_ex_thres}".$reg_ex;
            }



	    $motif = substr($key2,0,$offset).$motif;
	    $path_weight +=$edge{$key2}{$key1};
	    $edge_dup{$key2}{$key1}=1;
	    $max_inner=0;
	    $max_inner_key = "";
	    $reg_ex_thres = "";
            $reg_ex_all = "";

	    foreach $alpha (keys %matrix){
		$key=$alpha.substr($key2,0,$overlap);
		if (exists($edge{$key}{$key2})){
		    if ($edge{$key}{$key2} >= $threshold){
			$reg_ex_thres= $alpha.$reg_ex_thres;
		    }
		   
		    $reg_ex_all= $alpha.$reg_ex_all;
		   
		    if (($edge{$key}{$key2} > $max_inner) && ($edge_dup{$key}{$key2} == 0) && (!(($key eq $key2) && ($key2 eq $key1)))){
			$max_inner=$edge{$key}{$key2};
			$max_inner_key=$key;
		    }
		}
	    }

	}
	if ($max_inner_key eq "") {$break_flag = 1;}
	else {
	    $key1=$key2;
	    $key2=$max_inner_key;
	    if($edge{$key2}{$key1} < $threshold){
                if($jump_count < $jump){
                    $jump_count++;
		    $look_ahead_flag = 1;
                }
		else{
		    $reg_ex = "($reg_ex_look_ahead)?".$reg_ex;
		    $motif = $motif_look_ahead.$motif;
                    $look_ahead_flag =0;
                    $reg_ex_look_ahead = "";
                    $motif_look_ahead = "";
                }
            }

	}
	if($edge{$key2}{$key1} >= $threshold && $look_ahead_flag == 1){
            $reg_ex = $reg_ex_look_ahead.$reg_ex;
            $motif = $motif_look_ahead.$motif;
            $look_ahead_flag = 0;
            $jump_count = 0;
            $reg_ex_look_ahead = "";
            $motif_look_ahead = "";
        }

	
    }
    
    $reg_ex =~ s/\{\}//g;

    if ($motif =~ /$triplet/){
        print TE "\t*$motif";
    }
    else{
        print TE "\t$motif";
    }


    print DI "The motif is: $motif , Path weight: $path_weight\n";
    print "The motif is: $motif , Path weight: $path_weight<br>";
    print CON ">Consensus_$motif_count \n$motif\n";


## The next iteration

    $max=0;
    foreach $key1 (keys %edge){
	foreach $key2 (keys %{$edge{$key1}}){
	    if (($edge{$key1}{$key2}>$max) && ($edge_dup{$key1}{$key2} == 0)){ 
		$max=$edge{$key1}{$key2}; $maxkey1=$key1;$maxkey2=$key2;}
	}
    }
  

    if ($max >= $max_threshold)  { 
		
	@char_maxkey1 = split("", $maxkey1);

	$reg_ex_all="";
	$reg_ex_thres="";

	foreach $key (keys %{$edge{$maxkey1}}){
	    if ($edge{$maxkey1}{$key}>= $threshold) {
		$reg_ex_thres .= substr($key,$node_len-$offset);
	    }
	    $reg_ex_all .= substr($key,$node_len-$offset);
	}
		
	if($reg_ex_thres eq ""){
	    $look_ahead_flag  = 1;
	}
	else{
	    $look_ahead_flag = 0;
	}

	$key1=$maxkey1;
	$key2=$maxkey2;
	$motif=$key1;
	$path_weight=0;
    }

    $motif_count++;

}


print TE "\n";
close TE;

my $script_end_time = new Benchmark;

my $diff = timediff($script_end_time, $start_time);

my $running_time = timestr($diff);

open (TM, ">>time.txt");

print TM "$protein_file\t$running_time\n"; 

close(TM);


close(DI);


exit;



########################

sub open_file

{ my $file_name = $_[0];  

  if ($file_name)
  {if(!open(FD,"$file_name")) {print "\nCouldn't open $file_name \n"; exit;}
  }
  else
  {print "Sorry...Please specify the file name \n";exit;}
    
   my @file_data=<FD>;
   close(FD);   

   return(@file_data);

}

######################

sub process_file_data
  
{ my ($i,$j,$temp);
  my @strings=();
  my @file_data = @_;
 
  $j=-1;
  for($i=0;$i<=$#file_data;$i++)
  {
   $temp=$file_data[$i];
   #Check for FASTA format
   if (substr($temp,0,1) ne ">")
   {
    # Remove any newline character
    chomp $temp;
    $strings[$j]=$strings[$j].$temp;
   }
   else {$j++; $strings[$j]="";}
  }  

  for($i=0;$i<=$j;$i++)
{

  # Remove any spaces
  $strings[$i] =~ s/\s+//g;

  # convert the strings to upper case
  $strings[$i] = uc $strings[$i];
  
  # Remove long runs of XX's
  $strings[$i] =~ s/XXXX//g;

  


  # Error check  
  if($strings[$i] =~ /[^ARNDCQEGHILKMFPSTWVYXZB]/)
  { print "<font face=\"verdana\" size=3 color= #FF2266 >";
    print "Sorry ! Your input file $protein_file ($i) does not contain a valid protein sequence!!!\n";
    print "</font>"; exit 2;
  }
   
}
return @strings;
}

##################


sub read_matrix_file{

# read the substution matrix into a hash 
 
    my @matrix_data = @_;
    my %matrix = ();

    my $m=-1;
    do # Skip the comment and blank lines
    {
	$m++;
	chomp $matrix_data[$m];
      
    } until  ((substr($matrix_data[$m],0,1) ne "#") && ($matrix_data[$m] ne ""));
    
    # Read the actual matrix into array
    
    my $k = $m+1;
    
    my @matrix_letters = split(' ',$matrix_data[$m]);
    
    for($i=0;$i<=23;$i++)
    {  
	chomp $matrix_data[$k];
	my @row = split(/\s+/,$matrix_data[$k]);
	$alpha = shift @row; 
	for($j=0;$j<=23;$j++){
	    $matrix{$alpha}{$matrix_letters[$j]} = $row[$j];
	}
	$k++;
    }
    return(%matrix);

}


##########


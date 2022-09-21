#!/usr/bin/perl -w 

# copied this from cgi-bin/deBruijn on 04/17/04

use Getopt::Long;
#use strict;

my ($i,$j,$m,$k,$alpha,$temp,$max,$max2,$max_main,$start,$key,$key1,$key2,$seq_len,$sequence,$maxkey1,$maxkey2,$motif,$path_weight, $reg_ex, $add, $break_flag);
my ($max_inner,$max_inner_key);
my ($subs_matrix, $protein_file, $seq_file);
my (@config_data,@str1,@str2,@seq,@len,%edge_dup,%edge_rec,%edge_back);
my ($node_len,$overlap,$offset,$threshold,$max_threshold, $graph_threshold);
my ($scale_threshold,$scale_max_threshold, $scale_graph_threshold,$nsites);
my ($score_messy, $score_per, $score_neg, $messy_match, $score);
my ($look_ahead_flag, $jump_count, $reg_ex_all, $reg_ex_thres, $reg_ex_look_ahead, $motif_look_ahead,$jump);
my ($reg_ex_temp, $reg_ex_middle, $option_count,@temp_maxkey1,$exact,$strict);

GetOptions ("sub=s" => \$subs_matrix , "fam=s" => \$protein_file,  "seq=s" => \$seq_file,  "nlen=s" => \$node_len, "ov=s" =>\$overlap, "th=s" =>\$scale_threshold,  "mth=s" =>\$scale_max_threshold, "gth=s" =>\$scale_graph_threshold, "n=s" =>\$nsites, "jump=s" => \$jump);
   
# "sub" needs the name of substitution matrix file as the argument
# "fam" needs the name of the input file containing the protein family sequences in FASTA format as the argument
# "seq" needs the name of the file containing the protein sequence to be tested

$|=1;

#### Set the default values

if (!(defined($node_len))) {$node_len = 4;}
if (!(defined($overlap))) {$overlap = 3;}
if (!(defined($scale_threshold))) {$scale_threshold = 90;}
if (!(defined($scale_max_threshold))) {$scale_max_threshold = 99;}
if (!(defined($scale_graph_threshold))) {$scale_graph_threshold = 60;}
if (!(defined($nsites))) {$nsites = 1;}
if (!(defined($jump))){$jump = 0;}


$overlap = $node_len - 1;
$offset = $node_len - $overlap;
$padding = 4;

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


my $number_of_seq = $#string +1;
#print "<br>\nThe number of sequences in the family is $number_of_seq<br><br>\n";

### Create the de Bruijn graph 

my %edge = ();
my %edge_original = (); 
my %nodes_array =();
my %smaller_nodes_array = ();
my %node_cluster=();
my %num_node_cluster =();
my %cluster = ();


#print "Constructing the graph ...<br>";
#print "Processing Sequence ";
for ($i=0;$i<$number_of_seq;$i++){
    $start = 0;
    my $seq_no = $i +1;
   # print ".. ".$seq_no;
    my $pre_node = "";

    while($start < $len[$i]- $node_len ){

	my $node1 = substr($string[$i],$start,$node_len);
	my $node2 = substr($string[$i],$start+$offset,$node_len);

        # Reduce size of alphabet using Sun's classes
#	my $con_node1 = &collapse_node($node1);
#	my $con_node2 = &collapse_node($node2);
	
	if ((!(exists($edge_rec{$node1}{$node2}{$i})) || ($edge_rec{$node1}{$node2}{$i} <= $nsites)) && ($pre_node ne $node2 || $pre_node ne $node1) ){

	    &create_masked_graph($node1,$node2,$i);

#	    $edge{$node1}{$node2}++;
	    $edge_original{$node1}{$node2}++;
#	    $edge_rec{$node1}{$node2}{$i}++;

	}
	 
	$start++;  
        $pre_node = $node1; 	
    }

}

print "<br>Done creating graph<br>";


=head
### To make it fuzzy

foreach $node1 (keys %edge){
        
    foreach $node2 (keys %{$edge{$node1}}) {
		
	### Find if any other node accessible from node 1 is very similar to node 2
	foreach $key2 (keys %{$edge{$node1}}){
	    my $char_1 = substr($node2,-1*$offset);
	    my $char_2 = substr($key2,-1*$offset);
	    
	    if (($matrix{$char_1}{$char_2} > 0) && ($key2 ne $node2)){
		$edge{$node1}{$node2} += $edge_original{$node1}{$key2}*(0.5+($matrix{$char_1}{$char_2})/(2*($matrix{$char_1}{$char_1})));
	    }
	
	}
    }
}

=cut


### Find the edge with maximum weight

my @rec_array=&first_find_max();

$maxkey1 = $rec_array[0];
$maxkey2 = $rec_array[1];
$max = $rec_array[2];


$max_main = $max;

$max_threshold = $max_main *  $scale_max_threshold / 100;
#$threshold = $max * $scale_threshold / 100;
$threshold = $max_main *  $scale_max_threshold * $scale_threshold / 10000;

$max_temp=0;
$look_ahead_flag = 0;



=head
#### Create text files as input to Graphviz

$graph_threshold = $max * $scale_graph_threshold / 100 ;
  
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

open(OP,">motif.txt");
open(DE, ">test_rel.txt");

open(DI, ">detailed_info.txt");    ## to write all the additional info like forward max, backward max, best greedy motif, etc

print OP "$protein_file\n";

### Traverse the graph starting from the maximum weight edge till edge weight is above threshold

$key1=$maxkey1;
$key2=$maxkey2;
#$motif=&best_of_two_nodes_forward($key1,$key2);;
$motif=substr($key1,0,$offset);
@motif_array = ();

$path_weight=0;

foreach my $local_key1 (keys %edge){ 
    foreach my $local_key2 (keys %{$edge{$local_key1}}){ 
	$edge_dup{$local_key1}{$local_key2} = 0;
    }
}

while ($max >= $max_threshold){ 
#   print  "The forward max is $max at $maxkey1 -> $maxkey2 <br>\n";
    print DI "The forward max is $max at $maxkey1 -> $maxkey2 \n"; 
    $break_flag = 0;
    $jump_count =0;
    $motif_look_ahead = "";

    while (($edge{$key1}{$key2}>= $threshold || $look_ahead_flag ==1) && ($break_flag == 0) ){
	if($look_ahead_flag ==1){
#	    $motif_look_ahead = $motif_look_ahead.substr($key2,$node_len-$offset);
	    $motif_look_ahead = $motif_look_ahead.&best_of_two_nodes_forward($key1,$key2);	
	    $max_inner=0;
	    $max_inner_key = "";
	}
	else{
#	    $motif = $motif.substr($key2,$node_len-$offset);
	    $motif = $motif.&best_of_two_nodes_forward($key1,$key2);
	    $path_weight += $edge{$key1}{$key2};
	    foreach $key (keys %{$edge{$key1}}){
		$edge_dup{$key1}{$key}=1;
		foreach $key_back (keys %{$edge_back{$key}}){
		    $edge_dup{$key_back}{$key}=1;
		}

	    }
	    foreach $key (keys %{$edge_back{$key2}}){
                $edge_dup{$key}{$key2}=1;
		foreach $key_back (keys %{$edge{$key}}){
                    $edge_dup{$key}{$key_back}=1;
                }

            }

	}
	
	$max_inner_key = &find_max($key1,$key2,1,\%edge);

	if ($max_inner_key eq "") {
	    $break_flag = 1;
	}
	else {
	    push(@motif_array, $key2);
	    $key1=$key2;
	    $key2=$max_inner_key;
	    if($edge{$key1}{$key2} < $threshold){
		if($jump_count < $jump){
		    $jump_count++;
		    $look_ahead_flag = 1;
		}
		else{
		    $motif .= $motif_look_ahead;
		    $look_ahead_flag =0;
		    $motif_look_ahead = "";
		}
	    }
	
	    if($edge{$key1}{$key2} >= $threshold && $look_ahead_flag == 1){
		$motif .= $motif_look_ahead;
		$look_ahead_flag = 0;
		$jump_count = 0;
		$motif_look_ahead = "";
	    }
	}
    }

    $motif .= substr($key1,$offset);
    
#### Looking backward

    unshift(@motif_array, $maxkey1);

    $key1 = $maxkey1;
    $max_inner_key="";
    $max2=0;

    foreach $key (keys %{$edge_back{$key1}}){
	if ($edge_back{$key1}{$key}>=$max2) {
	    if ($edge_back{$key1}{$key}>$max2) {
		$max2=$edge_back{$key1}{$key};
		@max_array = ();
		push(@max_array, $key);
	    }
	    else{
		push(@max_array, $key);
	    }
	}

    }


### If multiple edges have maximum weight, we want to pick the more specific one.


    $min_x =5;
    $key2= "";

#    print "the contents of max_array :\n";

    for(my $i = 0; $i<= $#max_array; $i++){
	if ($max_array[$i] =~ /x..x|x.x.|.x.x|.xx./){
	    $x=2;
	}
	else{
	    $x=1;
	}
	
	if(($x == $min_x)&& ($max_array[$i]=~ /...x/)){
            $key2 = $max_array[$i];
            $min_x = $x;

        }


	if($x < $min_x){
	    $key2 = $max_array[$i];
	    $min_x = $x;

	}
	
	#print "$max_array[$i] <- $key1\n";
    }

    if($key2 ne ""){

	print DI "The backward max is $max2 at $key2 <- $maxkey1\n";
	
	$break_flag = 0;
	
 	if( $edge_back{$key1}{$key2}>= $threshold ){
	    $look_ahead_flag = 0;
	}
	else{
	    $look_ahead_flag = 1;
	}
	
	$jump_count = 0;
	$motif_look_ahead = "";
	
	while (($edge_back{$key1}{$key2}>= $threshold || $look_ahead_flag ==1 ) && ($break_flag == 0)){

	    if($look_ahead_flag ==1){
		
		$motif_look_ahead = substr($key2,0,$offset).$motif_look_ahead;

#		$motif_look_ahead = &best_of_two_nodes_back($key1,$key2).$motif_look_ahead; 
	    }
	    else{   
		$motif = substr($key2,0,$offset).$motif;
#		$motif = &best_of_two_nodes_back($key1,$key2).$motif;
		$path_weight +=$edge_back{$key1}{$key2};

		foreach $key (keys %{$edge_back{$key1}}){
		    $edge_dup{$key}{$key1}=1;
		    foreach $key_back (keys %{$edge{$key}}){
			$edge_dup{$key}{$key_back}=1;
		    }


		}
		foreach $key (keys %{$edge{$key2}}){
                    $edge_dup{$key2}{$key}=1;
		    foreach $key_back (keys %{$edge_back{$key}}){
                        $edge_dup{$key_back}{$key}=1;
                    }

                }


	    }
	    
	    $max_inner_key = &find_max($key1,$key2,2,\%edge_back);

	    if ($max_inner_key eq "") {$break_flag = 1;}
	    else {
		unshift(@motif_array, $key2);
		$key1=$key2;
		$key2=$max_inner_key;
		if($edge_back{$key1}{$key2} < $threshold){
		    if($jump_count < $jump){
			$jump_count++;
			$look_ahead_flag = 1;
		    }
		    else{
			$motif = $motif_look_ahead.$motif;
			$look_ahead_flag =0;
			$reg_ex_look_ahead = "";
			$motif_look_ahead = "";
		    }
		}		
		if($edge_back{$key1}{$key2} >= $threshold && $look_ahead_flag == 1){
		
		    $motif = $motif_look_ahead.$motif;
		    $look_ahead_flag = 0;
		    $jump_count = 0;
		    $motif_look_ahead = "";
		}
	    }
		
	}
    }
    
#    $motif= substr($key1,0,$overlap).$motif;

    $newmotif = substr($motif_array[0],0,$offset);

    for($i=0;$i<$#motif_array;$i++){
	$newmotif .= &best_of_two_nodes_forward($motif_array[$i],$motif_array[$i+1]);
    }

    $newmotif .= substr($motif_array[$#motif_array],$offset);

    $newmotif =~ s/^x//g;
    $newmotif =~ s/x$//g;

#    print "<p>The motif is $newmotif <br>";

    my $find_motif =$newmotif;
    $find_motif =~ s/x/\./g;
    $find_motif .= ".{$padding}";
    $find_motif = ".{$padding}".$find_motif;

    print DI "The motifs in sequences are\n";

    foreach my $seq (@string){
	my $temp = $seq;
	while($temp =~ s/(...$find_motif...)//){
	    print DE "$1 1\n";
	}
    }
    
    print DE "\n";
    
## The next iteration

    $max=0;
    foreach $key1 (keys %edge){
	foreach $key2 (keys %{$edge{$key1}}){
	    if (($edge{$key1}{$key2}>$max) && ($edge_dup{$key1}{$key2} == 0)){ 
		$max=$edge{$key1}{$key2}; $maxkey1=$key1;$maxkey2=$key2;}
	}
    }
  
    if ($max >= $max_threshold)  { 
		
	if($max < $threshold ){
	    $look_ahead_flag  = 1;
	}
	else{
	    $look_ahead_flag = 0;
	}

	$key1=$maxkey1;
	$key2=$maxkey2;
#	$motif=&best_of_two_nodes_forward($key1,$key2);
	$motif=substr($key1,0,$offset);
	@motif_array = ();
	$path_weight=0;
    }
}

system("perl rank_motifs.pl $protein_file $subs_matrix");

close(OP);
close(DI);
close(DE);

##### If seq specified, check if it belongs to the family

if ($seq_file){

    @str2 = &open_file($seq_file); 
    @seq = &process_file_data(@str2); 
    $sequence = $seq[0]; 
    $seq_len = length $sequence; 
    $score = 0;
    $start = 0;
    $score_messy = 0;
    $score_neg=0;
    $score_per = 0;
   
    while($start < $seq_len - $node_len ){
	my $node1 = substr($sequence,$start,$node_len);
	my $node2 = substr($sequence,$start+$offset,$node_len);
	if (exists($edge{$node1}{$node2})) {$score_per++; $score += $edge{$node1}{$node2};}
#if (!(exists $edge{$node1}{$node2})) {$score_neg++;}
  
	### To make it messy
  
	$messy_match = 0; 
	foreach $key2 (keys %{$edge{$node1}}){  
	    
	    ### Find if any other node accessible from node 1 is very similar to node 2

	    my $char_1 = substr($node2,-1*$offset);   
	    my $char_2 = substr($key2,-1*$offset);
 
	    if (($matrix{$char_1}{$char_2} > 0) && ($key2 ne $node2)){ 
		$messy_match=1;
		$add = $edge{$node1}{$key2}*(1.5+($matrix{$char_1}{$char_2})/(2*($matrix{$char_1}{$char_1})));
		$score_messy += $add;   
	    }
	}
	if (($messy_match == 0) &&(!(exists $edge{$node1}{$node2}))){$score_neg++;}
	$start++;
    }  

    print "The score is $score \n The messy score is $score_messy \n Perfect matches $score_per \n No matches $score_neg \n Length $seq_len";

    #if (($score_per/($seq_len/$offset) > 0.2) && ($score+$score_messy > $score_per * 0.1 * ($number_of_seq)))
#if ((($seq_len-$score_neg)/($seq_len/$offset) > 0.15) && (($score)/($score_per) >  0.1 * ($number_of_seq)))    
if ((($seq_len-$score_neg)/($seq_len/$offset) > 0.15) || (($score)/($score_per) > 3))

{print "\nThe protein seems to belong to the family \n";}
    else {print "\nThe protein does not seem to belong to the family \n";}
	
}

#print "</font>";
exit;


=head


print "<p>Here is a link to the output file ";
print "<a href='http://biokdd.informatics.indiana.edu/rpatward/L519/writable/output.txt'> Output </a><br>";

=cut

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

sub collapse_node{

# Reduce size of alphabet using Sun's classes  
 
  #  W F Y 
  #  E K Q R
  #  L I V M  
  #  A S T                                    
 
    my  $coll_node = $_[0];

    $coll_node =~ s/[WFY]/F/g;
    $coll_node =~ s/[EKQR]/E/g;
    $coll_node =~ s/[LIVM]/I/g;
    $coll_node =~ s/[AST]/A/g;
    
    return($coll_node);

}

###############


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

sub cluster_node{

 # Node clustering part
    my %hash = %{$_[0]};
    my $node = $_[1];
    
    my @split_node = split("", $node);
    my $cluster_flag = 0;
    my $simi_max = $clus_threshold;
    my $simi_max_node = "";
    
    foreach my $key (keys %hash){
	my @temp_split = split("",$key);
	my $simi_score =0;
	for($letter = 0;$letter <= $#split_node;$letter++){
	    $simi_score += $matrix{$split_node[$letter]}{$temp_split[$letter]} ;
	}
	if($simi_score >= $simi_max){
	    $simi_max = $simi_score;
	    $simi_max_node = $key;
	    $cluster_flag = 1;
	}
    }
    
    if($cluster_flag ==0){  ## Create new cluster
	$node_cluster{$node}{$node}++;
	$num_node_cluster{$node}++;
	$cluster{$node} = $node;
    }
    else{
	$node_cluster{$cluster{$simi_max_node}}{$node}++;
	$num_node_cluster{$cluster{$simi_max_node}}++;
	$cluster{$node} = $cluster{$simi_max_node};
    }
    
    $nodes_array{$node}++;
    
return; 
 
}

################

sub switch_to_smaller_array{

# Copy the keys of clusters with more than 2 members 

    foreach my $key (keys %num_node_cluster){
	if($num_node_cluster{$key} > 2){
	    $smaller_nodes_array{$key}++;
	}
    }

    return;
}


#################

sub sun_cluster{

# Find the highest weight edge between sun's clusters

    $max_edge_sun_clus = 0;

    foreach my $key1 (keys %edge_sun_clus){
	foreach my $key2 (keys %{$edge_sun_clus{$key1}}){
	    if($edge_sun_clus{$key1}{$key2} > $max_edge_sun_clus){
		$max_edge_sun_clus = $edge_sun_clus{$key1}{$key2}; 
		$max_key1_sun_clus = $key1;
		$max_key2_sun_clus = $key2;
	    }
	}
    }
    
    
    print " max_edge_sun_clus: $max_edge_sun_clus\n max_key1_sun_clus $max_key1_sun_clus \n  max_key2_sun_clus $max_key2_sun_clus \n ";
    
    $max_sun_clus =0;

    foreach $key1 (keys %{$node_sun_clus{$max_key1_sun_clus}}){
	if($node_sun_clus{$max_key1_sun_clus}{$key1} > $max_sun_clus){
	    $max_sun_clus = $node_sun_clus{$max_key1_sun_clus}{$key1};
	    $max_node1_sun_clus = $key1;
	}
    }

    $max=0;

    $maxkey1=$max_node1_sun_clus;

    foreach $key1 (keys %{$edge{$max_node1_sun_clus}}){
        if ($edge{$max_node1_sun_clus}{$key1}>$max) {
	    $max=$edge{$max_node1_sun_clus}{$key1};
	    
	    $maxkey2=$key1;
	}
    }
    
    return;    
    
}

#######

sub cluster{
    $clus_max=0;
    foreach $key1 (keys %edge_clus){
	foreach $key2 (keys %{$edge_clus{$key1}}){
	    if ($edge_clus{$key1}{$key2}>$clus_max) {
		$clus_max=$edge_clus{$key1}{$key2}; 
		$clus_maxkey1=$key1;
		$clus_maxkey2=$key2;
	    }	   
	    if ($edge_clus{$key1}{$key2}>45) {
		print "clus_max is $edge_clus{$key1}{$key2} , $key1 , $key2 \n";
	    }
	}
    }
    
    print "clus_max is $clus_max, $clus_maxkey1 , $clus_maxkey2 \n";
 
    foreach $key (keys %{$node_cluster{$clus_maxkey1}}){
	print "$key $node_cluster{$clus_maxkey1}{$key}\n";
    }

    return;   
}

######

sub create_masked_graph{

  my  $n1 = $_[0];
  my  $n2= $_[1];

  my $seq_number =$_[2];

#  print "ni and n2 $n1,$n2 \n";

  my @n1_array = split("",$n1);
  my @n2_array = split("",$n2);

  my %node1_masks = ();
  my %node2_masks = ();

### Applying the masks

# The xAAx mask
  my $n1_xAAx = 'x'.$n1_array[1].$n1_array[2].'x';
  my $n2_xAAx = 'x'.$n2_array[1].$n2_array[2].'x';

  $node1_masks{$n1_xAAx}++;
  $node2_masks{$n2_xAAx}++;
 
#The AxAA mask
  my $n1_AxAA = $n1_array[0].'x'.$n1_array[2].$n1_array[3];
  my $n2_AxAA = $n2_array[0].'x'.$n2_array[2].$n2_array[3];

  $node1_masks{$n1_AxAA}++;
  $node2_masks{$n2_AxAA}++;

#The AAxA mask
  my $n1_AAxA = $n1_array[0].$n1_array[1].'x'.$n1_array[3];
  my $n2_AAxA = $n2_array[0].$n2_array[1].'x'.$n2_array[3];

  $node1_masks{$n1_AAxA}++;
  $node2_masks{$n2_AAxA}++;


#The AxAx mask
  my $n1_AxAx = $n1_array[0].'x'.$n1_array[2].'x';
  my $n2_AxAx = $n2_array[0].'x'.$n2_array[2].'x';

  $node1_masks{$n1_AxAx}++;
  $node2_masks{$n2_AxAx}++;

#The xAxA mask
  my $n1_xAxA = 'x'.$n1_array[1].'x'.$n1_array[3];
  my $n2_xAxA = 'x'.$n2_array[1].'x'.$n2_array[3];

  $node1_masks{$n1_xAxA}++;
  $node2_masks{$n2_xAxA}++;


#The AxxA mask

  my $n1_AxxA = $n1_array[0].'xx'.$n1_array[3];
  my $n2_AxxA = $n2_array[0].'xx'.$n2_array[3];

  $node1_masks{$n1_AxxA}++;
  $node2_masks{$n2_AxxA}++;


#The AAAx mask

  my $n1_AAAx = $n1_array[0].$n1_array[1].$n1_array[2].'x';
  my $n2_AAAx = $n2_array[0].$n2_array[1].$n2_array[2].'x';

  $node1_masks{$n1_AAAx}++;
  $node2_masks{$n2_AAAx}++;

#The xAAA mask

  my $n1_xAAA = 'x'.$n1_array[1].$n1_array[2].$n1_array[3];
  my $n2_xAAA = 'x'.$n2_array[1].$n2_array[2].$n2_array[3];

  $node1_masks{$n1_xAAA}++;
  $node2_masks{$n2_xAAA}++;


## Creating the graph using masked nodes
  
  foreach my $key1 (keys %node1_masks){
      foreach my $key2 (keys %node2_masks){
	  if ((!(exists($edge_rec{$key1}{$key2}{$seq_number})) || ($edge_rec{$key1}{$key2}{$seq_number} < $nsites)) ) {
	      $edge{$key1}{$key2}++;
	      $edge_rec{$key1}{$key2}{$seq_number}++;
	      $edge_back{$key2}{$key1}++;

	  }
      }
  }
  
  return;
  
}


#####

sub find_max{

    my $local_key1 = $_[0];
    my $local_key2 = $_[1];
    my %array= %{$_[3]};

    my $local_max=0;
    my @max_array = ();

    my $direction = $_[2];

    if ($direction == 1){

	foreach my $key (keys %{$array{$local_key2}}){
	    if (($array{$local_key2}{$key} >= $local_max) && ($edge_dup{$local_key2}{$key} == 0) && (!(($key eq $local_key2) && ($local_key2 eq $local_key1)))){
		if ($array{$local_key2}{$key}>$local_max) {
		    $local_max= $array{$local_key2}{$key};
		    @max_array = ();
		    push(@max_array, $key);
		}
		else{
		    push(@max_array, $key);
		}
	    }
	    
	}
    }
    else{
	foreach my $key (keys %{$array{$local_key2}}){
            if (($array{$local_key2}{$key} >= $local_max) && ($edge_dup{$key}{$local_key2} == 0) && (!(($key eq $local_key2) && ($local_key2 eq $local_key1)))){
                if ($array{$local_key2}{$key}>$local_max) {
                    $local_max= $array{$local_key2}{$key};
                    @max_array = ();
                    push(@max_array, $key);
                }
                else{
                    push(@max_array, $key);
                }
            }

	}
    }
    
	

### If multiple edges have maximum weight, we want to pick the more specific one.


    my $min_x =5;
    $key= "";
    my $x;

    for(my $i = 0; $i<= $#max_array; $i++){
        if ($max_array[$i] =~ /x..x|x.x.|.x.x|.xx./){

            $x=2;
        }
        else{
            $x=1;
        }

	if($x == $min_x && !($max_array[$i] =~ /...x/)){
	    $key = $max_array[$i];
	    $min_x = $x;
	}

        if($x < $min_x){
            $key = $max_array[$i];
	    $min_x = $x;
        }
        
    }


    return($key);
}

sub best_of_two_nodes_forward{

    my $one = $_[0];
    my $two = $_[1];

    my @one_array = split("",$one);
    my @two_array = split("",$two);

#    my $string = substr($one,0,$offset);

    my $string = "";

#    for(my $i=$offset;$i<=$#one_array;$i++){

    my $i = $offset;

	if (($one_array[$i] eq "x" && $two_array[$i-1] eq "x")||($one_array[$i] ne "x" && $two_array[$i-1] ne "x")){
	    $string .= $one_array[$i];	    
	}
	elsif($one_array[$i] ne "x" && $two_array[$i-1] eq "x"){
	    
	    $string .= $one_array[$i];

	}
	elsif($one_array[$i] eq "x" && $two_array[$i-1] ne "x") {
	    $string .= $two_array[$i-1];
	}
#    }


    return($string);    
}

sub best_of_two_nodes_back{

    my $one = $_[0];
    my $two = $_[1];

    my @one_array = split("",$one);
    my @two_array = split("",$two);

#    my $string = substr($one,0,$offset);                                                                                                     

    my $string = "";

    for(my $i=$offset;$i<=$#one_array;$i++){                                                                                                 

#    my $i = $offset;

	if (($one_array[$i] eq "x" && $two_array[$i-1] eq "x")||($one_array[$i] ne "x" && $two_array[$i-1] ne "x")){
	    $string .= $one_array[$i];
	}
	elsif($one_array[$i] ne "x" && $two_array[$i-1] eq "x"){

	}
    }
}



sub first_find_max{
### Find the edge with maximum weight

    my @max_array = ();
    my @ret_array = ();
    my $max=0;
    my $x;

    foreach $key1 (keys %edge){
	foreach $key2 (keys %{$edge{$key1}}){
	    if ($edge{$key1}{$key2}>=$max) {
		if ($edge{$key1}{$key2}>$max) {
		    $max=$edge{$key1}{$key2};
		    @max_array = ();
		    push(@max_array, $key1, $key2);
		}
		else{
		    push(@max_array, $key1, $key2);
		}
	    }

	}
    }

### If multiple edges have maximum weight, we want to pick the more specific one.
    my $min_x =5;

#print "the contents of max_array :\n";

    for(my $i = 0; $i< $#max_array; $i+=2){
	if ($max_array[$i] =~ /x..x|x.x.|.x.x|.xx./){
	    if($max_array[$i+1] =~ /x..x|x.x.|.x.x|.xx./){
		$x=4;
	    }
	    else{
		$x=3;
	    }
	}
	else{
	    if($max_array[$i+1] =~ /x..x|x.x.|.x.x|.xx./){
		$x=3;
	    }
	    else{
		$x=2;
	    }
	}

	if(($x == $min_x)&&($max_array[$i+1])){
	    $maxkey1 = $max_array[$i];
	    $maxkey2 = $max_array[$i+1];
	    $min_x = $x;
	}

	if($x < $min_x){
	    $maxkey1 = $max_array[$i];
	    $maxkey2 = $max_array[$i+1];
	    $min_x = $x;
	}

    }

    $ret_array[0]=$maxkey1;
    $ret_array[1]=$maxkey2;
    $ret_array[2]=$max;

    return(@ret_array);

}

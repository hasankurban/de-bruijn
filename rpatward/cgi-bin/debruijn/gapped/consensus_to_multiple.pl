## Takes in the pair wise alignment of each input sequence against the consensus and converts it to a multiple alignment

use strict;
use warnings;

my $pairwise_alignment_file = $ARGV[0];  ## alignment.out
my $consensus_seq_file = $ARGV[1];  ## consensus_seq.fasta

my $similarity_matrix_file = $ARGV[2]; ## BLOSUM62

#my $writable_dir_path = "/var/www/cgi-bin/rpatward/debruin/approx/writable5";

open(PAL, "$pairwise_alignment_file");
open(CN, "$consensus_seq_file");


my @consensus_seqs = ();

## read in the consensus sequence file

while( my $line = <CN>){
    my $line2 = <CN>;
    chomp $line2;
    push(@consensus_seqs,$line2);
}

close CN;

open(OP, ">final_output.out");

my %aln_seq = ();
my %aln_con = ();

my %con_start = ();

my $consensus;

while (my $line = <PAL>){
    
    chomp $line;

    if($line ne "" && $line ne "//"){

	if ($line =~ /\*Consensus_(\d+)\*/) {  ## get corr consensus seq

	    $consensus = $consensus_seqs[$1];
	}

	else{

	    ### Process concensus line
	    
	    my @temp = split(/\s+/, $line);
	    
	    my $con_part = $temp[1];
	    
	    $con_part =~ s/-//g;
	    
	    my $start = index($consensus, $con_part);
	    
	    ### Process input seq line
	    
	    my $line2 = <PAL>;
	    chomp $line2;
	    
	    my @temp2 = split(/\s+/, $line2);
	    
	    $aln_seq{$temp2[0]} = $temp2[1];
	    
	    $aln_con{$temp2[0]} = $temp[1];
	    
	    $con_start{$temp2[0]} = $start;
	    
	    ## Left pad
	    $aln_seq{$temp2[0]} = "." x ($start).$aln_seq{$temp2[0]};
	    $aln_con{$temp2[0]} = "." x ($start).$aln_con{$temp2[0]};
	}   
    }
    else{

	if($line eq "//"){	

	    my $consensus_new = $consensus;
	    
	    my %gaps = ();
	    
	    my %seq_gap_count =();
	    
	    foreach my $seq_id2 (keys %aln_con){
		my $count = 0;
		my @positions = split("", $aln_con{$seq_id2});
#		for(my $i=$con_start{$seq_id2}; $i<$#positions;$i++){
		for(my $i=0; $i<$#positions;$i++){
		    if($positions[$i] eq "-"){
			$gaps{$count}{$seq_id2}++;
		    }
		    else{
			$count++;
		    }
		}
		$seq_gap_count{$seq_id2} = 0;
	    }
	    
	    my %pos_gap_count = ();
	    
	    foreach my $pos (sort{$a<=>$b} keys %gaps){
		$pos_gap_count{$pos}=0;
		foreach my $seq_id (keys %{$gaps{$pos}}){
		    if($gaps{$pos}{$seq_id} > $pos_gap_count{$pos}){
			$pos_gap_count{$pos}=$gaps{$pos}{$seq_id};
		    }
		}
	    }
	    
	    
	    
	    my $con_gap_count = 0;
	    
	    foreach my $pos (sort{$a<=>$b} keys %gaps){
		
		$consensus_new=insert($consensus_new,$pos+$con_gap_count,$pos_gap_count{$pos});
		
		foreach my $seq_id3 (sort keys %aln_seq){

		    if(length($aln_seq{$seq_id3}) > $pos + $seq_gap_count{$seq_id3}){
			if (!(exists($gaps{$pos}{$seq_id3}))){
			    $aln_seq{$seq_id3} = insert($aln_seq{$seq_id3},$pos + $con_gap_count , $pos_gap_count{$pos});
#			    $aln_seq{$seq_id3} = insert($aln_seq{$seq_id3},$pos + $seq_gap_count{$seq_id3} , $pos_gap_count{$pos});

			    $seq_gap_count{$seq_id3}+= $pos_gap_count{$pos};
			}
			else{
			    if($gaps{$pos}{$seq_id3} < $pos_gap_count{$pos}){
				$aln_seq{$seq_id3} = insert($aln_seq{$seq_id3},$pos + $con_gap_count , $pos_gap_count{$pos}-$gaps{$pos}{$seq_id3});
				$seq_gap_count{$seq_id3}+= $pos_gap_count{$pos}-$gaps{$pos}{$seq_id3};
			    }
			}
		    }
		}
		$con_gap_count +=$pos_gap_count{$pos};
	    }

	    my $length = length($consensus_new);

	    ## Right pad 

	    foreach my $seq_id4 (keys %aln_seq){
		my $length_aln = length $aln_seq{$seq_id4};
		$aln_seq{$seq_id4} = $aln_seq{$seq_id4}. "-" x ($length-$length_aln)  ;

	    }

	    print OP "consensus           \t$consensus_new\n\n";
	    
	    my $i =0;

	    my @seq_array = ();

	    foreach my $seq_id (sort keys %aln_seq){
		$aln_seq{$seq_id} =~ s/\./-/g;
		
		## pad spaces to make all seq_ids same length

		my $length_id = length $seq_id;

		my $id;

		if($length_id > 20){
		    $id = substr($seq_id,0,20);
		}
		else{
		    $id = $seq_id. " " x (20 - $length_id)  ;
		}

		print OP "$id\t$aln_seq{$seq_id}\n";
		my @temp = split("", $aln_seq{$seq_id});
		push(@{$seq_array[$i]}, @temp);
		$i++;
    
	    }

	    ### calculate Sum-of-Pairs score for multiple alignment

	    my $score = sum_of_pairs($similarity_matrix_file, \@seq_array);

	    print OP "\nSOP score is $score\n";
	 
	    %aln_seq = ();
	    %aln_con = ();
	    print OP "\n================\n";
   
	}
    }
}

close PAL;
close OP;

exit;

### subroutines ###

sub insert {

    my ($sequence, $offset,$count) = @_ ;
    my $prefix = substr($sequence, 0, $offset);
    my $postfix = substr($sequence, $offset);
    my $new_sequence  = $prefix.("-" x ($count).$postfix);
    return($new_sequence);
}


sub sum_of_pairs {

    my $sub_matrix_file = shift @_;

    my $array_name = shift @_;

    my @array = @{$array_name};

    open(MX, "$sub_matrix_file");

    my @temp_matrix = <MX>;

    close MX;

    ## read substitution matrix

    my %sub_score = read_matrix_file(@temp_matrix);

    my $gap = "-";

    my $gap_penalty = -1;

   $sub_score{$gap}{$gap}=0;

    foreach my $aa (keys %sub_score){
	$sub_score{$aa}{$gap} = $gap_penalty;
	$sub_score{$gap}{$aa} = $gap_penalty;
    }

    $sub_score{$gap}{$gap}=0;



    ## comparisons

    my $num_cols = $#{$array[0]} + 1;

    my $sum_total = 0;

    for (my $i = 0; $i < $#array; $i++){
	for (my $j = $i+1 ; $j <= $#array; $j++){
	    for (my $k = 0 ; $k < $num_cols ; $k++){
		$sum_total += $sub_score{$array[$i][$k]}{$array[$j][$k]};
	    }
	}
    }
 
    return($sum_total);
	
}


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

    my $alpha;

    my @matrix_letters = split(' ',$matrix_data[$m]);

    for(my $i=0;$i<=23;$i++)
    {
        chomp $matrix_data[$k];
        my @row = split(/\s+/,$matrix_data[$k]);
        $alpha = shift @row;
        for(my $j=0;$j<=23;$j++){
            $matrix{$alpha}{$matrix_letters[$j]} = $row[$j];
        }
        $k++;
    }
    return(%matrix);

}

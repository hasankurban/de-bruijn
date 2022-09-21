#!/usr/bin/perl

#####
# This is the main program that is called from the html form, and in turn calls all other programs
# It first saves the uploaded files

use CGI qw(:standard);

my $uploaded_files_path = "/var/www/cgi-bin/rpatward/debruijn/writable";
my $program_path = "/var/www/cgi-bin/rpatward/debruijn";
 
my $motif_program = "new_0310.pl";        #o/p: motif.txt, graph_input.txt, weights.txt
my $motif_output = "motif.txt";
my $graph_program = "Read_nodes.pl";
my $search_program = "search_motif.pl";     # uses $motif_output, o/p: output.txt, allmotifs.txt 
my $search_output = "output.txt";        
my $display_program = "display_motif.pl";  # uses $search_output


print header,start_html(-title =>'Results',-bgcolor =>'#DDDD99'), "<font face = \"verdana\" size = 3 color= #669933><center><B>Motif Discovery Results</B></center>";
print "</font>";
print "<font face = \"verdana\" size = 2 ><p>";
&main();
print "</font></input>";    
print end_html;

#### Subroutines 

# the results are printed if there has been parameters submitted

sub main {
    my $file1 = param('file1');
    my $matrix = param('matrices');
    my $node_length = param('node_len');
    my $scale_threshold = param('scale_threshold');	
    my $scale_max_threshold = param('scale_max_threshold');
#    my $nsites = param('nsites');
    my $nsites = 100;
    
    if (!$file1) {
	print "The protein family file not uploaded correctly.";
	return;
    }
    
    ## Read the input files and write the data to file1.txt and file2.txt
       
    my @file1data = &open_file($file1);
    &save_file(@file1data, "file1.txt");
    system("chmod 755 $uploaded_file_path/file1.txt");

    system("perl $program_path/$motif_program -fam $uploaded_files_path/file1.txt -sub $program_path/$matrix.txt -nlen $node_length  -th $scale_threshold -mth $scale_max_threshold  -n $nsites");
  

#      system("perl $program_path/$graph_program");
      
#      system("perl $program_path/$search_program");
    
#      system("perl $program_path/$display_program"); 
      

}


# simply opens a file and returns an array of the data
  sub open_file {
     my $file = shift;
     my @data = ();

     my $line = "";

     while(<$file>) {
        $line = $_;
        chomp($line);
	push(@data, $line);
     }

     return @data;
  }

# saves the file to a writable directory 
  sub save_file {
      my $name = pop;
      my @saved = @_;

      open(OUT, ">$uploaded_files_path/$name");
     
      foreach my $row(@saved) {
	 chomp($row);
	 print OUT "$row\n";
         
      }

      close OUT;
  }



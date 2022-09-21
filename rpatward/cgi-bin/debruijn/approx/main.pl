#!/usr/bin/perl

#####
# This is the main program that is called from the html form, and in turn calls all other programs
# It first saves the uploaded files

use CGI qw(:standard);

my $uploaded_files_path = "/var/www/cgi-bin/rpatward/debruijn/approx/writable";
my $program_path = "/var/www/cgi-bin/rpatward/debruijn/approx";
my $matrix_path = "/var/www/cgi-bin/rpatward/debruijn";

my $random = int(rand()*10000);

my $rpath = "r_".$random;

my $dir_path = $program_path."/$rpath";

mkdir $dir_path, 0777 || die "Server busy !!";

chmod 0777, $dir_path;

#my $motif_program = "sim_batch_testing_3.pl";        #o/p: motif.txt, graph_input.txt, weights.txt
my $motif_program = "backup.pl";        #o/p: motif.txt, graph_input.txt, weights.txt
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
exit;


#### Subroutines 

# the results are printed if there has been parameters submitted

sub main {
    my $file1 = param('file1');
    my $matrix = param('matrices');
    my $sim_constant = param('sim_constant');
    my $scale_threshold = param('scale_threshold');	
    my $scale_max_threshold = param('scale_max_threshold');
    
    open(CNT, ">>counter");
    print CNT "1\t$ENV{REMOTE_ADDR}\n";
    close CNT;

    if (!$file1) {
	print "The protein family file not uploaded correctly.";
	return;
    }
    
    ## Read the input files and write the data to file1.txt and file2.txt
       
    my @file1data = &open_file($file1);
    &save_file(@file1data, "file1.txt");
    system("chmod 755 $uploaded_file_path/file1.txt");

    system("perl $program_path/$motif_program -fam $dir_path/file1.txt -sub $matrix_path/$matrix.txt -sim_constant $sim_constant -th $scale_threshold -mth $scale_max_threshold -dir $dir_path");
  
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

      open(OUT, ">$dir_path/$name");
     
      foreach my $row(@saved) {
	 chomp($row);
	 print OUT "$row\n";
         
      }

      close OUT;
  }



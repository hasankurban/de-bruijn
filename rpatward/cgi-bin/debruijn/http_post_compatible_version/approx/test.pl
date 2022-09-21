#!/usr/bin/perl

use CGI qw(:standard);


#print "Content-type: text/xml\n\n";

my $program_path = "/var/www/cgi-bin/rpatward/debruijn/http_post_compatible_version/approx";
my $matrix_path = "/var/www/cgi-bin/rpatward/debruijn/";

=head
my $random = int(rand()*10000);
my $rpath = "r_".$random;
my $dir_path = $program_path."/$rpath";

mkdir $dir_path, 0777 || die "<error>Server busy</error>";

chmod 0777, $dir_path;

my $uploaded_files_path = $dir_path;

my $motif_program = "backup.pl";
my $graph_program = "Read_nodes.pl";

=cut



#print "<?xml version='1.0' encoding='UTF-8'?>";


=head
## The Internal DTD


print "<!DOCTYPE motifs [";
print "<!ELEMENT motifs (motif+)>";
print "<!ELEMENT motif (protein+)>";
print "<!ELEMENT protein (id,region,index)>";
print "<!ELEMENT id (#CDATA)>";
print "<!ELEMENT region (#CDATA)>";
print "<!ELEMENT index (#CDATA)>";
print "<!ATTLIST motif consensus CDATA #REQUIRED>";
print "]>";

=cut


print header,start_html(-title =>'Results',-bgcolor =>'#DDDD99');

print "<motif>hi</motif>";



print end_html;
exit;


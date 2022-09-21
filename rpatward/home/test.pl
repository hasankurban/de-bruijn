$h{one} = 1; 
$h{two} = 2;

#$a[0]= \%h;
print "hi";
#print $a[];

push ( @a, {%h} );

print $a[0]{one};

exit;

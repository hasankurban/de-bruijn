package Time::Out ;
@ISA = qw(Exporter) ;
@EXPORT = qw(timeout affects) ;

use strict ;
use Exporter ;
use Carp ;


$Time::Out::VERSION = '0.05' ;


sub timeout($@){
    my $secs = shift ;
    carp("Timeout value evaluates to 0: no timeout will be set") if ! $secs ;
    my $code = pop ;
    usage() unless ((defined($code))&&(UNIVERSAL::isa($code, 'CODE'))) ;
    my @other_args = @_ ;

    my @ret = eval {
	local $SIG{ALRM} = sub { die $code } ;
	if (Time::HiRes->can('alarm')){
	    Time::HiRes::alarm($secs) ;
	  }
	else{
	    alarm($secs) ;
	}
	$code->(@other_args) ;
    } ;
    if ($@){
	alarm(0) ;
	if ((ref($@))&&($@ eq $code)){
	    $@ = "timeout" ;
	}
	else {
	    if (! ref($@)){
		chomp($@) ;
		die("$@\n") ;
	    }
	    else {
		croak $@ ;
	    }
	}
    }
    else{
	alarm(0) ;
    }

    return wantarray ? @ret : $ret[0] ;
}


sub affects(&){
    my $code = shift ;
    usage() unless ((defined($code))&&(UNIVERSAL::isa($code, 'CODE'))) ;

    return $code ;
}


sub usage {
    croak("Usage: timeout \$nb_secs => affects {\n  #code\n} ;\n") ;
    }



1 ;

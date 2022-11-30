#!/usr/local/bin/perl -w
use strict;
use Cwd;
use Getopt::Long;
use FindBin qw($Bin $Script);
my %opts;
GetOptions(\%opts,"i=s","o=s","ref=i","h" );
if(! defined($opts{i}) ||! defined($opts{o})||! defined($opts{'ref'}) || defined($opts{h})){
	&USAGE;
}


# open file
open (IN,"$opts{i}") || die "Can't open file $opts{i}\n";
open (OUT,">$opts{o}") || die "open or create file $opts{o} is failed\n";

# get individual info
my ($chr,$pos,$id,$ref,$alt,$qual,$filter,$info,$format) = ();
my $get_indi_ok = 0;
my @indi = ();
my $indi_num = 0;
while(<IN>){
	chomp;
	next if (/^$/);		# ignore the blank lines
	# check individual info
	if(/^\#CHROM/){
		($chr,$pos,$id,$ref,$alt,$qual,$filter,$info,$format,@indi)=(split/\s+/,$_);
		$indi_num = @indi;
		$get_indi_ok = 1;
		last;
	}
	# check commend
	next if(/^\#/); 	# ignore the commend
	# error: get individual info is failed
	print "error1: get individual info is failed\n";
	exit;
}


# check get indi result
if($get_indi_ok==0){
	print "error2: get individual info is failed\n";
	exit;
}


# print header
print OUT "#Chr\tPos\t";
if($opts{'ref'}){print OUT "Ref\tAlt";}
for (my $i=0; $i<@indi; $i++){
	print OUT "\t$indi[$i]\tDepth\tAlleDp" ;
}
print OUT "\tEffect\tCodon_change\tGene_id\n";


# abstract snplist line by line
my @flag = ();
my @twobase = ("","",""); # for sample alleles not the same as ref
while(<IN>){
	chomp;
	next if (/^$/);		# ignore the blank lines
	# abstract info from one line
	($chr,$pos,$id,$twobase[0],$twobase[1],$qual,$filter,$info,$format,@flag)=(split/\s+/,$_);
	# ignore indel
	#next if( (length($twobase[0])!=1) || (length($twobase[1])!=1) );
	#next if (length($twobase[0]) != 1);
	my $raw_alt = $twobase[1] ;
	#if (length($twobase[1]) != 1){
	if ($twobase[1] =~ /\,/){
		my @bases = split /\,/, $twobase[1] ;
		#my $flag = 0 ;
		#for my $base (@bases){
		#	if (length $base != 1){
		#		$flag = 1 ;
		#		last ;
		#	}
		#}
		#next if ($flag == 1) ;
		for (my $i=0; $i<@bases; $i++){
			$twobase[$i+1] = $bases[$i] ;
		}
	}
	# output chr and pos
	print OUT "$chr\t$pos";
	# output ref if "-ref is show"
	if( $opts{'ref'} ){
		print OUT "\t$twobase[0]\t$raw_alt";
	}
	# iter process each indi
	foreach my $f (@flag){
		my ($gt, $ad, $dp) = (split /\:/, $f)[0, 1, 2] ;
		my @bases = split /\//, $gt ;
		#my @ads = split /,/, $ad ;
		my %hbase = ();
		my @new_base = ();
		for (my $i=0; $i<@bases; $i++){
			next if ($bases[$i] eq '.');
			$hbase{$twobase[$bases[$i]]} ++ ;
			push @new_base, $twobase[$bases[$i]] ;
		}
		my @f_base = keys %hbase ;
		# for loss
		if (@f_base == 0){
			print OUT "\tN\t.\t." ;
			next ;
		}
		# for normal base
		#my $t = join ',', sort @f_base ;
		my $t = join ',', @new_base ;
		print OUT "\t$t\t$dp\t$ad";
		#if(!exists $type{$t}){
		#	print "key($t) not exists\n";
		#	exit;
		#}else{
		#	print OUT "\t$type{$t}\t$dp\t$ad";
		#}
	}
	if (m/SNPEFF_EFFECT=(.*?);/){
		print OUT "\t$1" ;
		if (m/SNPEFF_CODON_CHANGE=(\D+?);/){
			print OUT "\t$1" ;
			if (m/SNPEFF_GENE_NAME=(.*?);/){
				#(my $name = $1) =~ s/\.\d+$// ;
				print OUT "\t$1" ;
			}
		}
	}
	else{
		print OUT "\tUNKNOWN" ;
	}
	# print \n
	print OUT "\n";
}

close OUT;
close IN;


###############Time

sub sub_format_datetime {#Time calculation subroutine
    my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	$wday = $yday = $isdst = 0;
    sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}


sub Runtime{ # &Runtime($BEGIN);
	my ($t1)=@_;
	my $t=time()-$t1;
	print "Total elapsed time: ${t}s\n\n";
}




sub USAGE{
print <<"Usage End.";
Usage: perl vcf_to_snplist.pl -i test.vcf -o out.hete -ref 0
	-i   SNP rawdata file(vcf format)		must be given

	-o   the filename of snplist			must be given

	-ref output ref base (0:no; 1:yes)		must be given

	-h   Help document
Usage End.
exit;
}



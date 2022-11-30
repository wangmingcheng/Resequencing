#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
#===============================================================

my %opts;
GetOptions(\%opts,"i=s","l=s","o=s","h" );

#&help()if(defined $opts{h});
if(!defined($opts{i}) || !defined($opts{l}) || !defined($opts{o}) || defined($opts{h}))
{
	print <<"	Usage End.";

	Usage:

		-i           depth file               <infile>               must be given

		-l           ref length file          <infile>               must be given

		-o           outfile prefix           <outfile>              must be given

		-h           Help document

	Usage End.

	exit;
}

###############Time
my $Time_Start;
$Time_Start = sub_format_datetime(localtime(time()));
print "\nStart Time :[$Time_Start]\n\n";
################

# get parameters
my $infile  = $opts{i} ;
my $reflenfile = $opts{l} ;
my $outfile = $opts{o} ;

# get ref length
my $reflen = &get_ref_length($reflenfile);

# stat depth and coverage
&stat_depth_and_coverage($infile, $reflen, $outfile);

###############Time
my $Time_End;
$Time_End = sub_format_datetime(localtime(time()));
print "\nEnd Time :[$Time_End]\n\n";

###############Subs
sub sub_format_datetime {#Time calculation subroutine
    my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	$wday = $yday = $isdst = 0;
    sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

sub ABSOLUTE_DIR
{ #$pavfile=&ABSOLUTE_DIR($pavfile);
	my $cur_dir=`pwd`;
	$cur_dir =~ s/\n$//;
	my ($in)=@_;
	my $return="";
	
	if(-f $in)
	{
		my $dir=dirname($in);
		my $file=basename($in);
		chdir $dir;$dir=`pwd`;
		$dir =~ s/\n$// ;
		$return="$dir/$file";
	}
	elsif(-d $in)
	{
		chdir $in;$return=`pwd`;
		$return =~ s/\n$// ;
	}
	else
	{
		warn "Warning just for file and dir\n";
		exit;
	}
	
	chdir $cur_dir;
	return $return;
}

# &show_log("txt")
sub show_log()
{
	my ($txt) = @_ ;
	my $time = time();
	my $Time = &sub_format_datetime(localtime($time));
	print "$Time:\t$txt\n" ;
	return ($time) ;
}

#&run_or_die($cmd);
sub run_or_die()
{
	my ($cmd) = @_ ;
	&show_log($cmd);
	my $flag = system($cmd) ;
	if ($flag != 0){
		&show_log("Error: command fail: $cmd");
		exit(1);
	}
	&show_log("done.");
	return ;
}

#my $reflen = &get_ref_length($reflenfile);
sub get_ref_length()
{
	my ($reflenfile) = @_ ;
	open (IN, $reflenfile) || die "Can't open $reflenfile, $!\n" ;
	my $total_base ;
	while (<IN>){
		chomp ;
		if (m/^\#Total_Len/){
			$total_base = (split /\s+/)[1] ;
			last ;
		}
	}
	close(IN);
	return ($total_base);
}

#&stat_depth_and_coverage($infile, $reflen, $outfile);
sub stat_depth_and_coverage()
{
	my ($infile, $reflen, $outfile) = @_ ;
	my @base_depth = () ;
	# record base depth
	my ($cov_base, $total_depth) = &reading_depth_file(\@base_depth, $infile);

	# stat ave depth and coverage
	my $cov_ratio = $cov_base/$reflen ;
    print "#######################$cov_base\t$reflen\t$cov_ratio\n";
	$base_depth[0] = $reflen - $cov_base ;  # unmapped base (with coverage 0) count by reflen - covered_base
	my @cumulatively_base_depth = ();
	&get_cumulatively_base_depth(\@base_depth, \@cumulatively_base_depth) ;
	my ($ave_depth) = &change_depth_to_percent(\@base_depth, \@cumulatively_base_depth, $total_depth, $cov_base, $reflen) ;

	# output
	open (ST, ">$outfile.stat") || die "Can't creat $outfile.stat, $!\n" ;
	print ST "ave_depth:\t$ave_depth\n" ;
	print ST "cov_ratio_1X:\t$cov_ratio\n" ;
	print ST "cov_ratio_5X:\t$cumulatively_base_depth[5]\n" ;
	print ST "cov_ratio_10X:\t$cumulatively_base_depth[10]\n" ;
	close(ST) ;
	open (OUT, ">$outfile.dis") || die "Can't creat $outfile.dis, $!\n" ;
	print OUT "#depth\tbase_frac\tcumulatively_base_frac\n" ;
	for (my $i=0; $i<@base_depth; $i++){
		print OUT "$i\t", $base_depth[$i],"\t", $cumulatively_base_depth[$i],"\n" ;
	}
	close(OUT);
	open (CUT, ">$outfile.dis.cut") || die "Can't creat $outfile.dis.cut, $!\n" ;
	print CUT "#depth\tbase_frac\tcumulatively_base_frac\n" ;
	my $end_number = &select_end_number($ave_depth);
	for (my $i=1; $i<=$end_number && $i<@base_depth; $i++){
		print CUT "$i\t", $base_depth[$i],"\t", $cumulatively_base_depth[$i],"\n" ;
	}
	close(CUT);

	return ;
}

#my $end_number = &select_end_number($ave_depth);
sub select_end_number()
{
	my ($ave_depth) = @_ ;
	my $end_number ;
	if ($ave_depth < 7){
		$end_number = 30 ;
	}
	elsif ($ave_depth < 15){
		$end_number = 50 ;
	}
	elsif ($ave_depth < 30){
		$end_number = 70 ;
	}
	elsif ($ave_depth < 40){
		$end_number = 90 ;
	}
	elsif ($ave_depth < 80){
		$end_number = 160 ;
	}
	elsif ($ave_depth < 110){
		$end_number = 200 ;
	}
	else {
		$end_number = 300 ;
	}

	return ($end_number);
}

#&reading_depth_file()
sub reading_depth_file()
{
	my ($abase_depth, $infile) = @_ ;
	my $cov_base = 0 ;
	my $total_depth = 0 ;
	open (IN, $infile) || die "Can't open $infile, $!\n" ;
	while(<IN>){
		chomp ;
		next if (m/^\#/ || m/^\s*$/);
		my ($chr, $pos, $depth) = split ;
		$total_depth += $depth ;
		$cov_base++ if ($depth > 0) ;
		if ($depth > 1000){
			$abase_depth->[1000] ++ ;
		}
		else{
			$abase_depth->[$depth] ++ ;
		}
	}
	close(IN);

	return($cov_base, $total_depth) ;
}

#&get_cumulatively_base_depth(\@base_depth, \@cumulatively_base_depth) ;
sub get_cumulatively_base_depth()
{
	my ($abase_depth, $acumulatively_base_depth) = @_ ;
	my $len = @{$abase_depth} ;
	$acumulatively_base_depth->[$len - 1] = $abase_depth->[-1] ;
	for (my $i=@{$abase_depth}-2; $i>=0; $i--){
		$acumulatively_base_depth->[$i] = $abase_depth->[$i]+$acumulatively_base_depth->[$i+1] ;
	}

	return ;
}

#my ($ave_depth) = &change_depth_to_percent(\@base_depth, \@cumulatively_base_depth, $total_depth, $cov_base, $reflen) ;
sub change_depth_to_percent()
{
	my ($abase_depth, $acumulatively_base_depth, $total_depth, $cov_base, $reflen) = @_ ;
	my $ave_depth ;
	my $cum_depth = 0 ;
	for (my $i=1; $i<@{$abase_depth}; $i++){
		$cum_depth += $abase_depth->[$i] ;
		if ($cum_depth > $cov_base/2){
			$ave_depth = $i ;
			last ;
		}
	}
	for (my $i=0; $i<@{$abase_depth}; $i++){
		$abase_depth->[$i] /= $reflen ;
		$acumulatively_base_depth->[$i] /= $reflen ;
	}

	return($ave_depth);
}


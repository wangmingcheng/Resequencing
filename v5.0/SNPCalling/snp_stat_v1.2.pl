#!/usr/bin/perl -w
# 
use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);

my %opts;
GetOptions(\%opts,"i=s","o=s","h" );

if(!defined($opts{i}) || !defined($opts{o}) || defined($opts{h}))
{
	print <<"	Usage End.";

	Usage:

		-i           snp list file                 <infile>          must be given

		-o           snp stat result               <outfile>         must be given

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
my $outfile = $opts{o} ;
my %htrans = (
	'AG' => "Transition", 'CT' => "Transition", 'GA'=>"Transition", 'TC'=>"Transition",
	'AC' => "Transversion", 'AT' => "Transversion", 'CA'=>"Transversion", 'TA'=>"Transversion",
	'CG' => "Transversion", 'GT' => "Transversion", 'GC'=>"Transversion", 'TG'=>"Transversion",
	'R' => "Transition", 'Y'=>'Transition', 
	'K' => "Transversion", 'K'=>"Transversion",
	'M' => "Transversion", 'S'=>'Transversion', 'W'=>"Transversion"
);
# reading snp file
my %hsnp = ();
my ($total_snp, $total_ti, $total_tv, $asamples) = &reading_snp_file_and_stat($infile, \%hsnp) ;

# output stat result
&output_stat_result($total_snp, $total_ti, $total_tv, \%hsnp, $asamples, $outfile);

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


#&reading_snp_file_and_stat($infile, \%hsnp) ;
sub reading_snp_file_and_stat()
{
	my ($infile, $ahsnp) = @_ ;
	open (IN, $infile) || die "Can't open $infile, $!\n" ;
	my $header = <IN> ;
	chomp($header) ;
	my ($chr, $pos, $ref, @sample) = split /\s+/, $header ;
	my $total_snp = 0 ;
	my $total_ti = 0 ;
	my $total_tv = 0 ;
	while (<IN>){
		chomp ;
		my ($chr, $pos, $ref_snp, @snp) = split ;
		$total_snp++ ;
		my $is_ti = &stat_snp($ref_snp, \@snp, \@sample, $ahsnp) ;
		if ($is_ti){
			$total_ti++ ;
		}
		else{
			$total_tv++ ;
		}
	}
	close(IN);

	return ($total_snp, $total_ti, $total_tv, \@sample);
}

#my $is_ti = &stat_snp($ref_snp, \@snp, \@sample, $ahsnp) ;
sub stat_snp()
{
	my ($ref_snp, $asnp, $asample, $ahsnp) = @_ ;
	my $is_ti = 0 ;
	for(my $i=0; $i<@{$asnp}; $i++){
		next if ($asnp->[$i] eq 'N' || $asnp->[$i] eq $ref_snp);
		$ahsnp->{$asample->[$i]}->{"SNPnum"}++ ;
		if ($asnp->[$i] =~ /[ACGT]/){
			$ahsnp->{$asample->[$i]}->{"Homo"}++ ;
			my $trans_type = $htrans{$ref_snp.$asnp->[$i]} ;
			if ($trans_type eq "Transition"){
				$is_ti = 1 ;
			}
			$ahsnp->{$asample->[$i]}->{$trans_type}++ ;
		}
		elsif ($asnp->[$i] =~ /[BDHV]/){
			$ahsnp->{$asample->[$i]}->{"multiploidy"} ++ ;
		}
		else{
			$ahsnp->{$asample->[$i]}->{"Heter"}++ ;
			my $trans_type = $htrans{$asnp->[$i]} ;
			if ($trans_type eq "Transition"){
				$is_ti = 1 ;
			}
			$ahsnp->{$asample->[$i]}->{$trans_type}++ ;
		}
	}

	return($is_ti) ;
}

#&output_stat_result($total_snp, $total_ti, $total_tv, \%hsnp, $asamples, $outfile);
sub output_stat_result()
{
	my ($total_snp, $total_ti, $total_tv, $ahsnp, $asamples, $outfile) = @_ ;
	open (OUT, ">$outfile") || die "Can't creat $outfile, $!\n" ;
	print OUT "BMK ID\tSNPnumber\tTransition\tTransversion\tTi/Tv\tHeterozygosity\tHomozygosity\tHet-ratio\tmultiploidy\n" ;
	for my $sample (@{$asamples}){
		my $snp_num = defined $ahsnp->{$sample}->{"SNPnum"} ? $ahsnp->{$sample}->{"SNPnum"} : 0 ;
		my $transition = defined $ahsnp->{$sample}->{"Transition"} ? $ahsnp->{$sample}->{"Transition"} : 0 ;
		my $transversion = defined $ahsnp->{$sample}->{"Transversion"} ? $ahsnp->{$sample}->{"Transversion"} : 0 ;
		my $Ti_Tv = ($transversion == 0) ? "NA" : (int(100*$transition/$transversion))/100 ;
		my $homo_num = defined $ahsnp->{$sample}->{"Homo"} ? $ahsnp->{$sample}->{"Homo"} : 0 ;
		my $heter_num = defined $ahsnp->{$sample}->{"Heter"} ? $ahsnp->{$sample}->{"Heter"} : 0 ;
		my $heter_ratio = (int(10000*$heter_num/$snp_num))/100 ;
		my $multiploidy = defined $ahsnp->{$sample}->{"multiploidy"} ? $ahsnp->{$sample}->{"multiploidy"} : 0 ;
		print OUT $sample,"\t", $snp_num,"\t", $transition,"\t", $transversion,"\t", $Ti_Tv,"\t", $heter_num,"\t", $homo_num,"\t", $heter_ratio,"\%\t",$multiploidy,"\n" ; 
	}
	my $total_titv = int(100*($total_ti/$total_tv))/100 ;
	print OUT "Total\t$total_snp\t$total_ti\t$total_tv\t$total_titv\n" ;
	close(OUT);

	return ;
}



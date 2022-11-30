#!/usr/bin/perl -w

use strict;
use Cwd;
use Getopt::Long;
use Data::Dumper;
use File::Basename qw(basename dirname);
use FindBin qw($Bin $Script);

my $programe = &ABSOLUTE_DIR($0);
my $path     = dirname($0);

#######################################################################################

my ($in, $od, $key, $chr_id, $extra_fa,$link);
GetOptions(
			"in:s"	=>	\$in,
			"od:s"	=>	\$od,
			"of:s"	=>	\$key,
			"chr:s" =>	\$chr_id,
			"ex:s"  =>      \$extra_fa,
			"link:s" =>     \$link,
			"h|?"	=>	\&help,
			) || &help;
&help unless ($in && $od && $key && $chr_id);

sub help
{
	print <<"	Usage End.";
    Description:
    Usage:
        -in         infile                    must be given
        -od         outdir                    must be given
        -of         out fa file               must be given
        -chr        chr_id table              optional
	-ex         extra fasta file          optional
	-link       link or unlink scaffold   optional
        -h          Help document
	Usage End.
	exit;
}

#######################################################################################
mkdir $od if (!-d $od);
$od = &ABSOLUTE_DIR($od);
$in = &ABSOLUTE_DIR($in);
if(defined $extra_fa){
	$extra_fa = &ABSOLUTE_DIR($extra_fa);
}

##################################################################add by wangyt in 2019.3.11
open IN,$in or die $!;
my %link_hash;
my $chr;
while(<IN>){
	chomp;
	next if /^$/;
	if(/^>(\S+)/){
		$chr = $1;	
	}
	$link_hash{$chr} = 1;
}
close IN;

my %chr_hash;
$chr_id = &ABSOLUTE_DIR($chr_id);
open(IN, $chr_id) || die $!;
while (<IN>) {
	chomp;
	next if(/^$/ || /^\#/);
	my ($seq_id, $others) = split(/\s+/, $_, 2);
	$chr_hash{$seq_id} = $others;
	delete $link_hash{$seq_id};
}
close(IN);


###################

open(IN, $in) || die $!;
open(OUT, ">$od/$key");
open(OUT2, ">$od/id_table.txt");
$/ = ">";
my $m = 1;
my $outhead;
my $outline;
my $start = 0; 
my %scaffold;
my $new_scaffold_len = 0;
my $total_temp_len = 0;
my $scaffoldnum = 0;
while (<IN>) {
	chomp;
	next if (/^$/);
	my ($id, $seq) = split(/\n/, $_, 2);
	$id = (split(/\s+/, $id))[0];
	$seq =~ s/\n//g;
	delete $link_hash{$id} if exists $link_hash{$id};
	if (exists $chr_hash{$id}){
		my $len = length($seq);
		$seq=~s/(.{100})/$1\n/g; 
		$outhead = ">$id";
		$seq = "$seq\n";
		$seq =~ s/\n\n/\n/g;
		print OUT ">$id\n";
		print OUT "$seq";
		print OUT2 "$id\t$id\t0\t$len\n";
	} elsif($link eq 'link') {
		my $len = length($seq);
		$total_temp_len += $len;
		if(!$scaffold{"Scaffold\_$m"}){
			$scaffold{"Scaffold\_$m"}++;
			$start = 0;
		}
		if($total_temp_len <= 80000000){
			print OUT2 "Scaffold\_$m\t$id\t$start\t$len\n";
			$start = $start + $len + 2000;
			if (eof) {
				$outhead = ">Scaffold\_$m";
				$outline .= "$seq";
				$outline=~s/(.{100})/$1\n/g;
				$outline=~"$outline\n";
				$outline=~s/\n\n/\n/g;
				print OUT "$outhead\n$outline";
			}elsif (!eof and !%link_hash ){
				$outhead = ">Scaffold\_$m";
				$outline .= "$seq";
				$outline=~s/(.{100})/$1\n/g;
				$outline=~"$outline\n";
				$outline=~s/\n\n/\n/g;
				print OUT "$outhead\n$outline\n";		
			}else {
				$outhead = ">Scaffold\_$m";
				$outline .= "$seq"."N" x 2000;
				$total_temp_len += 2000;
			}
		}elsif($total_temp_len > 80000000){
			print OUT2 "Scaffold\_$m\t$id\t$start\t$len\n";
			$outhead = ">Scaffold\_$m";
			$outline .= "$seq";
			$outline=~s/(.{100})/$1\n/g;
			$outline="$outline\n";
			$outline=~s/\n\n/\n/g;
			print OUT "$outhead\n$outline";
			$total_temp_len = 0;
			$outhead = "";
			$outline = "";
			$m++;
		}
			
	}elsif($link eq 'unlink'){
		my $len = length($seq);
		$seq =~s/(.{100})/$1\n/g;
		$outhead = ">$id";
		$outline = "$seq\n";
		$outline =~ s/\n\n/\n/g;
		print OUT ">$id\n";
		print OUT "$outline";
		print OUT2 "$id\t$id\t0\t$len\n";
	}
}
close(IN);
$/ = "\n";
close(OUT);
close(OUT2);

#######################################################################################

sub sub_format_datetime #Time calculation subroutine
{
	my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	$wday = $yday = $isdst = 0;
	sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

sub Runtime # &Runtime($BEGIN);
{
	my ($t1)=@_;
	my $t=time()-$t1;
	print "Total $programe elapsed time : [",&sub_time($t),"]\n";
}

sub sub_time
{
	my ($T)=@_;chomp $T;
	my $s=0;my $m=0;my $h=0;
	if ($T>=3600) {
		my $h=int ($T/3600);
		my $a=$T%3600;
		if ($a>=60) {
			my $m=int($a/60);
			$s=$a%60;
			$T=$h."h\-".$m."m\-".$s."s";
		}else{
			$T=$h."h-"."0m\-".$a."s";
		}
	}else{
		if ($T>=60) {
			my $m=int($T/60);
			$s=$T%60;
			$T=$m."m\-".$s."s";
		}else{
			$T=$T."s";
		}
	}
	return ($T);
}

sub ABSOLUTE_DIR #$pavfile=&ABSOLUTE_DIR($pavfile);
{
	my $cur_dir=`pwd`;chomp($cur_dir);
	my ($in)=@_;
	my $return="";
	if(-f $in){
		my $dir=dirname($in);
		my $file=basename($in);
		chdir $dir;$dir=`pwd`;chomp $dir;
		$return="$dir/$file";
	}elsif(-d $in){
		chdir $in;$return=`pwd`;chomp $return;
	}else{
		warn "Warning just for file and dir\n";
		exit;
	}
	chdir $cur_dir;
	return $return;
}


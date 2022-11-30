#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);

my %opts;
GetOptions(\%opts,"i=s","o=s","s=s","h","c=s","m=s" );

if(!defined($opts{i}) || !defined($opts{o}) || !defined($opts{c}) || defined($opts{h}))
{
	print <<"	Usage End.";

	Usage:

		-i           infile                          must be given

		-o           outfile                         must be given

		-c           chromosome id file              must be given

		-s           filter score, default 20        optional

		-h           Help document

	Usage End.

	exit;
}

## get parameters
my $infile  = $opts{i} ;
my $outfile = $opts{o} ;
my $score = defined $opts{s} ? $opts{s} : 20 ;
my $chr_id = $opts{c} ;
$chr_id = &ABSOLUTE_DIR($chr_id);
my %chr_hash;
open(IN, $chr_id);
while(<IN>){
	chomp;
	next if (/^$/);
	next if (/^\#/);
	my ($chr, $others) = split(/\s+/, $_, 2);
	$chr_hash{$chr} = $others;
}
close(IN);
## filter SV
open (IN,"$infile") || die "Can't open $infile,$!\n" ;
open (OUT,">$outfile") || die "Can't creat $outfile,$!\n" ;
while (<IN>){
	chomp ;
	next if (m/^\s*$/) ;
	if (m/^\#/){
		print OUT $_, "\n" ;
		next ;
	}
	my ($chr1,$chr2,$type,$sv_score) = (split)[0,3,6,8] ;
	next if (!exists $chr_hash{$chr1} || !exists $chr_hash{$chr2});
	next if ($sv_score < $score) || ( $chr1 ne $chr2 && ($type eq "INS" || $type eq "DEL" ||$type eq "ITX" || $type eq "INV") );
	next if ($chr1 eq $chr2 and $type eq "CTX");
	print OUT $_, "\n" ;
}
close(IN);
close(OUT);

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
		warn "Warning just for file and dir [$in]\n";
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
	my $start_time = &show_log($cmd);
	my $flag = system($cmd) ;
	if ($flag != 0){
		my $end_time = &show_log("Error: command fail: $cmd");
		&show_log("Elaseped time: ".($end_time - $start_time)."s\n");
		exit(1);
	}
	my $end_time = &show_log("done.");
	&show_log("Elaseped time: ".($end_time - $start_time)."s\n");
	return ;
}



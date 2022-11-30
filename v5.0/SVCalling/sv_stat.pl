#!/usr/bin/perl -w
# 
use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);

my %opts;
GetOptions(\%opts,"id=s","i=s","o=s","m=s","h" );

if((!defined($opts{i}) && !defined($opts{id})) || !defined($opts{o}) || defined($opts{h}))
{
	print <<"	Usage End.";

	Usage:

		-id          indir of sv files         <indir>          -id/-i must have one
		-i           infile                    <infile>         -id/-i must have one
		-o           outfile                   <outfile>        must be given
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
my $indir = defined $opts{id} ? $opts{id} : "" ;
my $infile = defined $opts{i} ? $opts{i} : "" ;
my $outfile = $opts{o} ;
my $maxproc = defined $opts{m} ? $opts{m} : 20 ;
# reading infile and stat
my %hsv = ();
if (defined $opts{id}){
	$indir = &ABSOLUTE_DIR($indir);
	my @infiles = glob("$indir/*.max") ;
	for my $infile (@infiles){
		&reading_break_dancer_file($infile, \%hsv);
	}
}
if (defined $opts{i}){
	&reading_break_dancer_file($infile, \%hsv);
}

# output stat result
&output_stat_result($outfile, \%hsv) ;




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


#&reading_break_dancer_file($infile, \%hsv);
sub reading_break_dancer_file()
{
	my ($infile, $ahsv) = @_ ;
	open (IN, $infile) || die "Can't open $infile, $!\n" ;
	my $header = <IN> ;
	chomp($header);
	my $sample ;
	$sample = $1 if ($header=~/library:(.*?)L\d+(\s+|$)/ || $header =~ /library:(.*?)(\s+|$)/) ;
	while(<IN>){
		chomp ;
		next if (m/^\#/) ;
		my ($chr1,$pos1,$ori1,$chr2,$pos2,$ori2,$type,$size,$score,$num_reads,$num_read_lib,$alle_freq) = split ;
		$ahsv->{$sample}->{$type}++ ;
		$ahsv->{$sample}->{'all'}++ ;
	}
	close(IN);

	return ;
}

#&output_stat_result($outfile, \%hsv) ;
sub output_stat_result()
{
	my ($outfile, $ahsv) = @_ ;
	open (OUT, ">$outfile") || die "Can't creat $outfile, $!\n" ;
	print OUT "BMK ID\tSV\tINS\tDEL\tINV\tITX\tCTX\tUN\n" ;
	for my $sample (sort keys %{$ahsv}){
		my $sv = defined $ahsv->{$sample}{'all'} ? $ahsv->{$sample}{'all'} : 0 ;
		my $INS = defined $ahsv->{$sample}{'INS'} ? $ahsv->{$sample}{'INS'} : 0 ;
		my $DEL = defined $ahsv->{$sample}{'DEL'} ? $ahsv->{$sample}{'DEL'} : 0 ;
		my $INV = defined $ahsv->{$sample}{'INV'} ? $ahsv->{$sample}{'INV'} : 0 ;
		my $ITX = defined $ahsv->{$sample}{'ITX'} ? $ahsv->{$sample}{'ITX'} : 0 ;
		my $CTX = defined $ahsv->{$sample}{'CTX'} ? $ahsv->{$sample}{'CTX'} : 0 ;
		my $UN = defined $ahsv->{$sample}{'UN'} ? $ahsv->{$sample}{'UN'} : 0 ;
		print OUT join("\t", ($sample, $sv, $INS, $DEL, $INV, $ITX, $CTX, $UN)), "\n" ;
	}
	close(OUT);

	return ;
}



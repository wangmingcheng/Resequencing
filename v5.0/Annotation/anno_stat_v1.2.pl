#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Cwd qw(abs_path);

my %opts;
GetOptions(\%opts,"i=s","id=s","o=s","m=s","h" );

if((!defined($opts{i}) && !defined($opts{id})) || !defined($opts{o}) || defined($opts{h}) and !defined($opts{m}))
{
	print <<"	Usage End.";

	Usage:

		-id          indir of annotation file        <indir>            -id/-i must have one
		-i           annotation file of GATK         <infile>           -id/-i must have one
		-o           statistic result                <outfile>          must be given
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
my$mode=$opts{m};
my ($infile, $indir) ;
if (defined $opts{i}){
	$infile = $opts{i} ;
}
if (defined $opts{id}){
	$indir = $opts{id} ;
	$indir = abs_path($indir);
}
my $outfile = $opts{o} ;
my %hanno = ();
my %htype = ();
my %hinfo = (
	'NONE'=>'--', 'CHROMOSOME'=>'--', 'CUSTOM'=>'--', 'CDS'=>'CDS', 'INTERGENIC'=>'--',
	'INTERGENIC_CONSERVED'=>'--', 'UPSTREAM'=>'--', 'UTR_5_PRIME'=>'--', 'UTR_5_DELETED'=>'--',
	'START_GAINED'=>'--', 'SPLICE_SITE_ACCEPTOR'=>'--', 'SPLICE_SITE_DONOR'=>'--', 'SPLICE_SITE_REGION'=>'--',
	'INTRAGENIC'=>'--', 'START_LOST'=>'CDS', 'SYNONYMOUS_START'=>'CDS', 'NON_SYNONYMOUS_START'=>'CDS',
	'GENE'=>'--', 'TRANSCRIPT'=>'--', 'EXON'=>'CDS', 'EXON_DELETED'=>'CDS', 'NON_SYNONYMOUS_CODING'=>'CDS',
	'SYNONYMOUS_CODING'=>'CDS', 'FRAME_SHIFT'=>'CDS', 'CODON_CHANGE'=>'CDS', 'CODON_INSERTION'=>'CDS',
	'CODON_CHANGE_PLUS_CODON_INSERTION'=>'CDS', 'CODON_DELETION'=>'CDS', 'CODON_CHANGE_PLUS_CODON_DELETION'=>'CDS',
	'STOP_GAINED'=>'CDS', 'SYNONYMOUS_STOP'=>'CDS', 'STOP_LOST'=>'CDS', 'RARE_AMINO_ACID'=>'CDS', 'INTRON'=>'--',
	'INTRON_CONSERVED'=>'--', 'UTR_3_PRIME'=>'--', 'UTR_3_DELETED'=>'--', 'DOWNSTREAM'=>'--',
	'REGULATION'=>'--', 'Other'=>'--'
);
my %horder = (
	'INTERGENIC'=>'1', 'INTERGENIC_CONSERVED'=>'2', 'INTRAGENIC'=>'3', 'INTRON'=>'4', 'INTRON_CONSERVED'=>'5', 
	'UPSTREAM'=>'6', 'DOWNSTREAM'=>'7', 'UTR_5_PRIME'=>'8', 'UTR_5_DELETED'=>'9','UTR_3_PRIME'=>'10', 
	'UTR_3_DELETED'=>'11', 'SPLICE_SITE_ACCEPTOR'=>'12', 'SPLICE_SITE_DONOR'=>'13', 'SPLICE_SITE_REGION'=>'14', 
	'START_GAINED'=>'15', 'START_LOST'=>'16', 'SYNONYMOUS_START'=>'17', 'NON_SYNONYMOUS_START'=>'18',
	'SYNONYMOUS_CODING'=>'19', 'NON_SYNONYMOUS_CODING'=>'20', 'SYNONYMOUS_STOP'=>'21','STOP_GAINED'=>'22', 
	'STOP_LOST'=>'23', 'FRAME_SHIFT'=>'17', 'CODON_CHANGE'=>'18', 'CODON_INSERTION'=>'18', 'CODON_DELETION'=>'18',
	'EXON_DELETED'=>'18', 'CODON_CHANGE_PLUS_CODON_INSERTION'=>'19', 'CODON_CHANGE_PLUS_CODON_DELETION'=>'19', 'RARE_AMINO_ACID'=>'20',
	'Other'=>'24', 'GENE'=>'25', 'TRANSCRIPT'=>'26', 'EXON'=>'27', 
	'REGULATION'=>'28', 'NONE'=>'29', 'CHROMOSOME'=>'30', 'CUSTOM'=>'31', 'CDS'=>'32'
);
#my @region_type = ('--', 'CDS');
my%trans=(
	'TG'=>"T:A>G:C",
	'AC'=>"T:A>G:C",
	'TC'=>"T:A>C:G",
	'AG'=>"T:A>C:G",
	'TA'=>"T:A>A:T",
	'AT'=>"T:A>A:T",
	'CT'=>"C:G>T:A",
	'GA'=>"C:G>T:A",
	'CG'=>"C:G>G:C",
	'GC'=>"C:G>G:C",
	'CA'=>"C:G>A:T",
	'GT'=>"C:G>A:T"
);
my$outdir=dirname($outfile);
# reading vcf file
if (defined $opts{id}){
	my @vcffiles = glob("$indir/*.vcf");
	if ($mode=~/snp/i){
		
		open QUA ,">$outdir/$mode.quality.stat.xls" or die "$!";
		print QUA "#Id\tvalue\tsample\n";
		for my $vcf_file(@vcffiles){
			&stat_vcf_snp($vcf_file);
		}
		
		close(QUA);
		print "Rscript $Bin/plot_snp_quality.r -i $outdir/$mode.quality.stat.xls -o $outdir -n SNP";

			system("Rscript $Bin/plot_snp_quality.r -i $outdir/$mode.quality.stat.xls -o $outdir -n SNP");
		
	}

	
	
	if ($mode=~/indel/i){
		#my$outdir=dirname($outfile);
		open QUA ,">$outdir/$mode.vennID.stat.xls" or die "$!";
		print QUA "#Id\tvalue\tsample\n";
		for my $vcf_file(@vcffiles){
			&stat_vcf_indel($vcf_file);
		}
		close(QUA);
		print "Rscript $Bin/plot_indel_venn.r -i $outdir/$mode.vennID.stat.xls -o $outdir -n INDEL\n";

		system("Rscript $Bin/plot_indel_venn.r -i $outdir/$mode.vennID.stat.xls -o $outdir -n INDEL");
		
		
	}

	
	
	
	
	for my $vcf_file(@vcffiles){
		&stat_vcf_file($vcf_file, \%hanno, \%htype);
		
	}
	
}
if (defined $opts{i}){
	&stat_vcf_file($infile, \%hanno, \%htype);
	
}

# output result
&output_result($outfile, \%hanno, \%htype, \%hinfo, \%horder);

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

sub stat_vcf_file()
{
	my ($vcf_file, $ahanno, $ahtype) = @_ ;
	open (IN, $vcf_file) || die "Can't open $vcf_file, $!\n" ;
	my ($chr,$pos,$id,$ref,$alt,$qual,$filter,$info,$format) = ();
	my $get_indi_ok = 0;
	my @indi = ();
	my $indi_num = 0;
	my $sample_name = "all" ;
	while(<IN>){
		chomp ;
		next if (m/^\s*$/);
		if (m/^\#CHROM/){
			($chr,$pos,$id,$ref,$alt,$qual,$filter,$info,$format,@indi)=(split/\s+/,$_);
			$indi_num = @indi;
			if ($indi_num == 1){
				$sample_name = $indi[0] ;
			}
		}
		next if (m/^\#/);
		if (m/SNPEFF_EFFECT=(.*?)\;/){
			$ahanno->{$sample_name}->{$1}++ ;
			$ahtype->{$1}++ ;
		}
		else{
			$ahanno->{$sample_name}->{'Other'}++ ;    # un-annotation snp number
		}
	}
	close(IN);
	
	return ;
	
	
	
}




sub stat_vcf_snp(){
	my$vcf=shift;
	open  IN ,"$vcf" or die "$!";
	my$filename=basename($vcf);
	my($sample)=($filename=~/\.(R\d+)\./);
	if ($sample eq ""){die "can't find sample name in file name $vcf,please check R\\d+\n";}
	my$upPos=0;
	my$upChr=0;
	my%tranStat=();
	while(<IN>){
		next if (/^#/);
		my($chr,$pos,$id,$ref,$alt,$qual,$filter,$info,$format,@indi)=(split/\s+/,$_);
		#next if $alt=~/,/;
		my$len=length($alt);
		my($dp)=($info=~/DP=(\d+);/);
		print QUA "id\t$chr-$pos\t$sample\n";
		print QUA "depth\t$dp\t$sample\n";
		if($upChr ne $chr){
			$upPos=0;
			$upChr=$chr;
		}
		
		if($upPos==0){
			$upPos=$pos;
		}else{
			my$gap=$pos-$upPos;
			$upPos=$pos;
			print QUA "distance\t$gap\t$sample\n";
		}
		if( $alt=~/,/){
			
		}else{
			$tranStat{$trans{"$ref$alt"}}++;
		}
		
	}
	
	for my$i(keys %tranStat ){
		print QUA "trans\t$i-$tranStat{$i}\t$sample\n";
	}
	
}
sub stat_vcf_indel(){
	my$vcf=shift;
	open  IN ,"$vcf" or die "$!";
	my$filename=basename($vcf);
	my($sample)=($filename=~/\.(R\d+)\./);
	if ($sample=~/^\s*$/){die "can't find sample name in file name $vcf,please check R\\d\n";}
	my$upPos=0;
	my$upChr=0;
	while(<IN>){
		next if (/^#/);
		my($chr,$pos,$id,$ref,$alt,$qual,$filter,$info,$format,@indi)=(split/\s+/,$_);
		#next if $alt=~/,/;
		my$len=length($alt);
		my($dp)=($info=~/DP=(\d+);/);
		print QUA "id\t$chr-$pos\t$sample\n";

	}
	
	
}
#&output_result($outfile, \%hanno, \%htype, \%hinfo, \%horder);
sub output_result()
{
	my ($outfile, $ahanno, $ahtype, $ahinfo, $ahorder) = @_ ;
	open (OUT, ">$outfile") || die "Can't creat $outfile, $!\n" ;
	my @types = sort {$ahorder->{$a} <=> $ahorder->{$b}} keys %{$ahtype} ;
	push @types, 'Other' ;
	my @samples = sort keys %{$ahanno} ;
	print OUT "#Region\tType\t", join ("\t", @samples), "\n" ;
	for my $type (@types){
		print OUT $ahinfo->{$type}, "\t$type" ;
		for my $sample (@samples){
			my $num = defined $ahanno->{$sample}->{$type} ? $ahanno->{$sample}->{$type} : 0 ;
			print OUT "\t", $num ;
		}
		print OUT "\n" ;
	}
	close(OUT);

	return ;
}



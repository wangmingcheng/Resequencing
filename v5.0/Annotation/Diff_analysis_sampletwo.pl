#!/usr/bin/perl -w
use strict;
use Cwd;
use Getopt::Long;
use Data::Dumper;
use File::Basename qw(basename dirname);
use FindBin qw($Bin $Script);
use Cwd qw(abs_path);
my $writer = "Liangsh<liangsh\@biomarker.com.cn>";
###########################################
my $Time_Start;
$Time_Start = sub_format_datetime(localtime(time()));
print "\nStart Time :[$Time_Start]\n\n";
##############################################
my ($fin,$sam,$od);
GetOptions(
        "h|?"=>\&help,
        "i:s"=>\$fin,
        "od:s"=>\$od,
        "sample:s"=>\$sam,
    ) or &help;
&help unless ($fin && $sam && $od);
mkdir $od unless (-d $od);
$od = ABSOLUTE_DIR($od);
my ($sample1,$sample2) = split /,/,$sam;
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
open IN,"<$fin" or die $!;
my $vcf = "$od/$sample1\_$sample2.snp.anno.gatk.vcf";
open OUT,">$vcf" or die $!;
my ($index1,$index2);
while (<IN>) {
    chomp;
    if (/^##/) {
        print OUT "$_\n";
        next;
    }
    my ($chr,$pos,$id,$ref,$alt,$qual,$filter,$info,$format,@indi)=split/\s+/,$_;
    if (/^#CHROM/) {
        for (my $i=0; $i < @indi; $i++){
            if ($indi[$i] eq $sample1) {
                $index1 = $i;
            }
            if ($indi[$i] eq $sample2) {
                $index2 = $i;
            }
        }
        print OUT "$chr\t$pos\t$id\t$ref\t$alt\t$qual\t$filter\t$info\t$format\t$indi[$index1]\t$indi[$index2]\n";
        next;
    }
    my $geno1 = (split /:/,$indi[$index1])[0];
    my $geno2 = (split /:/,$indi[$index2])[0];
    next if ($geno1 eq $geno2);
    print OUT "$chr\t$pos\t$id\t$ref\t$alt\t$qual\t$filter\t$info\t$format\t$indi[$index1]\t$indi[$index2]\n";
}
close IN;
close OUT;

&stat_vcf_file($vcf, \%hanno, \%htype);
my $outfile = "$od/$sample1\_$sample2.anno.stat.xls";
&output_result($outfile, \%hanno, \%htype, \%hinfo, \%horder);

my $cmd = "python $Bin/plot_ann.py -i $outfile -o $od/ -p SNP.anno";
print "$cmd\n";
`$cmd`;






##-----------------------sub function
#&stat_vcf_file($vcf_file, \%hanno, \%htype);
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
			#if ($indi_num == 1){
			$sample_name = "$indi[0]\_$indi[1]" ;
			#}
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





#########################################3
my $Time_End;
$Time_End = sub_format_datetime(localtime(time()));
print "\nEnd Time :[$Time_End]\n\n";
###########-----------------
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

sub sub_format_datetime {#Time calculation subroutine
    my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	$wday = $yday = $isdst = 0;
    sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

#######################################
sub help{
    my $help =<< "help";
    Writer: $writer
    Description:  diff snp analysis among two samples 
    Usage:
        -i  gatk anno vcf file 
        -od  output dir
        -sample  the name of two sample [exp: R01,R02]
    
    
help
    print $help;
    exit;

}








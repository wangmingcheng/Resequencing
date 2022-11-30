#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
#===============================================================

my %opts;
GetOptions(\%opts,"i=s","o=s","h" );

if(!defined($opts{i}) || !defined($opts{o}) || defined($opts{h}))
{
	print <<"	Usage End.";

	Usage:
		-i           snpEff annotated vcf format file                <infile>     must be given
		-o           selected one most effective-anno result file    <outfile>    must be given
		-h           Help document

	Usage End.

	exit;
}

###############Time
my $Time_Start;
$Time_Start = sub_format_datetime(localtime(time()));
print "\nStart Time :[$Time_Start]\n\n";
################
## get parameters
my $infile  = $opts{i} ;
my $outfile = $opts{o} ;

## reading snpEff annotated file and select
open (IN,"$infile") || die "Can't open $infile,$!\n" ;
open (OUT,">$outfile") || die "Can't creat $outfile,$!\n" ;
while (<IN>){
	chomp ;
	next if (m/^\s*$/);
	if (m/^\#/){
		if (m/^\#CHROM/){
			print OUT <<"INFO.end" ;
##INFO=<ID=SNPEFF_AMINO_ACID_CHANGE,Number=1,Type=String,Description="Old/New amino acid for the highest-impact effect resulting from the current variant (in HGVS style)">
##INFO=<ID=SNPEFF_CODON_CHANGE,Number=1,Type=String,Description="Old/New codon for the highest-impact effect resulting from the current variant">
##INFO=<ID=SNPEFF_EFFECT,Number=1,Type=String,Description="The highest-impact effect resulting from the current variant (or one of the highest-impact effects, if there is a tie)">
##INFO=<ID=SNPEFF_EXON_ID,Number=1,Type=String,Description="Exon ID for the highest-impact effect resulting from the current variant">
##INFO=<ID=SNPEFF_FUNCTIONAL_CLASS,Number=1,Type=String,Description="Functional class of the highest-impact effect resulting from the current variant: [NONE, SILENT, MISSENSE, NONSENSE]">
##INFO=<ID=SNPEFF_GENE_BIOTYPE,Number=1,Type=String,Description="Gene biotype for the highest-impact effect resulting from the current variant">
##INFO=<ID=SNPEFF_GENE_NAME,Number=1,Type=String,Description="Gene name for the highest-impact effect resulting from the current variant">
##INFO=<ID=SNPEFF_IMPACT,Number=1,Type=String,Description="Impact of the highest-impact effect resulting from the current variant [MODIFIER, LOW, MODERATE, HIGH]">
##INFO=<ID=SNPEFF_TRANSCRIPT_ID,Number=1,Type=String,Description="Transcript ID for the highest-impact effect resulting from the current variant">
INFO.end
		}
		print OUT "$_\n" ;
		next ;
	}
	my ($chr, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, $genotype) = split /\s+/, $_, 10 ;
	my @anno = split /\;/, $info ;
	for (my $i=0; $i<@anno; $i++){
		if ($anno[$i] =~ /EFF=(.*)/){
			my $snpeff_anno_info = $1 ;
			my ($flag, $select_anno) = &select_best_anno($snpeff_anno_info);
			if ($flag == 1){
				$anno[$i] = $anno[$i-1] ;
				$anno[$i-1] = $select_anno ;
			}
			else{
				splice(@anno, $i, 1) ;
			}
			last ;
		}
	}
	print OUT "$chr\t$pos\t$id\t$ref\t$alt\t$qual\t$filter\t", join(";", @anno), "\t$format\t$genotype\n" ;
}
close(IN);
close(OUT);


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


sub select_best_anno()
{
	my ($snpeff_anno_info) = @_ ;
	my ($flag, $select_anno) ;
	my @snpeff_annos = split /\,/, $snpeff_anno_info ;
	my %hanno = ();
	for my $anno (@snpeff_annos){
		my ($effect, $anno_info) = ($anno=~/(.*)\((.*)\)/) ;
		my ($impact, $functional_class, $coden_change, $amino_acid_change, $gene_name, $gene_biotype, $gene_coding, $transcript_id, $exon_id, $error_info) = split /\|/, $anno_info ;
		push @{$hanno{$impact}{$effect}}, [$impact, $functional_class, $coden_change, $amino_acid_change, $gene_name, $gene_biotype, $gene_coding, $transcript_id, $exon_id, $error_info] ;
	}

	for my $type ('HIGH', 'MODERATE', 'LOW', 'MODIFIER'){
		if (defined $hanno{$type}){
			($flag, $select_anno) = &get_one_anno($type, $hanno{$type});
			last ;
		}
	}

	return($flag, $select_anno) ;
}

sub get_one_anno()
{
	my ($type, $ahanno_info) = @_ ;
	my ($flag, $select_anno) ;

	if ($type eq 'HIGH' || $type eq 'MODERATE' || $type eq 'LOW'){
		for my $eff (keys %{$ahanno_info}){
			for my $aanno (@{$ahanno_info->{$eff}}){
				my ($impact, $functional_class, $coden_change, $amino_acid_change, $gene_name, $gene_biotype, $gene_coding, $transcript_id, $exon_id, $error_info) = @{$aanno} ;
				$functional_class = 'NONE' if ($functional_class eq '');
				$flag = 1 ;
				$select_anno = "SNPEFF_EFFECT=$eff;SNPEFF_IMPACT=$impact;SNPEFF_FUNCTIONAL_CLASS=$functional_class;" ;
				$select_anno .= "SNPEFF_CODON_CHANGE=$coden_change;" if($coden_change ne '');
				$select_anno .= "SNPEFF_AMINO_ACID_CHANGE=$amino_acid_change;" if ($amino_acid_change ne '');
				$select_anno .= "SNPEFF_GENE_NAME=$gene_name;" if ($gene_name ne '');
				$select_anno .= "SNPEFF_GENE_BIOTYPE=$gene_biotype;" if ($gene_biotype ne '') ;
				$select_anno .= "SNPEFF_TRANSCRIPT_ID=$transcript_id;" if ($transcript_id ne '') ;
				$select_anno .= "SNPEFF_EXON_ID=$exon_id;" if ($exon_id ne '') ;
				$select_anno =~ s/\;$// ;
				last ;
			}
			return($flag, $select_anno);
		}
	}elsif ($type eq 'MODIFIER'){
		for my $eff ('CDS', 'EXON', 'INTRON', 'INTRON_CONSERVED', 'TRANSCRIPT', 'GENE',
			'INTRAGENIC', 'TRANSCRIPT', 'MICRO_RNA', 'UTR_3_PRIME', 'UTR_5_PRIME', 'REGULATION'){
			if (defined $ahanno_info->{$eff}){
				for my $aanno (@{$ahanno_info->{$eff}}){
					my ($impact, $functional_class, $coden_change, $amino_acid_change, $gene_name, $gene_biotype, $gene_coding, $transcript_id, $exon_id, $error_info) = @{$aanno} ;
					$functional_class = 'NONE' if ($functional_class eq '');
					$flag = 1 ;
					$select_anno = "SNPEFF_EFFECT=$eff;SNPEFF_IMPACT=$impact;SNPEFF_FUNCTIONAL_CLASS=$functional_class;" ;
					$select_anno .= "SNPEFF_CODON_CHANGE=$coden_change;" if($coden_change ne '');
					$select_anno .= "SNPEFF_AMINO_ACID_CHANGE=$amino_acid_change;" if ($amino_acid_change ne '');
					$select_anno .= "SNPEFF_GENE_NAME=$gene_name;" if ($gene_name ne '');
					$select_anno .= "SNPEFF_GENE_BIOTYPE=$gene_biotype;" if ($gene_biotype ne '') ;
					$select_anno .= "SNPEFF_TRANSCRIPT_ID=$transcript_id;" if ($transcript_id ne '') ;
					$select_anno .= "SNPEFF_EXON_ID=$exon_id;" if ($exon_id ne '') ;
					$select_anno =~ s/\;$// ;
					last ;
				}
				return($flag, $select_anno);
			}
		}
		my $min_dis = 500000 ;
		my $return_flag = 0 ;
		for my $eff ('UPSTREAM', 'DOWNSTREAM'){
			if (defined $ahanno_info->{$eff}){
				for my $aanno (@{$ahanno_info->{$eff}}){
					my ($impact, $functional_class, $coden_change, $amino_acid_change, $gene_name, $gene_biotype, $gene_coding, $transcript_id, $exon_id, $error_info) = @{$aanno} ;
					if ($coden_change < $min_dis){
						$return_flag = 1 ;
						$min_dis = $coden_change ;
						$functional_class = 'NONE' if ($functional_class eq '');
						$flag = 1 ;
						$select_anno = "SNPEFF_EFFECT=$eff;SNPEFF_IMPACT=$impact;SNPEFF_FUNCTIONAL_CLASS=$functional_class;" ;
						$select_anno .= "SNPEFF_CODON_CHANGE=$coden_change;" if($coden_change ne '');
						$select_anno .= "SNPEFF_AMINO_ACID_CHANGE=$amino_acid_change;" if ($amino_acid_change ne '');
						$select_anno .= "SNPEFF_GENE_NAME=$gene_name;" if ($gene_name ne '');
						$select_anno .= "SNPEFF_GENE_BIOTYPE=$gene_biotype;" if ($gene_biotype ne '') ;
						$select_anno .= "SNPEFF_TRANSCRIPT_ID=$transcript_id;" if ($transcript_id ne '') ;
						$select_anno .= "SNPEFF_EXON_ID=$exon_id;" if ($exon_id ne '') ;
						$select_anno =~ s/\;$// ;
					}
				}
			}
		}
		return($flag, $select_anno) if ($return_flag == 1) ;

		for my $eff ('INTERGENIC', 'INTERGENIC_CONSERVED'){
			if (defined $ahanno_info->{$eff}){
				for my $aanno (@{$ahanno_info->{$eff}}){
					my ($impact, $functional_class, $coden_change, $amino_acid_change, $gene_name, $gene_biotype, $gene_coding, $transcript_id, $exon_id, $error_info) = @{$aanno} ;
					$functional_class = 'NONE' if ($functional_class eq '');
					$flag = 1 ;
					$select_anno = "SNPEFF_EFFECT=$eff;SNPEFF_IMPACT=$impact;SNPEFF_FUNCTIONAL_CLASS=$functional_class;" ;
					$select_anno .= "SNPEFF_CODON_CHANGE=$coden_change;" if($coden_change ne '');
					$select_anno .= "SNPEFF_AMINO_ACID_CHANGE=$amino_acid_change;" if ($amino_acid_change ne '');
					$select_anno .= "SNPEFF_GENE_NAME=$gene_name;" if ($gene_name ne '');
					$select_anno .= "SNPEFF_GENE_BIOTYPE=$gene_biotype;" if ($gene_biotype ne '') ;
					$select_anno .= "SNPEFF_TRANSCRIPT_ID=$transcript_id;" if ($transcript_id ne '') ;
					$select_anno .= "SNPEFF_EXON_ID=$exon_id;" if ($exon_id ne '') ;
					$select_anno =~ s/\;$// ;
					last ;
				}
				return ($flag, $select_anno) ;
			}
		}
	}

}



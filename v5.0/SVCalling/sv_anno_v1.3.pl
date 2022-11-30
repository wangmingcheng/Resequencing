#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Cwd qw(abs_path);
######################

my %opts;
GetOptions(\%opts,"id=s","i=s","g=s","o=s","annotype=s","h" );

if((!defined($opts{i}) && !defined($opts{id})) || !defined($opts{g}) || !defined($opts{o}) || defined($opts{h}))
{
	print <<"	Usage End.";

	Usage:

		-id          indir of sv files(*.max)         <indir>          -id/-i must have one
		-i           infile                           <infile>         -id/-i must have one
		-g           gff file                         <infile>         must be given
		-o           outfile                          <outfile>        must be given
		-annotype    gene/mRNA                        default gene
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
my $gff_file = $opts{g} ;
my $outfile = $opts{o} ;
my $annotype = defined $opts{annotype} ? $opts{annotype} : "gene";
# reading sv result file and its position
my %hsv = ();
if (defined $opts{id}){
	$indir = abs_path($indir);
	my @infiles = glob("$indir/*.max") ;
	for my $infile (@infiles){
		&reading_break_dancer_file($infile, \%hsv);
	}
}
if (defined $opts{i}){
	&reading_break_dancer_file($infile, \%hsv);
}

# reading gff file and store each element
my %hgff = ();
&reading_gff_file($gff_file, \%hgff) ;
# statistic each type and output
&stat_and_output($outfile, \%hsv, \%hgff) ;

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
		if ($type eq 'INS' || $type eq 'DEL' || $type eq 'INV'){
			if ($pos1 < $pos2){
				push @{$ahsv->{$sample}->{$chr1}}, [$pos1, $pos2, $type, $_] ;
			}
			else{
				push @{$ahsv->{$sample}->{$chr1}}, [$pos2, $pos1, $type, $_] ;
			}
		}
	}
	close(IN);

	return ;
}

#&reading_gff_file($gff_file, \%hgff) ;
sub reading_gff_file()
{
	my ($gff_file, $ahgff) = @_ ;
	open (GFF, $gff_file) || die "Can't open $gff_file, $!\n" ;
	my $name;
	while(<GFF>){
		chomp ;
		next if (m/^\#/ || m/^\s*$/) ;
		my ($seqid, $source, $type, $star, $end, $score, $strand, $phase, $attributes) = split /\t/ ;
		if ($type eq $annotype){
			if ($attributes =~ /ID=(.*?)\;/){
				$name = $1 ;
			}elsif($attributes =~ /ID=(.*)/){
				$name = $1;
			}

		}
		next if ($type eq 'gene' and $annotype eq 'mRNA');
		push @{$ahgff->{$type}->{$seqid}}, [$star, $end, $name] ;
	}
	close(GFF);
	return ;
}

#&stat_and_output($outfile, \%hsv, \%hgff) ;
sub stat_and_output()
{
	my ($outfile, $ahsv, $ahgff) = @_ ;
	open (OUT, ">$outfile.info") || die "Can't creat $outfile.info, $!\n" ;
	my %hstat = ();
	for my $sample (sort keys %{$ahsv}){
		for my $chr (sort keys %{$ahsv->{$sample}}){
			my @sv_region = sort {$a->[0]<=>$b->[0] || $a->[1]<=>$b->[1]} @{$ahsv->{$sample}{$chr}} ;
			my $index_cds = 0 ;
			my $index_gene = 0 ;
			my ($type, $name);
			for my $arr (@sv_region){
				my ($star, $end) = ($arr->[0], $arr->[1]) ;
				($index_cds, $index_gene, $type, $name) = &judge_anno_type($index_cds, $index_gene, $chr, $star, $end, $ahgff);
				if ($type ne "intergenic"){
					print OUT $sample,"\t", $arr->[3],"\t", $type,"\t", $name,"\n" ;
				}
				else{
					print OUT $sample,"\t", $arr->[3],"\t", $type,"\n" ;
				}
				$hstat{$sample}{$arr->[2]}{$type}++ ;
			}
		}
	}
	close(OUT);

	open (OUT, ">$outfile.stat") || die "Can't creat $outfile.stat, $!\n" ;
	print OUT "BMK ID\tType\tExon\tIntron\tIntergenic\n" ;
	for my $sample (sort keys %hstat){
		my $total_exon_num = 0 ;
		my $total_intron_num = 0 ;
		my $total_intergenic_num = 0 ;
		for my $type ("DEL", "INS", "INV"){
			my $exon_num = defined $hstat{$sample}{$type}{'exon'} ? $hstat{$sample}{$type}{'exon'} : 0 ;
			$total_exon_num += $exon_num ;
			my $intron_num = defined $hstat{$sample}{$type}{'intron'} ? $hstat{$sample}{$type}{'intron'} : 0 ;
			$total_intron_num += $intron_num ;
			my $intergenic_num = defined $hstat{$sample}{$type}{'intergenic'} ? $hstat{$sample}{$type}{'intergenic'} : 0 ;
			$total_intergenic_num += $intergenic_num ;
			print OUT "$sample\t$type\t$exon_num\t$intron_num\t$intergenic_num\n" ;
		}
		print OUT "$sample\tTotal\t$total_exon_num\t$total_intron_num\t$total_intergenic_num\n" ;
	}
	close(OUT);

	return ;
}

#($index_cds, $index_gene, $type, $name) = &judge_anno_type($index_cds, $index_gene, $chr, $star, $end, $ahgff);
sub judge_anno_type()
{
	my ($index_cds, $index_gene, $chr, $star, $end, $ahgff) = @_ ;
	my $type = "intergenic" ;
	my $gene_name = "" ;
	if (defined $ahgff->{'gene'} && defined $ahgff->{'gene'}{$chr}){
		my ($flag, $index, $name) = &find_overlap($star, $end, $index_gene, $ahgff->{'gene'}->{$chr});
		$index_gene = $index ;
		if ($flag == 1){
			$type = "intron" ;
			$gene_name = $name ;
		}
	}
	elsif (defined $ahgff->{'mRNA'} && defined $ahgff->{'mRNA'}{$chr}){
		my ($flag, $index, $name) = &find_overlap($star, $end, $index_gene, $ahgff->{'mRNA'}->{$chr});
		$index_gene = $index ;
		if ($flag == 1){
			$type = "intron" ;
			$gene_name = $name ;
		}
	}
	if (defined $ahgff->{'CDS'} && defined $ahgff->{'CDS'}{$chr}){
		my ($flag, $index, $name) = &find_overlap($star, $end, $index_cds, $ahgff->{'CDS'}->{$chr});
		$index_cds = $index ;
		if ($flag == 1){
			$type = "exon" ;
			$gene_name = $name ;
		}
	}
	elsif (defined $ahgff->{'exon'} && defined $ahgff->{'exon'}{$chr}){
		my ($flag, $index, $name) = &find_overlap($star, $end, $index_cds, $ahgff->{'exon'}->{$chr});
		$index_cds = $index ;
		if ($flag == 1){
			$type = "exon" ;
			$gene_name = $name ;
		}
	}

	return ($index_cds, $index_gene, $type, $gene_name);
}

#my ($flag, $index, $name) = &find_overlap($star, $end, $index_gene, $ahgff->{'gene'}->{$chr});
sub find_overlap()
{
	my ($star, $end, $index, $aregion) = @_ ;
	my $flag = 0 ;
	my $return_name = "" ;
	for (my $i=$index; $i<@{$aregion}; $i++){
		my ($ref_star, $ref_end, $name) = @{$aregion->[$i]} ;
		if ($star > $ref_end){
			$index = $i ;
			next ;
		}
		elsif($end < $ref_star){
			last ;
		}
		else{
			$flag = 1 ;
			$return_name = $name ;
			if (!defined $name){
				print Dumper @{$aregion->[$i]} ;
				die ;
			}
			last ;
		}
	}

	return ($flag, $index, $return_name);
}



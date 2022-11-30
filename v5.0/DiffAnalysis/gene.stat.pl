#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($cfg);
GetOptions(
	"help|?" =>\&USAGE,
	"cfg:s"=>\$cfg,
) or &USAGE;

&USAGE unless ($cfg);
my %cfg_hash;
my %project_hash;
open(IN,$cfg) or die $!;
while (<IN>){
	chomp;
    next if /^$/ || /^\#/;
	my @data=split/\s+/,$_;
	$cfg_hash{$data[0]}=$data[1] if ($_=~/^analysis_dir/);
	$project_hash{$data[0]}=$data[1] if ($_=~/^Project/);
}
close(IN);

foreach my $pro (keys %project_hash){
	my %hash;
	my $indir="$cfg_hash{analysis_dir}";
	unless (-d "$indir/Diff_analysis"){
		`mkdir $indir/Diff_analysis`;
	}
	if (-d "$cfg_hash{analysis_dir}/Analysis/SV/") {
		my @svinfo=glob "$cfg_hash{analysis_dir}/Analysis/SV/*.info";
		for (my $i=0;$i<@svinfo ;$i++)
		{
			open (SV,$svinfo[$i]) or dir $!;
			while (<SV>)
			{
				chomp;
				next if ($_=~/\#/);
				next if ($_=~/^\s*$/);
				my ($sample,$exon,$gene)=(split/\t/,$_)[0,15,16];
				next if !defined $exon;
				if ($exon=~/^exon/)
				{
					$hash{$sample}{$gene}=$gene;
				}
			}
			close(SV);
		}
	}
	if (-d "$cfg_hash{analysis_dir}/Analysis/Annotation/snp_anno/") 
	{
		my @sam=glob "$cfg_hash{analysis_dir}/Analysis/Annotation/snp_anno/R*";
		for (my $i=0;$i<@sam ;$i++)
		{
			my $bname=basename($sam[$i]);
			my @samlist=glob "$sam[$i]/*.list";
			for (my $j=0;$j<@samlist ;$j++)
			{
				open(DEG,"$samlist[$j]") or die $!;
				my $header=<DEG>;
				my $totalnum=(split/\t/,$header)[1];
				$hash{$bname}{SNP}=$totalnum;
				close(DEG);
			}
		}
	}
	if (-d "$cfg_hash{analysis_dir}/Analysis/Annotation/indel_anno/") 
	{
		my @sam=glob "$cfg_hash{analysis_dir}/Analysis/Annotation/indel_anno/R*";
		for (my $i=0;$i<@sam ;$i++)
		{
			my $bname=basename($sam[$i]);
			my @samlist=glob "$sam[$i]/*.list";
			for (my $j=0;$j<@samlist ;$j++)
			{
				open(DEG,"$samlist[$j]") or die $!;
				my $header=<DEG>;
				my $totalnum=(split/\t/,$header)[1];
				$hash{$bname}{INDEL}=$totalnum;
				close(DEG);
			}
		}
	}
	
	if (-d "$cfg_hash{analysis_dir}/Analysis/SV/") {
		
		open(OUT,">$indir/Diff_analysis/Diff_gene.stat");
		print "$indir/Diff_analysis/Diff_gene.stat\n";
		print OUT "BMK_ID\tGenes with Non-synonymous SNP\tGenes with InDel\tGenes with SV\n";
		foreach my $sample (sort keys %hash)
		{
			my $b=keys %{$hash{$sample}};
			my $a=$b-2;
			chomp($hash{$sample}{SNP});
			chomp($hash{$sample}{INDEL});
			
			print OUT "$sample\t$hash{$sample}{SNP}\t$hash{$sample}{INDEL}\t$a\n";
			
			open(DESV,">$indir/Analysis/SV/$sample.Diff.gene.list");
			print DESV "#total_num\t$a\n";
			foreach my $ge (sort keys %{$hash{$sample}})
			{
				next if ($ge=~/INDEL/);
				next if ($ge=~/SNP/);
				print DESV "$hash{$sample}{$ge}\n";
			}
			close(DESV);
			
		}
		close(OUT);
	
	}else{
		open(OUT,">$indir/Diff_analysis/Diff_gene.stat");
		print "$indir/Diff_analysis/Diff_gene.stat\n";
		print OUT "BMK_ID\tGenes with Non-synonymous SNP\tGenes with InDel\n";
		foreach my $sample (sort keys %hash){
			chomp($hash{$sample}{SNP});
			chomp($hash{$sample}{INDEL});
			print OUT "$sample\t$hash{$sample}{SNP}\t$hash{$sample}{INDEL}\n";		
		}
		close(OUT);
	}


}

#############################################################################################

sub USAGE {#
	my $usage=<<"USAGE";
Usage:

  -cfg <file>  input file,fasta format,forced 
 
  -h         Help

USAGE
	print $usage;
	exit;
}

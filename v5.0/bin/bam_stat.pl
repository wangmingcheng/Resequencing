use strict;
use File::Basename qw(basename dirname);

print "#sample\tgenome_coverage(%)\taverage_depth\tmapping_rate(%)\n";

my $dir=shift;
my %f;
my($genome_size,$all_covbases,$all_covdepth);
my @coverage=<$dir/*coverage>;
for my $coverage(@coverage){
    open I,"$coverage";
    my $name = basename $coverage;
    my $sample=(split/\./,$name)[0];
    while(<I>){
        chomp;
        next if /^#/;
        my ($start, $end, $covbases, $meandepth)=(split)[1,2,4,6];
        $genome_size+=$end-$start+1;
        $all_covbases+=$covbases;
        $all_covdepth+=$covbases*$meandepth;
    }
    close I;
    $f{$sample}{"coverage"}=$all_covbases*100/$genome_size;
    $f{$sample}{"depth"}=$all_covdepth/$genome_size;
    #print "$sample\t", $all_covbases/$genome_size,"\t",$all_covdepth/$genome_size,"\n";
}
my @in=<$dir/*stat>;
for my $in (@in){
    my $name = basename $in;
    my $sample=(split/\./,$name)[0];
    open IN, "$in";
    while(<IN>){
        chomp;
        #print "$sample\t$ref\t\n" if /mapped \((\S+)/;
        $f{$sample}{"rate"}=$1 if /mapped \((\S+)\%/;
    }   
    close IN;
}
for my $sample (sort keys %f){
    #print "$sample\t$f{$sample}{coverage}\t$f{$sample}{depth}\t$f{$sample}{rate}\n";
    print "$sample\t",sprintf("%.2f",$f{$sample}{coverage}),"\t",sprintf("%.2f",$f{$sample}{depth}),"\t",sprintf("%.2f",$f{$sample}{rate}),"\n";
}

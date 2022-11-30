use strict;
use SVG;
use warnings;
open LEN,$ARGV[0] or die $!;
open CNV,$ARGV[1] or die $!;
my %cnv;my %chr;
my $maxlen=0;
while (<LEN>){
    chomp;
    next if /^#/;
    my @in=split /\t/,$_;
   $maxlen=$in[1] if($maxlen<$in[1]); 
}

seek(LEN,0,0);

my $w=50;
my $h=300;
my $chrnum=0;
while (<LEN>){
    chomp;
    next if /^#/;
    my @in=split /\t/,$_;
    my $height=$in[1]/$maxlen*300;
    next if($height<5);
    $chrnum++;
    $w+=50;
    if ($w>=900){
        $w=50;
        $h+=350;
    }
}

if ($h>=600) {
	$w=850;
}
$h+=$chrnum*15;

my $svg=SVG->new(width=>$w+50,height=>$h+50);

seek(LEN,0,0);

my $x=50;
my $y=50;
my %scale;
my $id=1;
my %legend;
while (<LEN>){
    chomp;
    next if /^#/;
    my @in=split /\t/,$_;
    # $in[0] =~ s/chr// ;
    $chr{$in[0]}=$in[1];
    my $height=$in[1]/$maxlen*300;
    next if($height<5);
    $legend{$id}=$in[0];
    $svg->rect(x=>$x,y=>$y,rx=>5,ry=>5,width=>12,height=>$height,fill=>"#A9A9A9",stroke=>'black');
    $svg->text(x=>$x,y=>$y-5,"-cdata" => "$id");
    $id++;
    $scale{$in[0]}{x}=$x;
    $scale{$in[0]}{y}=$y;
    $scale{$in[0]}{h}=$height;
    $x+=50;
    if ($x>=900){
        $x=50;
        $y+=350;
    }
}
close LEN;
<CNV>;
while (<CNV>){
    chomp;
    my @in=split /\t/,$_;
	#my $cnv=$in[0];
    my $chr=$in[0];
    #$chr=~s/chr//;
    next if (!exists $chr{$chr} or $chr{$chr}/$maxlen*300<5);
    my $pos=($in[1]+$in[2])/2;
    #$cnv{$cnv}{$chr}=$pos;
    $pos=$pos/$chr{$chr}*$scale{$chr}{h}+$scale{$chr}{y};
    my $pos_1=$pos+5;
    my $pos_2=$pos-5;
    my ($x1,$x2,$x3)=($scale{$chr}{x},$scale{$chr}{x}-5,$scale{$chr}{x}-5);
    my ($x4,$x5,$x6)=($scale{$chr}{x}+12,$scale{$chr}{x}+17,$scale{$chr}{x}+17);
    if ($in[4] eq 'loss'){
        $svg->polygon(points=>"$x1,$pos $x2,$pos_1 $x3,$pos_2",fill=>"blue");
    }else{
        $svg->polygon(points=>"$x4,$pos $x5,$pos_1 $x6,$pos_2",fill=>"red");
    }
}
close CNV;

my $num=1;
my $keylen=&lenkey(%legend);
$x=$w-$keylen*11;
$y=$h;
foreach my $id(sort{ $b <=> $a } keys(%legend)){
    $svg->text(x=>$x,y=>$y,"-cdata" => "$id:$legend{$id}");
    $y-=15;
    $num++;
    if ($num>4){
        $x=$x-$keylen*11;
        $y=$h;
        $num=1;
    }
}

open OUT,">$ARGV[2]" or die $!;
my $out=$svg->xmlify;
print OUT $out;

sub lenkey{
my (%x)=@_;
my $max=0;
    foreach my $key(keys(%x)){
        my $len=length($x{$key});
        $max=$len if($len>$max);
    }
    return $max;
}

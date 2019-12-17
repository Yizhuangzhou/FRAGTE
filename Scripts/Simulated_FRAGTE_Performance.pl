#!usr/bin/perl -w
use strict;
die "perl $0 [Data/Simulated_IntraSpecies.xls][indir][output|Sensitivity][output|NumStat][output|specificity]" unless@ARGV == 5;
open(IF,$ARGV[0])||die;
my $indir=$ARGV[1];
open(OC,">$ARGV[2]")||die;
open(ON,">$ARGV[3]")||die;
open(OS,">$ARGV[4]")||die;
my $i=10;
my @name=();
while($i<= 100){
	push @name,$i;
	$i+=10;
}
my %type=();
<IF>;
my $IntraNum=0;
while(<IF>){
	chomp;
	my @a=split /\t/;
	$type{$a[0]}{$a[1]}="Intra";
	$IntraNum++;
}
close IF;
my $total=469656;
my $InterNum=$total-$IntraNum;
my %num=();
my %sensitivity=();
my %specificity=();
foreach my $p1(@name){
	foreach my $p2(@name){
		my $file="$indir/Ref_P$p1\_Query_P$p2\_Pairs_byFRAGTE.xls";
		my $TP=0;
		my $NP=0;
		my $c=0;
		open(TMP,$file)||die;
		while(<TMP>){
			chomp;
			$c++;
			my @a=split /\t/;
			if(exists $type{$a[0]}{$a[1]}){
				$TP++;
			}
			else{
				$NP++;
			}
		}
		close TMP;
		my $frac1=$TP*100/$IntraNum;
		$frac1=sprintf("%.2f",$frac1);
		$sensitivity{$p1}{$p2}=$frac1;
		my $frac2=$c*100/$total;
		$frac2=sprintf("%.2f",$frac2);
		$num{$p1}{$p2}=$frac2;
		my $frac3=($InterNum-$NP)*100/$InterNum;
		$frac3=sprintf("%.2f",$frac3);
		$specificity{$p1}{$p2}=$frac3;
	}
}
print OC join("\t",@name),"\n";
print ON join("\t",@name),"\n";
foreach my $k1(@name){
	print OC "$k1";
	foreach my $k2(@name){
		print OC "\t$sensitivity{$k1}{$k2}";
	}
	print OC "\n";
}
close OC;
foreach my $k1(@name){
	print ON "$k1";
	foreach my $k2(@name){
		print ON "\t$num{$k1}{$k2}";
	}
	print ON "\n";
}
close ON;
print OS join("\t",@name),"\n";
foreach my $k1(@name){
	print OS "$k1"; 	
	foreach my $k2(@name){
		print OS "\t$specificity{$k1}{$k2}";
	}
	print OS "\n";
}
close OS;

#!usr/bin/perl -w
use strict;
die "perl $0 [Data/MAG_IntraSpecies.xls][MAG_Pairs_byFRAGTE.xls][output]" unless @ARGV ==3;
open(IF,$ARGV[0])||die;
open(IN,$ARGV[1])||die;
open(OUT,">$ARGV[2]")||die;
my %hash;
my $n=0;
while(<IF>){
	chomp;
	my @a=split /\t/;
	$n++;
	$hash{$a[0]}{$a[1]}=1;
	$hash{$a[1]}{$a[0]}=1;
}
close IF;
my $IntraNum=$n*2;
my $total=9189992;
my $InterNum=$total-$IntraNum;
my $TP=0;
my $FP=0;
my $sum=0;
while(<IN>){
	chomp;
	my @a=split /\t/;
	next if($a[0] eq $a[1]);
	$sum++;
	if($hash{$a[0]}{$a[1]}){
		$TP++;
	}
	else{
		$FP++;
	}
}
close IN;
my $sensitivity=sprintf("%.2f",$TP*100/$IntraNum);
my $NumStat=sprintf("%.2f",$sum*100/$total);
my $specificity=sprintf("%.2f",($InterNum-$FP)*100/$InterNum);
my $TN=$InterNum-$FP;
print OUT "Sensitivity\t$TP\t$IntraNum\t$sensitivity\n";
print OUT "Specificity\t$TN\t$InterNum\t$specificity\n";
print OUT "NumStat\t$sum\t$total\t$NumStat\n";
close OUT;
close OUT;

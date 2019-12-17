#!usr/bin/perl -w
use strict;
die "perl $0 [Pairs byFRAGTE][Query GenomeInfo][Ref GenomeInfo][output]" unless @ARGV ==4;
open(IN,$ARGV[0])||die;
open(QF,$ARGV[1])||die;
open(RF,$ARGV[2])||die;
open(OUT,">$ARGV[3]")||die;
my $TP=0;
my $FP=0;
my $sum=0;
my %Qsp;
my %Rsp;
my $Qnum=0;
my $Rnum=0;
while(<QF>){
	chomp;
	my @a=split /\t/;
	$Qsp{$a[0]}=$a[1];
	$Qnum++;
}
close QF;
while(<RF>){
	chomp;
	my @a=split /\t/;
	$Rsp{$a[0]}=$a[1];
	$Rnum++;
}
close RF;
my $IntraNum=0;
foreach my $k1(keys %Qsp){
	foreach my $k2(keys %Rsp){
		if($Qsp{$k1} eq $Rsp{$k2}){
			$IntraNum++;
		}
	}
}
my $total=$Qnum*$Rnum;
my $InterNum=$total-$IntraNum;
while(<IN>){
	chomp;
	my @a=split /\t/;
	$sum++;
	if($Qsp{$a[0]} eq $Rsp{$a[1]}){
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

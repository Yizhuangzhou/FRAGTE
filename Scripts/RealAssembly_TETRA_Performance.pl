#!usr/bin/perl -w
use strict;
die "perl $0 [TETRA][Query_GenomeInfo.xls][Ref_GenomeInfo.xls][Num][output]" unless @ARGV == 5;
open(LIST,$ARGV[0])||die;
open(IF,$ARGV[1])||die;
open(IN,$ARGV[2])||die;
my $num=$ARGV[3];
open(OUT,">$ARGV[4]")||die;
my %Qsp;
my $Qnum=0;
while(<IF>){
	chomp;
	my @a=split /\t/;
	$Qsp{$a[0]}=$a[1];
	$Qnum++;
}
close IF;
my %Rsp;
my $Rnum=0;
while(<IN>){
	chomp;
	my @a=split /\t/;
	$Rsp{$a[0]}=$a[1];
	$Rnum++;
}
close IN;
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
my $n=0;# number of intraspecies pairs sieved by TETRA
my @val=();
while(<LIST>){
	chomp;
	my @a=split /\t/;
	if($Qsp{$a[0]} eq $Rsp{$a[1]}){
		push @val,$a[2];
		if($a[2] >= 0.99){
			$n++;
		}
	}
}
close LIST;

@val=sort {$b <=> $a} @val;
my $cutoff=$val[$num-1];
print STDERR "$num\t$cutoff\n";
open(LIST,$ARGV[0])||die;
my $TP=0;
my $FP=0;
my $sum=0;
open(FL,$ARGV[0])||die;
while(<FL>){
	chomp;
	my @a=split /\t/;
	if($a[2] >= $cutoff){
		$sum++;
		if($Qsp{$a[0]} eq $Rsp{$a[1]}){
			$TP++;
		}
		else{
			$FP++;
		}
	}
}
close FL;
my $sensitivity=sprintf("%.2f",$n*100/$IntraNum);
my $specificity=sprintf("%.2f",($InterNum-$FP)*100/$InterNum);
my $NumStat=sprintf("%.2f",$sum*100/$total);
my $TN=$InterNum-$FP;
print STDERR "$cutoff\n";
print OUT "Sensitivity\t$n\t$IntraNum\t$sensitivity\n";
print OUT "Specificity\t$TN\t$InterNum\t$specificity\n";
print OUT "NumStat\t$sum\t$total\t$NumStat\n";
close OUT;

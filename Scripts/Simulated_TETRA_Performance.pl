#!usr/bin/perl -w
use strict;
use File::Basename;
die "perl $0 [Data/Simulated_IntraSpecies.xls][TETRA.list][FRAGTE_performance.xls][output|sensitivity][output|numstat][output|specificity]" unless @ARGV == 6;
open(IF,$ARGV[0])||die;
open(LIST,$ARGV[1])||die;
open(IN,$ARGV[2])||die;
open(OP,">$ARGV[3]")||die;
open(ON,">$ARGV[4]")||die;
open(OS,">$ARGV[5]")||die;
my $IntraNum=0;
my $total=469656;
my %type=();
<IF>;
while(<IF>){
	chomp;
	my @a=split /\t/;
	$type{$a[0]}{$a[1]}="Intra";
	$IntraNum++;
}
close IF;
my $InterNum=$total-$IntraNum;
my %SelectIntra=();
my $f=<IN>;
chomp($f);
my @name=split /\t/,$f;
while(<IN>){
	chomp;
	my @a=split /\t/;
	my $id=shift @a;
	for(my $i=0;$i<=$#a;$i++){
		$SelectIntra{$id}{$name[$i]}=$a[$i]/100;
	}
}
close IN;
my ($len1,$len2);
my $threshold="";
while(<LIST>){
	chomp;
	my $file=$_;
	my $base=basename($file,"_TETRA.xls");
	if($base=~/Ref\_P(\d+)\_Query\_P(\d+)/){
		$len1=$1;
		$len2=$2;
		my $c=0;
		my $num=int($IntraNum*$SelectIntra{$len1}{$len2})-1;
		my @PCCD=();
		open(TMP,$file)||die;
		while(<TMP>){
			chomp;
			my @a=split /\t/; 
			if(exists $type{$a[0]}{$a[1]}){
				push @PCCD,$a[2];
			}
		}
		close TMP;
		@PCCD=sort {$b <=> $a} @PCCD;
		my $cutoff=$PCCD[$num];
		if($threshold eq ""){
			$threshold=$cutoff;
		}
		elsif($threshold >$cutoff){
			$threshold=$cutoff;
		}
	}
	else{
		print STDERR "$file\n";
	}
}
close LIST;
print STDERR "$threshold\n";
open(LIST,$ARGV[1])||die;
my %correct=();
my %num;
my %specificity=();
while(<LIST>){
	chomp;
	my $file=$_;
	my $base=basename($file,"_TETRA.xls");
	if($base=~/Ref\_P(\d+)\_Query\_P(\d+)/){
		$len1=$1;
		$len2=$2;
		my $c=0;
		my $correct=0;
		my $incorrect=0;
		my @PCCD=();
		open(TMP,$file)||die;
		while(<TMP>){
			chomp;
			my @a=split /\t/; 
			if($a[2] >= 0.99 && exists $type{$a[0]}{$a[1]}){
				$correct++;
			}
			if($a[2] >=$threshold){
				$c++;
				if(exists $type{$a[0]}{$a[1]}){
				}
				else{
					$incorrect++;
				}
			}
		}   
		close TMP;
		my $frac1=$correct*100/$IntraNum;
		$frac1=sprintf("%.2f",$frac1);
		$correct{$len1}{$len2}=$frac1;
		my $frac2=$c*100/$total;
		$frac2=sprintf("%.2f",$frac2);
		$num{$len1}{$len2}=$frac2;
		my $frac3=($InterNum-$incorrect)*100/$InterNum;
		$frac3=sprintf("%.2f",$frac3);
		 $specificity{$len1}{$len2}=$frac3;
	}
	else{
		print STDERR "$file\n";
	}
}
close LIST;
my @p=qw(10 20 30 40 50 60 70 80 90 100);
print OP join("\t",@p),"\n";
foreach my $k1(@p){
	print OP "$k1";
	foreach my $k2(@p){
		print OP "\t$correct{$k1}{$k2}";
	}
	print OP "\n";
}
close OP;
print ON join("\t",@p),"\n";
foreach my $k1(@p){
	print ON "$k1";
	foreach my $k2(@p){
		print ON "\t$num{$k1}{$k2}";
	}
	print ON "\n";
}
close ON;
print OS join("\t",@p),"\n";
foreach my $k1(@p){
	print OS "$k1";
	foreach my $k2(@p){
		print OS "\t$specificity{$k1}{$k2}";
	}
	print OS "\n";
}
close OS;

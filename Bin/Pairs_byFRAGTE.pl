#!usr/bin/perl -w
use strict;
use FindBin qw($Bin);
use File::Basename;
use POSIX qw(strftime);
use Time::Local;

die "perl $0 [Ref zvalue][Query zvalue][output]" unless @ARGV == 3;
open(IN,$ARGV[0])||die;
open(IF,$ARGV[1])||die;
open(OUT,">$ARGV[2]")||die;
my $start_time = strftime "%Y-%m-%d %H:%M:%S", localtime;
my $t1=timelocal(localtime());
print STDERR "$start_time: Start $0\n";

######################################################## Reading Cutoff  ##############################################
open(TMP,"$Bin/Cutoff.xls")||die;
my %LSC;
<TMP>;
while(<TMP>){
	chomp;
	my @a=split /\t/;
	my @b=split /\s*vs\s*/,$a[0];
	my @cutoff=($a[5],$a[6]);
	@cutoff=sort {$a <=> $b} @cutoff;
	my $cutoff=$cutoff[0];
	if($cutoff=~/(0\.\d{2})/){
		$LSC{$b[0]}{$b[1]}=$1;
		$LSC{$b[1]}{$b[0]}=$1;
	}
}
close TMP;

######################################################## MAIN  #######################################################
my %spid=();
my %len=();
my %GSC=();
my @ID=();
my %refyy=();
my %SSyy=();
my %Prefyy=();
my %PSSyy=();
my %plen=();
while(<IN>){
	chomp;
	my @array=split /\#/;
	my @a=split /\t/,$array[0];
	my ($id,$spid,$plen,$len,$GSC)=@a;
	push @ID,$id;
	$spid{$id}=$spid;
	$len{$id}=$len;
	$GSC{$id}=$GSC;
	my @b=split /\t/,$array[2];
	my $refxx=\@b;
	my $SSxx=&SS($refxx,$refxx,__LINE__);
	$refyy{$id}=$refxx;
	$SSyy{$id}=$SSxx;
	if($len <= 200){
		$plen{$id}=$len;
		$Prefyy{$id}=$refxx;
		$PSSyy{$id}=$SSxx;
	}
	else{
		$plen{$id}=$plen;
		my @c=split /\t/,$array[1];
		my $refxx2=\@c;
		my $SSxx2=&SS($refxx2,$refxx2,__LINE__);
		$Prefyy{$id}=$refxx2;
		$PSSyy{$id}=$SSxx2;
	}
}
close IN;

my $threshold="";
my $len;
my ($refxx,$SSxx);
while(<IF>){
	chomp;
	my @array=split /\#/;
	my @a=split /\t/,$array[0];
	my ($id,$spid,$Qplen,$Qlen,$GSC)=@a;
	my @b=split /\t/,$array[2];
	my $refxx2=\@b;
	my $SSxx2=&SS($refxx2,$refxx2,__LINE__);
	if($Qlen <= 200){
		$Qplen=$Qlen;
		$refxx=$refxx2;
		$SSxx=$SSxx2;
	}
	else{
		my @c=split /\t/,$array[1];
		$refxx=\@c;
		$SSxx=&SS($refxx,$refxx,__LINE__);
	}
	
	foreach my $k (@ID){
		next if($id eq $k);
		my $Rlen=$len{$k};
		my $Rplen=$plen{$k};
		my $refyy=$Prefyy{$k};
		my $SSyy=$PSSyy{$k};
		my $SSxy=&SS($refxx,$refyy,__LINE__);
		my $PCCD=&correlation($SSxx,$SSyy,$SSxy);
		$PCCD=sprintf("%.2f",$PCCD);
		if($PCCD >=$LSC{$Qplen}{$Rplen} ){
			my $PCCD2=$PCCD;
			if($Qlen != $Qplen || $Rlen != $Rplen){
				my $refyy2=$refyy{$k};
				my $SSyy2=$SSyy{$k};
				my $SSxy2=&SS($refxx2,$refyy2,__LINE__);
				my $TmpPCCD=&correlation($SSxx2,$SSyy2,$SSxy2);
				if($TmpPCCD >$PCCD2){
					$PCCD2=$TmpPCCD;
				}
				$PCCD2=sprintf("%.2f",$PCCD2);
			}
			if($PCCD2 >= $GSC || $PCCD2 >= $GSC{$k}){
				print OUT "$id\t$k\t$PCCD\t$PCCD2\t$Qplen\t$Rplen\t$LSC{$Qplen}{$Rplen}\t$Qlen\t$Rlen\t$GSC\t$GSC{$k}\n";
			}
		}
	}
	my $current_time = strftime "%Y-%m-%d %H:%M:%S", localtime;
	print STDERR "$current_time: $id\n";
}
close IF;
close OUT;
my $finish_time = strftime "%Y-%m-%d %H:%M:%S", localtime;
my $t2=timelocal(localtime());
print STDERR "$finish_time: Finish\n";
my $s=$t2-$t1;
my $hour=int($s / 3600);
my $ms=$s % 3600;
my $minutes=int($ms / 60);
my $second=$ms % 60;
print STDERR "Time elaspe in seconds: $s\n";
print STDERR "Time elaspe: $hour\:$minutes\:$second\n";

###################################################### Subroutine for computing mean##########################################################
sub mean{
	my ($ref)=@_;
	my @a=@{$ref};
	my $sum=0;
	foreach(@a){
		$sum+=$_;
	}
	my $mean=$sum/@a;
	return $mean;
}

###################################################### Subroutine for computing covariance or variance##########################################
sub SS{
	my ($ref1,$ref2,$line)=@_;
	if(!$ref1){
		print STDERR "$line\n";
	}
	my @a1=@{$ref1};
	if(!$ref2){
		print STDERR "$line\n";
	}
	my @a2=@{$ref2};
	my $sum=0;
	for (my $k=0;$k<=$#a1;$k++){
		$sum+=$a1[$k]*$a2[$k];
	}
	return $sum;
}

##################################################### Subroutine for computing Pearson correlation coefficient ##################################
sub correlation {
	my ($ssxx,$ssyy,$ssxy)=@_;
	my $correlation=0;
	if($ssxy !=0 && $ssxx !=0 && $ssyy!=0){
		my $sign=$ssxy/abs($ssxy);
		$correlation=$sign*sqrt($ssxy*$ssxy/($ssxx*$ssyy));
	}
	else{
		$correlation="0";
	}
	return $correlation;
}



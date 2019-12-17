#!/usr/bin/perl
use strict;
use warnings;
use POSIX qw(strftime);
use Time::Local;

die "perl $0 [Ref_Zvalue.xls][Query_Zvalue.xls][output]" unless (@ARGV==3);
my $start_time = strftime "%Y-%m-%d %H:%M:%S", localtime;
my $t1=timelocal(localtime());
print STDERR "$start_time: Start\n";

######################################################## Intializing for Zvalue ##############################################
my $km=4;
my @mono=("A","T","G","C");
our @oligo_kmer=();
our (@oligo_k_1mer,@oligo_k_2mer)=((),());
my @nuc=@mono;
my $index=1;
while(){
	my @tmpnuc=();
	$index++;
	foreach(@mono){
		my $word=$_;
		foreach (@nuc){
			my $nuc=$word."$_";
			push @tmpnuc,$nuc;
		}
	}
	@nuc=@tmpnuc;
	if($index == $km-2){
		@oligo_k_2mer=@nuc;
	}
	elsif($index == $km -1){
		@oligo_k_1mer=@nuc;
	}
	elsif($index == $km){
		@oligo_kmer=@nuc;
		last;
	}
}
my $current_time = strftime "%Y-%m-%d %H:%M:%S", localtime;
print STDERR "$current_time: Finish Intializing for zvalues\n";

##################################################################### Sequences  ##################################################
open(LIST,$ARGV[0])||die;
my @ID=();
my (%refxx,%SSxx);
while(<LIST>){
	chomp;
	my @a=split /\t/;
	my $id=shift @a;
	push @ID,$id;
	my $spid=shift @a;
	my $zvalueref=\@a;
	my $mean=&mean($zvalueref);
	my @d=();
	foreach(@a){
		my $tmp=$_-$mean;
		push @d,$tmp;
	}
	my $refxx=\@d;
	my $SSxx=&SS($refxx,$refxx,__LINE__);
	$refxx{$id}=$refxx;
	$SSxx{$id}=$SSxx;
}
close LIST;

#################################################### correlation for same genus  #################################################
open(IF,$ARGV[1])||die;
open(OUT,">$ARGV[2]")||die;
my $flag=0;
while(<IF>){
	chomp;
	my @a=split /\t/;
	my $id=shift @a;
	my $spid=shift @a;
	my $zvalueref=\@a;
	my $mean=&mean($zvalueref);
	my @d=();
	foreach(@a){
		my $tmp=$_-$mean;
		push @d,$tmp;
	}
	my $refxx=\@d;
	my $SSxx=&SS($refxx,$refxx,__LINE__);
	foreach my $k (@ID){
		next if($id eq $k);
		my $refyy=$refxx{$k};
		my $SSyy=$SSxx{$k};
		my $SSxy=&SS($refxx,$refyy,__LINE__);
		my $cor=&correlation($SSxx,$SSyy,$SSxy);
		$cor=sprintf("%.2f",$cor);
		print OUT "$id\t$k\t$cor\n";
	}
	$current_time = strftime "%Y-%m-%d %H:%M:%S", localtime;
	print STDERR "$current_time: $id\n";
}
close OUT;
close IF;
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
	my $correl=0;
	if($ssxy !=0 && $ssxx !=0 && $ssyy!=0){
		my $sign=$ssxy/abs($ssxy);
		$correl=$sign*sqrt($ssxy*$ssxy/($ssxx*$ssyy));
	}
	else{
		$correl="0";
	}
	return $correl;
}



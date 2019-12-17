#!/usr/bin/perl
use strict;
use warnings;
use POSIX qw(strftime);
use Time::Local;

die "perl $0 [*_GenomeSize.xls][output|zvalue]" unless (@ARGV==2);
my $start_time = strftime "%Y-%m-%d %H:%M:%S", localtime;
my $t1=timelocal(localtime());
print STDERR "$start_time: Start ...\n";

######################################################## Intializing for Zvalue ##############################################
my $km=4;
my @mono=("A","T","G","C");
my @oligo_kmer=();
my (@oligo_k_1mer,@oligo_k_2mer)=((),());
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
print STDERR "$current_time: Finish intialization\n";

##################################################################### Sequences  ##################################################
open(IF,$ARGV[0])||die;
open(OUT,">$ARGV[1]")||die;
my $n=0;
while(<IF>){
	chomp;
	my @a=split /\t/;
	open(FA,$a[6])||die;
	my $seq="";
	while(<FA>){
		chomp;
		if(/>(\S+)/){
			next;
		}
		else{
			s/[^ATGC]//gi;
			$seq.=uc($_);
		}
	}
	close FA;
	my $ref=&NormalizedZvalue($seq);
	my @zvalue=@{$ref};
	print OUT "$a[0]\t$a[1]\t",join("\t",@zvalue),"\n";
	$current_time = strftime "%Y-%m-%d %H:%M:%S", localtime;
	print STDERR "$current_time: $a[0]\n";
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

############################################################## Computing zvalue ##########################################################
sub NormalizedZvalue{
	my ($seq)=@_;
	my $tmpseq=reverse($seq);
	$tmpseq=~tr/ATGC/TACG/;
	$seq.=$tmpseq;
	my (%kmer,%k_1mer,%k_2mer)=((),(),());
	$seq=~s/[^ATGC]//gi;
	my $len=length $seq;
	my $weight=($len)**0.5;
	foreach(@oligo_kmer){
		$kmer{$_}=1;
	}
	foreach(@oligo_k_1mer){
		$k_1mer{$_}=1;
	}
	foreach(@oligo_k_2mer){
		$k_2mer{$_}=1;
	}
	for(my $i=0;$i<=$len-$km;$i++){
		my $sub=substr($seq,$i,$km);
		$kmer{$sub}++;
	}
	for(my $i=0;$i<=$len-($km-1);$i++){
		my $sub=substr($seq,$i,$km-1);
		$k_1mer{$sub}++;
	}
	for (my $i=0;$i<=$len-($km-2);$i++){
		my $sub=substr($seq,$i,$km-2);
		$k_2mer{$sub}++;
	}
	my @zvalue=();
	foreach(@oligo_kmer){
		my $N_koligo=$kmer{$_};
		my $N_former=$k_1mer{substr($_,0,$km-1)};
		my $N_latter=$k_1mer{substr($_,1,$km-1)};
		my $N_midder=$k_2mer{substr($_,1,$km-2)};
		my $denominator=($N_midder-$N_former)*($N_midder-$N_latter)*$N_former*$N_latter;
		my $sqrt=$denominator**0.5;
		if($sqrt !=0 and $N_midder !=0){
			my $zvalue=$N_midder**0.5*($N_koligo*$N_midder-$N_former*$N_latter)/($sqrt*$weight);
			push @zvalue,$zvalue;
		}
		else{
			push @zvalue,0;
		}
	}
	return \@zvalue;
}


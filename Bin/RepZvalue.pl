#!/usr/bin/perl
use strict;
use warnings;
use FindBin qw($Bin);
use POSIX qw(strftime);
use Time::Local;

die "perl $0 [GenomeInfo table][outdir][prefix]" unless (@ARGV==3);
my $start_time = strftime "%Y-%m-%d %H:%M:%S", localtime;
my $t1=timelocal(localtime());
print STDERR "$start_time: Start $0 $ARGV[0] $ARGV[1] $ARGV[2] ...\n";

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
my @UsedKmer=();
my %tag=();
foreach my $k1 (@oligo_kmer) {
	my $k2=reverse($k1);
	$k2=~tr/ATGC/TACG/;
	if($tag{$k2}){
		next;
	}
	else{
		push @UsedKmer,$k1;
		$tag{$k1}=1;
	}
}
my @UsedKmer1=();
%tag=();
foreach my $k1 (@oligo_k_1mer) {
	my $k2=reverse($k1);
	$k2=~tr/ATGC/TACG/;
	if($tag{$k2}){
		next;
	}
	else{
		push @UsedKmer1,$k1;
		$tag{$k1}=1;
	}
}
my @UsedKmer2=();
%tag=();
foreach my $k1 (@oligo_k_2mer) {
	my $k2=reverse($k1);
	$k2=~tr/ATGC/TACG/;
	if($tag{$k2}){
		next;
	}
	else{
		push @UsedKmer2,$k1;
		$tag{$k1}=1;
	}
}
my $current_time = strftime "%Y-%m-%d %H:%M:%S", localtime;
print STDERR "$current_time: Finish intialization\n";

##################################################################### Reading Cutoff #############################################
open(TMP,"$Bin/Cutoff.xls")||die;
<TMP>;
my %cutoff;
while(<TMP>){
	chomp;
	my @a=split /\t/;
	my @b=split /\s*vs\s*/,$a[0];
	my @cutoff=($a[5],$a[6]);
	@cutoff=sort{$a <=> $b} @cutoff;
	if($b[0] eq $b[1]){
		if($cutoff[0]=~/(0\.\d{2})/){
			$cutoff[0]=$1;
		}
		$cutoff{$b[0]}=$cutoff[0];
	}
}
close TMP;

##################################################################### Sequences  ##################################################
open(IF,$ARGV[0])||die;
my $dir=$ARGV[1];
mkdir $dir unless (-e $dir);
my $prefix=$ARGV[2];
open(OUT,">$dir/$prefix\_RepZvalue.xls")||die;
while(<IF>){
	chomp;
	my @a=split /\t/;
	my $id=$a[0];
	my $spid=$a[1];
	next if($a[5]<10000);
	my $file=$a[6];
	open(FA,$file)||die;
	my $seq="";
	while(<FA>){
		chomp;
		if(/^>(\S+)/){
			next;
		}
		else{
			s/[^ATGC]//gi;
			$seq.=uc($_);
		}
	}
	close FA;
	my $len=length $seq;
	next if($len<10000);

	if($len <40000){
		my $kb=int($len/10000)*10;
		my $threshold=$cutoff{$kb};
		my $ref=&Fragment1($seq,$len);
		print OUT "$id\t$spid\t$kb\t$kb\t$threshold";
		my @refxx=@{$ref};
		print OUT "#",join("\t",@refxx);
		print OUT "#",join("\t",@refxx),"\n";
	}
	else{
		my $sublen=int($len/4);
		if($sublen >200000){
			$sublen=200000;
		}
		my $kb1=int($sublen/10000)*10;
		my $kb2=int($sublen*4/10000)*10;
		my ($threshold,$ref1,$ref2)=&Fragment2($seq,$len,$sublen);
		print OUT "$id\t$spid\t$kb1\t$kb2\t$threshold";
		my @refxx1=@{$ref1};
		print OUT "#",join("\t",@refxx1);
		my @refxx2=@{$ref2};
		print OUT "#",join("\t",@refxx2),"\n";
	}
	$current_time = strftime "%Y-%m-%d %H:%M:%S", localtime;
	print STDERR "$current_time: $id\n";
}
close IF;
close OUT;
my $finish_time = strftime "%Y-%m-%d %H:%M:%S", localtime;
my $t2=timelocal(localtime());
print STDERR "$finish_time: Finish $0 $ARGV[0] $ARGV[1] $ARGV[2]\n";
my $s=$t2-$t1;
my $hour=int($s / 3600);
my $ms=$s % 3600;
my $minutes=int($ms / 60);
my $second=$ms % 60;
print STDERR "Time elaspe in seconds: $s\n";
print STDERR "Time elaspe: $hour\:$minutes\:$second\n";

############################################################## Computing zvalue ##########################################################
sub zvalue{
	my ($ref1,$ref2,$ref3)=@_;
	my (%kmer,%k_1mer,%k_2mer)=((),(),());
	%kmer=%{$ref1};
	%k_1mer=%{$ref2};
	%k_2mer=%{$ref3};
	my @Zvalue=();
	foreach my $k1 (@UsedKmer){
		my $k2=reverse($k1);
		$k2=~tr/ATGC/TACG/;
		my $N_koligo=$kmer{$k1};
		my $N_former=$k_1mer{substr($k1,0,$km-1)};
		my $N_latter=$k_1mer{substr($k1,1,$km-1)};
		my $N_midder=$k_2mer{substr($k1,1,$km-2)};
		my $denominator=($N_midder-$N_former)*($N_midder-$N_latter)*$N_former*$N_latter;
		my $sqrt=$denominator**0.5;
		if($sqrt !=0 and $N_midder !=0){
			my $zvalue=$N_midder**0.5*($N_koligo*$N_midder-$N_former*$N_latter)/$sqrt;
			push @Zvalue,$zvalue;
			if($k1 ne $k2){
				push @Zvalue,$zvalue;
			}
		}
		else{
			push @Zvalue,0;
			if($k1 ne $k2){
				push @Zvalue,0;
			}
		}
	}
	return \@Zvalue;
}

###################################################### Subroutine for computing mean##############################################
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

###################################################### Computing covariance or variance ##########################################
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

##################################################### Computing Pearson correlation coefficient ###################################
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

####################################################################### SD #####################################################
sub sd{
	my ($ref,$mean)=@_;
	my @val=@{$ref};
	my $sum=0;
	foreach(@val){
		$sum +=($_-$mean)**2;
	}
	$sum /=$#val;
	$sum =$sum **0.5;
	return $sum;
}

####################################################################### Fragment #####################################################
sub Fragment1{
	my ($seq,$len)=@_;
	my %TNF=();
	for(my $j=0;$j<=$len-$km;$j++){
		my $sub=substr($seq,$j,$km);
		$TNF{$sub}++;
	}

	my (%SUM,%SUM1,%SUM2)=((),(),());
	foreach my $k1 (@UsedKmer){
		my $k2=reverse($k1);
		$k2=~tr/ATGC/TACG/;
		if($TNF{$k1}){
			$SUM{$k1}+=$TNF{$k1};
		}
		if($TNF{$k2}){
			$SUM{$k1}+=$TNF{$k2};
		}
		if(!$SUM{$k1}){
			$SUM{$k1}=1;
		}
		$SUM{$k2}=$SUM{$k1};
	}
	foreach my $kmer(@oligo_kmer) {
		for(my $j=0;$j<=$km-($km-1);$j++){
			my $sub=substr($kmer,$j,$km-1);
			if($SUM{$kmer} >1){
				$SUM1{$sub}+=$SUM{$kmer};
			}
		}
		for(my $j=0;$j<=$km-($km-2);$j++){
			my $sub=substr($kmer,$j,$km-2);
			if($SUM{$kmer}>1){
				$SUM2{$sub}+=$SUM{$kmer};
			}
		}
	}
	foreach my $k (@oligo_k_1mer){
		if(!$SUM1{$k}){
			$SUM1{$k}=1;
		}
		else{
			$SUM1{$k}=int(($SUM1{$k}+1)/2);
		}
	}
	foreach my $k (@oligo_k_2mer){
		if(!$SUM2{$k}){
			$SUM2{$k}=1;
		}
		else{
			$SUM2{$k}=int(($SUM2{$k}+2)/3);
		}
	}

	my $zvalueref=&zvalue(\%SUM,\%SUM1,\%SUM2);
	my $mean=&mean($zvalueref);
	my @f=();
	foreach(@{$zvalueref}){
		my $tmp=$_-$mean;
		push @f,$tmp;
	}
	my $ref=\@f;
	return ($ref);
}

####################################################################### Fragment #####################################################
sub Fragment2{
	my ($seq,$len,$sublen)=@_;
	my $addlen=int($sublen/2);
	my $kb=int($sublen/10000)*10;
	#my $addlen=$kb*1000/2;
	###TNF
	my %TNF=();
	my $index=0;
	my @index=();

	for(my $i=0;$i<=$len-$addlen;$i+=$addlen){
		$index++;
		push @index,$index;
		my $fragment=substr($seq,$i,$addlen);
		for(my $j=0;$j<=$addlen-$km;$j++){
			my $sub=substr($fragment,$j,$km);
			$TNF{$index}{$sub}++;
		}
	}

	## To obtain the representative fragment
	my $num=0;
	my @num=();
	my (%refxx,%SSxx)=((),());
	my $index2="";
	my %Index=();
	for (my $ind=0;$ind<=$#index;$ind++){
		$num++;
		push @num,$num;
		my %kmer=();
		my %k_1mer=();
		my %k_2mer=();
		my $index=$index[$ind];
		$Index{$num}=$index;
		### To obtain km frequency
		foreach my $k1 (@UsedKmer){
			my $k2=reverse($k1);
			$k2=~tr/ATGC/TACG/;
			foreach (0 .. 1) {
				$index2=$index+$_;
				if($index2 == @index +1){#To link the start and end sequences
					$index2 =1;
				}
				if($TNF{$index2}{$k1}){
					$kmer{$k1}+=$TNF{$index2}{$k1};
				}
				if($TNF{$index2}{$k2}){
					$kmer{$k1}+=$TNF{$index2}{$k2};
				}
			}
			if(! $kmer{$k1}){
				$kmer{$k1}=1;
			}
			$kmer{$k2}=$kmer{$k1};
		}
		
		### To obtain km-1 and km-2 frequencies
		foreach my $kmer(@oligo_kmer) {
			for(my $j=0;$j<=$km-($km-1);$j++){
				my $sub=substr($kmer,$j,$km-1);
				if($kmer{$kmer} >1){
					$k_1mer{$sub}+=$kmer{$kmer};
				}
			}
			for(my $j=0;$j<=$km-($km-2);$j++){
				my $sub=substr($kmer,$j,$km-2);
				if($kmer{$kmer}>1){
					$k_2mer{$sub}+=$kmer{$kmer};
				}
			}
		}
		foreach my $k (@oligo_k_1mer){
			if(!$k_1mer{$k}){
				$k_1mer{$k}=1;
			}
			else{
				$k_1mer{$k}=int(($k_1mer{$k}+1)/2);
			}
		}
		foreach my $k (@oligo_k_2mer){
			if(!$k_2mer{$k}){
				$k_2mer{$k}=1;
			}
			else{
				$k_2mer{$k}=int(($k_2mer{$k}+2)/3);
			}
		}

		my $zvalueref=&zvalue(\%kmer,\%k_1mer,\%k_2mer);
		my $mean=&mean($zvalueref);
		my @f=();
		foreach(@{$zvalueref}){
			my $tmp=$_-$mean;
			push @f,$tmp;
		}
		my $refxx=\@f;
		$refxx{$num}=$refxx;
		my $SSxx=&SS($refxx,$refxx,__LINE__);
		$SSxx{$num}=$SSxx;
	}

	#Determine the genome-specific cutoff
	my %TotalPCCD=();
	my %PCCD=();
	for(my $i=0;$i <=$#num-2;$i++){
		my $n1=$num[$i];
		my $refxx=$refxx{$n1};
		my $SSxx=$SSxx{$n1};
		for(my $j=$i+2;$j <=$#num;$j++){
			my $n2=$num[$j];
			next if($n1 == 1 && $n2 == @num);
			my $refyy=$refxx{$n2};
			my $SSyy=$SSxx{$n2};
			my $SSxy=&SS($refxx,$refyy,__LINE__);
			my $PCCD=&correlation($SSxx,$SSyy,$SSxy);
			push @{$PCCD{$n1}},$PCCD;
			push @{$PCCD{$n2}},$PCCD;
			$TotalPCCD{$n1}+=$PCCD;
			$TotalPCCD{$n2}+=$PCCD;
		}
	}

	@num=sort {$TotalPCCD{$b} <=> $TotalPCCD{$a}} @num;
	my $RefIndex= $num[0];
	my @PCCD=@{$PCCD{$RefIndex}};
	my $MeanPCCD=&mean(\@PCCD);
	my $SD=&sd(\@PCCD,$MeanPCCD);
	my $threshold=$MeanPCCD-2*$SD;
	if($threshold=~/(0\.\d{2})/){
		$threshold=$1;
	}
	if($threshold <$cutoff{$kb}){
		$threshold=$cutoff{$kb};
	}
	elsif($threshold >0.92){
		$threshold=0.92;
	}
	my $ref1=$refxx{$RefIndex};
	
	my @SelectedIndex=();
	my %tag=();
	foreach my $num (@num) {
		my $index=$Index{$num};
		if(!$tag{$index}){
			push @SelectedIndex,$index;
			$tag{$index}=1;
			if(scalar @SelectedIndex == 8){
				last;
			}
		}
		my $index2=$index+1;
		if($index2 == @index +1){
			$index2=1;
		}
		if(!$tag{$index2}){
			push @SelectedIndex,$index2;
			$tag{$index2}=1;
			if(scalar @SelectedIndex == 8){
				last;
			}
		}
	}
	my (%SUM,%SUM1,%SUM2)=((),(),());
	foreach my $k1 (@UsedKmer){
		my $k2=reverse($k1);
		$k2=~tr/ATGC/TACG/;
		foreach my $index (@SelectedIndex) {
			if($TNF{$index}{$k1}){
				$SUM{$k1}+=$TNF{$index}{$k1};
			}
			if($TNF{$index}{$k2}){
				$SUM{$k1}+=$TNF{$index}{$k2};
			}
		}
		if(! $SUM{$k1}){
			$SUM{$k1}=1;
		}
		$SUM{$k2}=$SUM{$k1};
	}
	foreach my $kmer(@oligo_kmer) {
		for(my $j=0;$j<=$km-($km-1);$j++){
			my $sub=substr($kmer,$j,$km-1);
			if($SUM{$kmer} >1){
				$SUM1{$sub}+=$SUM{$kmer};
			}
		}
		for(my $j=0;$j<=$km-($km-2);$j++){
			my $sub=substr($kmer,$j,$km-2);
			if($SUM{$kmer}>1){
				$SUM2{$sub}+=$SUM{$kmer};
			}
		}
	}
	foreach my $k (@oligo_k_1mer){
		if(!$SUM1{$k}){
			$SUM1{$k}=1;
		}
		else{
			$SUM1{$k}=int(($SUM1{$k}+1)/2);
		}
	}
	foreach my $k (@oligo_k_2mer){
		if(!$SUM2{$k}){
			$SUM2{$k}=1;
		}
		else{
			$SUM2{$k}=int(($SUM2{$k}+2)/3);
		}
	}
	my $zvalueref=&zvalue(\%SUM,\%SUM1,\%SUM2);
	my $mean=&mean($zvalueref);
	my @f=();
	foreach(@{$zvalueref}){
		my $tmp=$_-$mean;
		push @f,$tmp;
	}
	my $ref2=\@f;
	return($threshold,$ref1,$ref2);
}




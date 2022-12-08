package My;

use 5.026002;
use strict;
use warnings;
require Exporter;

our @ISA = qw(Exporter);

# Items to export into callers namespace by default. Note: do not export
# names by default without a very good reason. Use EXPORT_OK instead.
# Do not simply export all your public functions/methods/constants.

# This allows declaration	use My ':all';
# If you do not need this, moving things directly into @EXPORT or @EXPORT_OK
# will save memory.
our %EXPORT_TAGS = ( 'all' => [ qw(
	
) ] );

our @EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );

our @EXPORT = qw(sum rev aln clt tran cf np ny erf mi nr nu nx ab panseq coreseq
	
);

our $VERSION = '0.01';


# Preloaded methods go here.
package My;
sub max
{
my @out=sort {$a<=>$b}@_;
return $out[-1];
}
sub uniq
{
my @t;
my %c;
@t=grep{++$c{$_}==1;}@_;
return @t;
}


sub aln
{
my @query=split(//,$_[0]);
unshift(@query,"q0");
my @subject=split(//,$_[1]);
unshift(@subject,"s0");
my %h;
my %rq;
my %rs;
my %check;
my $score;
my $n=0;
my $max=0;
my $qs="";
my $ss="";
my ($mq,$ms)=(0,0);

###打分矩阵
for(0..$#query)
	{
	my $t1=$_;
	$h{$t1}{0}=0;
    }
for(0..$#subject)
	{
	my $t2=$_;
     $h{0}{$t2}=0;
     }

#罚分
my $w1=2;


for(1..$#query)
{
  my $qn=$_;
   for(1..$#subject)
   {
    my $sn=$_;
    my $s;


	if($query[$qn] eq $subject[$sn])
      {
		###匹配得分
      $s=3;  
      } 
      else
      {
		###不匹配罚分
      $s=-3;
      }

#动态规划算法
	my @t=(($h{$qn-1}{$sn-1}+$s), ($h{$qn}{$sn-1}-$w1),($h{$qn-1}{$sn}-$w1),0);
	$h{$qn}{$sn}=&max(@t);
#	print "\t$h{$qn}{$sn}";

####追溯记录
		if($h{$qn}{$sn}==$t[3])
	   {
	$rq{$qn}{$sn}=0;
	$rs{$qn}{$sn}=0;
	$check{$qn}{$sn}=3;
       }
		if($h{$qn}{$sn}==$t[2])
	   {
	$rq{$qn}{$sn}=$qn-1;
	$rs{$qn}{$sn}=$sn;
	$check{$qn}{$sn}=2;
       }
	    if($h{$qn}{$sn}==$t[1])
	   {
	$rq{$qn}{$sn}=$qn;
	$rs{$qn}{$sn}=$sn-1;
	$check{$qn}{$sn}=1;
       }
	   	if($h{$qn}{$sn}==$t[0])
	   {
	$rq{$qn}{$sn}=$qn-1;
	$rs{$qn}{$sn}=$sn-1;
	$check{$qn}{$sn}=0;
       }
#print "#$rq{$qn}{$sn}#$rs{$qn}{$sn}#$check{$qn}{$sn}";
##记录max
   if($h{$qn}{$sn}>$max)
	   {
	   $max=$h{$qn}{$sn};
	   $mq=$qn;
	   $ms=$sn;
       }
   }

 }


######追溯

while(exists $check{$mq}{$ms})
	{
    
	$n++ if($query[$mq] eq $subject[$ms]);
	if($check{$mq}{$ms}==0)
		{
	$qs="$query[$mq]".$qs;
	$ss="$subject[$ms]".$ss;
	    }
	if($check{$mq}{$ms}==1)
		{
	$qs="-".$qs;
	$ss="$subject[$ms]".$ss;
	    }
	if($check{$mq}{$ms}==2)
		{
	$qs="$query[$mq]".$qs;
	$ss="-".$ss;
	    }


######回溯终止
	if(exists $rq{$mq}{$ms})
		{
		my $temp1=$mq;
		my $temp2=$ms;
		$mq=$rq{$temp1}{$temp2};
		$ms=$rs{$temp1}{$temp2};
		}
	else
		{
		$mq="NA";
		$ms="NA";
		}
    }

my $qid="NA";
$qid=$n/$#query if($#query>0);
my $sid="NA";
$sid=$n/$#subject if($#subject>0);
my $output ="Query Identity is :\n$qid\nSuject Identity is :\n $sid\nalignment:\n$qs\n$ss\n"; 
my @out=($qid,$sid,$qs,$ss);
return @out;
}

##正态函数值，输入至少三个数字（u值 西格玛值 x值1 ...x值n），输出（y值1..y值n）
sub ny
{
@_=grep{/^-?\d+(\.\d+)?(e[-|+]\d+)?$/}@_;
if($_[1]<=0)
	{
return "第二个xigema参数必须为正数！\n";
    }
if(@_<3)
	{
	return "Error!参数数目至少需要3个数字！\n";
    }
my $pi=3.14159265358979;
my $u=shift @_;
my $xgm=shift @_;
my @x=@_;
my @out;
for(@x)
	{
	my $f=($_-$u)**2/(2*($xgm**2));
    $f=exp(-$f);
	$f=$f/(2*$pi)**0.5;
	$f=$f/$xgm;
	push (@out,$f);
	}
return @out;
}

###统计正态分布左边概率（正态面积），输入至少3个数字（u值 西格玛值 x值1 ...x值n），输出（P值1..P值n）
sub np
{
@_=grep{/^-?\d+(\.\d+)?(e[-|+]\d+)?$/}@_;
my @sum;
if($_[1]<=0)
	{
return "第二个xigema参数必须为正数！\n";
    }
if(@_<3)
	{
	return "Error!参数数目至少需要3个数字！\n";
    }

my $u=shift @_;
my $xgm=shift @_;
my $pi=3.14159265358979;
my $t||=4;
for(@_)
	{
	my $n=$_-$u;
	my $m="r";
	if($n<0)
		{
		$n=abs($n);
	    $m="l";
		}
	$n=99 if($n>99);
##四舍五入
	$n=$n*10**$t;
	$n=int($n+0.5);
    my @in=map{$_=$_*10**-$t}(1..$n);
    if(@in==0)
		{
        push(@sum,0.5);
		next;
	    }
	@in=((0,$xgm),@in);
    my @out=&ny(@in);
    my $sum;
    map{$sum+=$_}@out;
	$sum=$sum*10**-$t;
	if($m eq "l")
		{
	    $sum=0.5-$sum;
		}
		elsif($m eq "r")
		{
		$sum=0.5+$sum; 
		}
	push(@sum,$sum);
	}
return @sum;
}


###误差函数 输入至少1个数字（等分t x值）
sub erf
{
my $pi=3.14159265358979;
my $x=shift @_;
my $sum;
for(0..10)
	{
	my $mi=&mi($_);
    $sum+=($x**(2*$_+1)*(-1)**$_)/$mi*(2*$_+1);
	}
$sum=2*$sum/$pi**0.5;
return $sum;
}

###阶乘函数 计算输入x的阶乘
sub mi
{
  my $n = int($_[0]);
  return 1 if($n == 0);
  return &mi($n-1) * $n;
}

###正态函数原函数
sub nr
{
@_=grep{/^-?\d+(\.\d+)?(e[-|+]\d+)?$/}@_;
if($_[1]<=0)
	{
return "第二个xigema参数必须为正数！\n";
    }
if(@_<3)
	{
	return "Error!参数数目至少需要3个数字！\n";
    }

my $u=shift @_;
my $xgm=shift @_;
my $x=shift @_;
my $pi=3.14159265358979;
my $nr;
my $erf=(($x-$u)/$xgm*2**0.5);
$erf=&erf(($erf));
$nr=0.5*$erf;
return $nr;
}

###根据概率推导正态函数输入至少3个数字（u值 西格玛值 P值），输出（X值）)
sub nx
{
@_=grep{/^-?\d+(\.\d+)?(e[-|+]\d+)?$/}@_;
return "Na" if(@_<3);
my $u=shift @_;
my $xgm=shift @_;
my $p=shift @_;
my $n=0;
my($out,$up,$down);
my $x;
if($p<0.5)
	{
    $down=$u-5*$xgm;
	$up=$u;
	my $dp=(&np(($u,$xgm,$down)))[0];
	if($p<$dp)
		{
	    return $down; 
	    }
	}
elsif($p>0.5)
	{
    $up=$u+5*$xgm;
	$down=$u;
    my $upp=(&np(($u,$xgm,$up)))[0];
	if($p>$upp)
		{
	    return $up; 
	    }
    }
else
	{
	return $u;
    }

while($n<20)
	{
    $x=$down+rand(1)*($up-$down);
	my @out=&np(($u,$xgm,$x));
	if(abs($out[0]-$p)<0.0001)
		{
		return $x;
	    }

	if($out[0]>$p)
		{
		$up=$x;
		$n++;
		}
	elsif($out[0]<$p)
		{
		$down=$x;
		$n++;
		}
    else
		{
		$out=$x;
		return $out;
	    }
	}
  return $out;  
}


####输入两组数据x y（，分割），使用最小二乘法计算线性拟合的a b 值
sub ab
{
my $x=shift @_;
my $y=shift @_;
my @x=split(/\,/,$x);
my @y=split(/\,/,$y);
@x=grep{/^-?\d+(\.\d+)?$/}@x;
@y=grep{/^-?\d+(\.\d+)?$/}@y;
return "x y数列至少需要2组，且数目相同\n" if(@x != @y || @x<2);
my $n=@x;
my $sx=&sum(@x);
my $sy=&sum(@y);
my $sxy;
map{$sxy+=$x[$_]*$y[$_];}(0..$#x);
my $sxx;
map{$sxx+=$x[$_]*$x[$_];}(0..$#x);
my $a=($sy*$sxx-$sx*$sxy)/($n*$sxx-($sx)**2);
my $b=($n*$sxy-$sx*$sy)/($n*$sxx-($sx)**2);
return "$a\t$b";
}

##加和
sub sum
{
my $sum;
map{$sum+=$_}@_;
return $sum;
}


###字符串动态比对结果的最佳结果
sub clt
{
my $q=shift @_;
my $max=shift @_;
my @s=@_;
my $r="NA";
for(@s)
	{
    my $query=$_;
	my @in=($q,$query);
	my @out=&aln(@in);
	my  $out=$out[0];
	if($out>$max)
		{
		$max=$out;
		$r=$query;
	    }
	elsif($out==$max)
		{
		$r.=",$query";	
	    }
    }

    my @o=($q,$r,$max);
	return @o;
}

###反向互补序列
sub rev
{
my $t=uc($_[0]);
my @s=split("",$t);
my %rev;
$rev{"A"}="T";
$rev{"T"}="A";
$rev{"C"}="G";
$rev{"G"}="C";
my $q="";
for(@s)
	{
	if(!$rev{$_})
		{
	     $rev{$_}=$_;
	    }
    $q="$rev{$_}"."$q";
   }
   return $q;
}

##翻译
my $cd="
GGG	Gly	G	Glycine
GGA	Gly	G	Glycine
GGC	Gly	G	Glycine
GGT	Gly	G	Glycine
GCG	Ala	A	Alanine
GCA	Ala	A	Alanine
GCC	Ala	A	Alanine
GCT	Ala	A	Alanine
GTA	Val	V	Valine
GTG	Val	V	Valine
GTC	Val	V	Valine
GTT	Val	V	Valine
TTG	Leu	L	Leucine
TTA	Leu	L	Leucine
CTG	Leu	L	Leucine
CTA	Leu	L	Leucine
CTC	Leu	L	Leucine
CTT	Leu	L	Leucine
ATA	Ile	I	Isoleucine
ATC	Ile	I	Isoleucine
ATT	Ile	I	Isoleucine
CCG	pro	P	Proline
CCA	pro	P	Proline
CCC	pro	P	Proline
CCT	pro	P	Proline
TTC	Phe	F	Phenylalanine
TTT	Phe	F	Phenylalanine
TAC	Tyr	Y	Tyrosine
TAT	Tyr	Y	Tyrosine
TGG	Trp	W	Tryptophan
TCG	Ser	S	Serine
TCA	Ser	S	Serine
TCC	Ser	S	Serine
TCT	Ser	S	Serine
AGC	Ser	S	Serine
AGT	Ser	S	Serine
ACG	Thr	T	Threonine
ACA	Thr	T	Threonine
ACC	Thr	T	Threonine
ACT	Thr	T	Threonine
TGC	Cys	C	Cysteine
TGT	Cys	C	Cysteine
ATG	Met	M	Methionine
AAC	Asn	N	Asparagine
AAT	Asn	N	Asparagine
CAA	Gln	Q	Glutamine
CAG	Gln	Q	Glutamine
GAC	Asp	D	Aspartic Acid
GAT	Asp	D	Aspartic Acid
GAG	Glu	E	Glutamate
GAA	Glu	E	Glutamate
AAG	Lys	K	Lysine
AAA	Lys	K	Lysine
CGG	Arg	R	Arginine
CGA	Arg	R	Arginine
CGC	Arg	R	Arginine
CGT	Arg	R	Arginine
AGG	Arg	R	Arginine
AGA	Arg	R	Arginine
CAC	His	H	Histidine
CAT	His	H	Histidine
TAA	e	e	e
TAG	e	e	e
TGA	e	e	e
---	.	.	.
";
my %amino;
my @cd=split(/\n/,$cd);
for(@cd)
	{
		next if(/^$/);
		my @a=split(/\t/,$_);
		$amino{$a[0]}=$a[2];
	}
my $base="
R	A/G
Y	C/T
M	A/C
K	G/T
S	G/C
W	A/T
H	A/T/C
B	G/T/C
V	G/A/C
D	G/A/T
N	A/T/C/G
A	A
T	T
C	C
G	G
";
my %base;
my @base=split(/\n/,$base);
for(@base)
	{
		next if(/^$/);
		my @a=split(/\t/,$_);
		$base{$a[0]}=$a[1];
	}

##翻译CDS，支持兼并碱基
sub mab
{
my @mab;
for(@_)
	{
    my @out;
	$_=uc($_);
	my @a=split("",$_);
    if(@a!=3)
		{
		push(@mab,"!");
        next;
	    }
    
    map{$_ ="N" if(!$base{$_});}@a;
	my @a1=split("\/",$base{$a[0]});
	my @a2=split("\/",$base{$a[1]});
	my @a3=split("\/",$base{$a[2]});
	
	for my $a1(@a1)
		{
		for my $a2 (@a2)
			{
			for my $a3(@a3)
				{
				push(@out,"$a1$a2$a3");
			    }
		    }
	    }
    my %c;
    map{$_=$amino{$_};$c{$_}++;}@out;
	@out=sort {$c{$b}<=>$c{$a}}keys %c;
	push (@mab,$out[0]);
	}
	return @mab;
}

##支持密码子表的翻译
sub tran
{
#####读取密码表
my @g=qw/A T C G/;
my @s=split(//,uc($_[0]));
my $out="";
if(($#s+1)%3 !=0)
	{
return "!$_[0]不是CDS\n";
last;
    }
else
	{
	my($m,$amino);
    for(0..$#s)
	  {
       $m++;
	   $amino.="$s[$_]";
	   if($_ ==2)
		  {
	   if($amino=~/^.TG$/ || $amino=~/^AT.$/)
		   {
		   $amino="ATG";
	       }
       else
		   {
		   return "!$_[0]不含编码框起始位置\n";
		   last;
		   }
	      }
	   if($m==3)
		{
		if(!$amino{$amino})
			{
			$amino{$amino}=(mab($amino))[0];
		    }

		last if($amino{$amino} eq "e");
		$out.=$amino{$amino};
		$m=0;
	    $amino="";
	    }
      }
return $out;
	}
}


###检查数据文件格式
sub cf
{
my @f=@_;
@f=map{$_=$1 if(/^\"(.+?)\"$/);$_}@f;
my @out;
for(@f)
	{
    my $out;
	next if(/^$/);
	if(!-s "$_")
		{
	print "file $_ no exists!\n";
	$out="NA";
	next;
	    }
	else
		{
	if(/\.gz$/)
		{
		open IN,"zcat $_|less|" or die "$!";
	    }
	else
		{
	open IN,"$_" or die "$!";
		}
	my $n=0;
	my (%for,%seq);
	my $id;
	while(<IN>)
			{
		chomp;
		next if(/^$/);
		last if($n>10);
		    if(/^\@(.+?)$/)
		  	   {
				$n++;
				$for{$n}="fq";
				if(/\s1\:/)
				   {
					$for{$n}="1.fq";
				   }
				elsif(/\s2\:/)
				   {
					$for{$n}="2.fq";
				   }
			   }
			elsif(/^\>(.+?)$/)
			   {
               $n++;
			   $for{$n}="fa";
			   }
			else
			   {
			   $seq{$n}.=$_;
			   }
			}
	close IN;
	my @t;
    for(1..10)
		  {
          next if(!$for{$_});
		  if($for{$_} eq "fa")
			  {
              $seq{$_}=uc($seq{$_});
			  if(length $seq{$_}>50000)
				  {
			  push(@t,"fna");
			  next;
			      }
			  my $pro=&tran($seq{$_});
			  unless($pro=~/^!/)
				  {
				  push(@t,"cds");
			      }
				  else
				  {
                  my $s=$seq{$_}=~s/[A|T|C|G]//g;
				  if($s>length($seq{$_}))
					  {
				  push(@t,"fna");
					  }
				  else
					  {
				  push(@t,"pro");
				      }
				  }
		      }
		  elsif($for{$_} =~/fq$/)
			  {
		      push(@t,"$for{$_}");
		      }
	      }
	  my %h;
	  map{$h{$_}++}@t;
      @t=sort {$h{$b}<=>$h{$a}}keys %h;
	  if(@t==1)
		{
		  push(@out,$t[0]);
	    }
		else
		{		  
		  push(@out,"fna");
		}
	 }
#push(@out,$out);
   }
   return @out;
}



sub panseq
{
my @t=@_;
my (@f,@e);
my @t1;
for(@t)
	{
	next if(/^$/);
	my @s=split(/\,/,$_);
    @t1=(@t1,@s);
	}
#去重复
@t1=&uniq(@t1);
for(@t1)
	{
	  my ($k1,$k2);
	  ($k1,$k2)=($1,$2)if(/^(\d+?)\-(\d+?)$/);
	   next if(!$k1||!$k2);
	   if($k1<$k2)
		    {
			push(@f,$k1);
		    push(@e,$k2);
		    }
	   else
			{
			push(@f,$k2);
			push(@e,$k1);
			}
	}
##去重复
my %f;
for(0..$#f)
	{
  	my $f=$f[$_];
	my $e=$e[$_];
	if(!$f{$f})
		{
		$f{$f}=$e;
	    }
	elsif($f{$f}<$e)
	    {
          $f{$f}=$e;
		}
	}

undef @t;
my($ft,$et)=("","");
for(sort{$a<=>$b} keys %f)
	{
	if(!$ft)
		{
		$ft=$_;
		$et=$f{$_};
	    }
    else
		{
	    if($_<=$et+1)
			{
			$et=$f{$_} if($f{$_}>$et)
			}
		else
			{
		    push(@t,"$ft\-$et");
		    ($ft,$et)=($_,$f{$_});
			}
		}
    }
push(@t,"$ft\-$et") if(!@t || $t[-1] ne "$ft\-$et");
return @t;
}

###输出样本队列 CoreCutoff 
sub coreseq
{
my $cut;
$cut=pop @_;
my @t=@_;
if($cut=~/(\d+)\-(\d+)/)
	{
	@t=(@t,$cut);
	$cut=@t;
	}
elsif($cut=~/^\d*\.*\d+$/)
	{
	$cut=int($cut*@t) if($cut<=1 &&$cut>=0);
    }
else
	{
	return "";
    }

my (%f,%e,%t);
my $n;
##样本
for(@t)
	{
    next if(/^$/);
	$n++;
	my @s=split(/\,/,$_);
	my (@f,@e);
	@s=panseq(@s);
##样本内位置
	for(@s)
		{
	    my ($k1,$k2);
	    ($k1,$k2)=($1,$2)if(/^(\d+?)\-(\d+?)$/);
	    next if(!$k1||!$k2);
		my($f,$e);
		if($k1<$k2)
		    {
			($f,$e)=($k1,$k2);
		    }
	    else
			{
            ($f,$e)=($k2,$k1);
			}
		$t{$f}{$n}=$e;
        $e{$e}++;
        $f{$f}++;
		}
	}
my @n=keys %f;
@n=(@n,keys %e);
@n=sort {$a<=>$b}uniq(@n);
my %c;
my $c;
$c=$f{$n[0]};
for(1..$#n)
	{
    if($f{$n[$_]})
		{
		$c{$n[$_-1]}{$n[$_]}+=$c;
		$c=$c+$f{$n[$_]};
        $c=$c-$e{$n[$_]} if($e{$n[$_]});
	    }
    elsif($e{$n[$_]})
		{
        $c{$n[$_-1]}{$n[$_]}+=$c;
		$c=$c-$e{$n[$_]};
	    }
	}
my @t1;
my ($av,$sum);
my %te;
my %tf;
my $stat;
for(sort {$a<=>$b}keys %c)
	{
my $k1=$_;
for(keys %{$c{$k1}})
		{
	    my $k2=$_;
		$stat.= "$k1\-$k2";
		for(1..@t)
			{
			if($t{$k1}{$_})
				{
				$te{$_}=$t{$k1}{$_};
				$tf{$_}=$k1;
				$stat.= "\t$c{$k1}{$k2}\|$tf{$_}\-$te{$_}";
			    }
				elsif($te{$_} && $te{$_}>=$k2)
				{
					$stat.= "\t$c{$k1}{$k2}\|$tf{$_}\-$te{$_}";
				    $te{$_}=() if($te{$_}==$k2);
				}
				else
				{$stat.= "\t0\|0";}
		    }
			$stat.= "\n";
		push(@t1,"$k1\-$k2")unless($c{$k1}{$k2}<$cut);
        $sum+=abs($k2-$k1);
		$av+=abs($k2-$k1) * $c{$k1}{$k2};
		}
    }
@t1=panseq(@t1);
if($sum>0)
	{$av=$av/$sum;}
else
	{$av=0;}
$av=$av/@t;
my $out=join("\,",@t1);
return ($out,$av,$stat);
}



1;
__END__
# Below is stub documentation for your module. You'd better edit it!

=head1 NAME

My - Perl extension for blah blah blah

=head1 SYNOPSIS

  use My;
  blah blah blah

=head1 DESCRIPTION

Stub documentation for My, created by h2xs. It looks like the
author of the extension was negligent enough to leave the stub
unedited.

Blah blah blah.

=head2 EXPORT

None by default.



=head1 SEE ALSO

Mention other useful documentation such as the documentation of
related modules or operating system documentation (such as man pages
in UNIX), or any relevant external documentation such as RFCs or
standards.

If you have a mailing list set up for your module, mention it here.

If you have a web site set up for your module, mention it here.

=head1 AUTHOR

A. U. Thor, E<lt>chenjing@(none)E<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2022 by A. U. Thor

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.26.2 or,
at your option, any later version of Perl 5 you may have available.


=cut

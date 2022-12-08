#!/workdir/staff/chenjing/Biosoft/bin/miniconda/bin/perl

############################################################
#	Copyright (C) 2020 matridx-hangzhou-researchgroup
#	Written by chenjing (cjwell@163.com)
#	Version: 2.0 (April 23th, 2020)
############################################################

#use strict;
use Getopt::Std;
use FindBin qw($Bin $RealBin $Script);
use List::Util qw/max min/;
use List::MoreUtils qw(firstidx);
use Cwd 'abs_path';
use Thread;
use lib "$RealBin/lib/";
use My qw/panseq coreseq/;

sub usage{
	print STDOUT "usage: $0  -i fasta文件路径 -r ref.fa -t -o ./outfile \n";
	print STDOUT "Options\n";
	print STDOUT "\t-i : 输入fasta文件夹路径，仅支持fna/ffn/fa/fasta后缀\n";
    print STDOUT "\t-o : ./SNPtree_out\n";
	print STDOUT "\t-r : 参考菌株的fasta文件路径 \n";
	print STDOUT "\t-m : 默认m，使用mummer（m）或者blat（b）计算CoreSNP并构建StrainTree，基因组比对推荐设置为m,core基因集或者小片段推荐设置为b\n";
	print STDOUT "\t-a : 仅支持blat方法，设置该参数将输出全部的比对碱基而不仅限于SNP，适合病毒序列比对，枝长更接近实际进化速率\n";
	print STDOUT "\t-f : 设置该参数将删除-o路径下已存在的结果文件,重新运行程序生成新的结果文件。
	可选参数r,m,a,t
	1=r=res文件(单样本SNP)
	2=m=merge_out文件(SNP合并结果)
	3=a=aln/matrix文件(根据-q参数计算得到的CoreSNP比对结果)
	4=t=tree文件(CoreSNP结果计算得到的最大似然tree)
\n";
	print STDOUT "\t-t : 设置该参数将aln文件转化为tree文件,并生成tree的可视化PDF文件\n";
	print STDOUT "\t-c : 配置文件路径 default:$Bin/straintree_config
	fasttree=path/fasttree
	Rscript=path/Rscript
	mummer:nucmer,show-snps,show-coords
	\n";
	print STDOUT "\t-l : 输入需要makre的区间，格式为 SeqID\t1-19\t200-221\n";
	print STDOUT "\t-p : 设置-l时生效，提取或者过滤maker区间，filter（f）或者pick(p)\n";
	print STDOUT "\t-d : 输入参考系路径或者使用预设参考系（未测试）
	1：预设covid2019参考系（$Bin/database/2019_covid）
	2：预设Coronaviridae参考系（$Bin/database/Coronaviridae）
	3：预设Influenza参考系（$Bin/database/Influenza）
	\n";
	print STDOUT "\t-n : 多线程运行程序,默认-n 5\n";
	print STDOUT "\t-q : 可设置0-1，大于1时根据样本数自动转换0-1范围，默认值1，表示只使用core为100%的SNP位置构建进化树，设置0.0时表示将使用所有碱基包括插入缺失的碱基，近亲建议设置为1，远亲设置为0.0\n";
	print STDOUT "\t-h : show the help message\n";
	exit(1);
}

our($opt_i,$opt_o,$opt_f,$opt_t,$opt_d,$opt_c,$opt_g,$opt_m,$opt_r,$opt_p,$opt_a,$opt_q,$opt_n,$opt_s,$opt_l,$opt_h);
getopts('i:o:d:f:t,c:r:p:g:m:a,q:n:s:l:h');

###检查输入
&usage if ($opt_h);
my ($in,$str);
unless((defined $opt_i)==1)
{
print "input file is no exist!!!please check -i parmeter!\n\n";
exit;
}
else
{
$in=(split("\/",$opt_i))[-1];
}

###获取Bin
if(-f $0 && -l $0)
{
$Bin=abs_path($0);
$Bin=`dirname $Bin`;
chomp($Bin);
}


####获取当前时间
my  ($sec,$min,$hour,$mday,$mon,$year) = (localtime)[0..5];
($sec,$min,$hour,$mday,$mon,$year) = (
    sprintf("%02d", $sec),
    sprintf("%02d", $min),
    sprintf("%02d", $hour),
    sprintf("%02d", $mday),
    sprintf("%02d", $mon + 1),
    $year + 1900
);
$nowtime="$year$mon$mday$hour$min";

###设定默认参数
###生成出输出文件夹
$opt_o ||="StrainTree_out";
$opt_o = abs_path ($opt_o);
$opt_m ||= "m";
$opt_q ||= 1;
$opt_s ||= 10;
$opt_n ||=5;
$opt_n=abs(int($opt_n));
$opt_n =1 if($opt_n<1);
my %ck;
my @q;
my %snp;
my %stat;
if($opt_m eq "m")
{
`mkdir -p $opt_o/mummer_out`;
}
elsif($opt_m eq "b")
{
`mkdir -p $opt_o/blat_out`;
}

$opt_c ||= "$Bin/straintree_config";
unless(-s "$opt_c")
{
print "请设置正确的配置文件$opt_c!\n\n";
exit;
}
if($opt_d eq "1")
{
$opt_d ="$Bin/database/2019_covid";
$opt_r||="$Bin/database/2019_covid/china-wuhan-EPI_ISL_406798-2019.fasta";
}
elsif($opt_d eq "2a")
{
$opt_d ="$Bin/database/Coronaviridae/Alpha";
$opt_r ||="$Bin/database/Coronaviridae/Alpha/Coronaviridae_Alphacoronavirus_Alphacoronavirus-1_79-1146_USA.fna"
}
elsif($opt_d eq "2b")
{
$opt_d ="$Bin/database/Coronaviridae/Beta";
$opt_r ||="$Bin/database/Coronaviridae/Beta/Coronaviridae_Betacoronavirus_Betacoronavirus-1__USA.fna"
}
elsif($opt_d eq "3")
{
$opt_d ="$Bin/database/Influenza";
$opt_r ||="$Bin/database/Influenza/I_A_H1N1_1_genomic.fna"
}
elsif(-f "$opt_d")
{
$opt_r ||=$opt_d;
}

####生成脚本日志
`mkdir -p $opt_o`;
open LOG,">$opt_o/run_log_$nowtime";

##前处理
########读取maker区间
my %p;
my %m;
if(-s "$opt_l")
{
open L,"$opt_l" or die "$!";
while(<L>)
	{
chomp;
next if(/^$/);
my @a=split(/\t/,$_);
$p{$a[0]}="all" if(@a==1);
next if($p{$a[0]} eq "all");
my $n=1;
while($n <= $#a)
	{
	if($a[$n]=~/(\d+)\-(\d+)$/)
		{
		$p{$a[0]}.=".$a[$n]";
	    $n++;
		}
	elsif($a[$n]=~/^\d+$/) 
		{
		if($a[$n+1] ne "")
			{$p{$a[0]}.=",$a[$n]\-$a[$n+1]";}
		else
			{$p{$a[0]}.=",$a[$n]\-$a[$n]";}	
		$n+=2;
		}
    else
		{$n++;}
    }
	}
close L;
##重新排序探针区间
for my $k(keys %p)
{
next if($p{$k} eq "all");
my @k;
@k=panseq($p{$k});
$p{$k}=join("\,",@k);
print  "Mark:\n$k\t$p{$k}\n";
}
}

#####是否删除结果文件
if($opt_f ne "")
{
my $c;
if($opt_f =~/[1|r]/)
{
$c.="rm $opt_o/*_out/*.res 2>/dev/null
";
}
if($opt_f =~/[2|m]/)
{
$c.="rm $opt_o/merge_out 2>/dev/null\n";
}
if($opt_f =~/[3|a]/)
{
$c.="rm $opt_o/CoreBase.aln 2>/dev/null
rm $opt_o/CoreBase.matrix 2>/dev/null
";
}
if($opt_f =~/[4|t]/)
{
$c.="rm $opt_o/CoreBase.tree* 2>/dev/null
";
}
print LOG $c;
system("$c");
}

###

my @in;
my @in1;
if(-d "$opt_i")
{
@in1=`ls $opt_i/* 2>/dev/null`;
if(-s "$opt_i/$opt_r")
	{
	$opt_r ="$opt_i/$opt_r";
	$opt_r=abs_path($opt_r);
	}
}
elsif(-f "$opt_i")
{
@in1[0]=$opt_i;
}
else
{
@in1=`ls -d $opt_i 2>/dev/null`;
}

if(@in1==0)
{
print "-i参数$opt_i不存在\n";
exit;
}
map{chomp;$in1[$_]=abs_path($in1[$_])}(0..$#in1);


my @in2;
if($opt_d)
{
@in2=`ls $opt_d/* 2>/dev/null`;
if(-s "$opt_d/$opt_r")
	{
	$opt_r ="$opt_d/$opt_r";
	$opt_r=abs_path($opt_r);
	}
}
@in=(@in1,@in2);
@in=grep{$_=~/ffn$|fa$|fna$|fasta$/}@in;
###############################################################read config
my %h;
open IN,"$opt_c";
while(<IN>)
{
chomp;
my $p;
next if(/^$/);
if(/^\s*(.+?)\s*=\s*(.+?)\s*$/)
	{
$p=$1;
$h{$p}=$2;
	}
else
	{
	$p=$_;
	}
chomp($h{$1}=`which $1`) if(!-s "$h{$p}");
if(!-s "$h{$p}")
	{
    print "请使用which $p确认$opt_c中的$p是否存在于PATH变量下，或者用$p=path/$p设置一个可用的路径\n";
	}
}
close IN;

##获取参考名
my $r=(split(/\//,$opt_r))[-1];
my $ref=$1 if($r=~/^(.+)\.(ffn|fa|fna|fasta)$/);

#获取输入队列
for(@in)
{
chomp;
my $q=(split("/",$_))[-1];
next if($q eq "$r");
$q=$1 if($q=~/^(.+)\.(ffn|fa|fna|fasta)$/i);
##记录query数目
push(@q,"$q");
push(@p,"$_");
}
@q=&uniname(@q);
for(0..$#p)
{
$q[$_]="$q[-2]$q[$_]\t$p[$_]";
}
pop @q;
pop @q;
##-q值根据输入统一为比例
if(abs($opt_q)>1)
	{
print "注意：\n-q 参数 $opt_q 绝对值值大于1\n";
$opt_q = (abs($opt_q)/@q);
print "使用比例值 $opt_q 进行core过滤\n";
	}


if($opt_m eq "m")
{
`mkdir -p $opt_o/mummer_shell/`;
open SH,">$opt_o/mummer_shell/multi.sh";
for(@q)
{
chomp;
my ($q,$path)=($1,$2)if(/^(.+?)\t(.+?)$/);
my $c;
$c= "$h{nucmer} --mum -o -p $opt_o/mummer_out/$q\_vs_$ref  $opt_r $path ;";
$c.="$h{'show-snps'} -C -r $opt_o/mummer_out/$q\_vs_$ref.delta > $opt_o/mummer_out/$q\_vs_$ref.snps \n";
print SH "$c" if(!-s "$opt_o/mummer_out/$q\_vs_$ref.delta");
}
close SH;
&msh("$opt_o/mummer_shell/",$opt_n,"$opt_o/mummer_shell/multi.sh");

for(@q)
	{
my ($q,$path)=($1,$2)if(/^(.+?)\t(.+?)$/);
#读取比对区间
open C,"$opt_o/mummer_out/$q\_vs_$ref.coords";
<C>;
<C>;
<C>;
<C>;
<C>;
while(<C>)
	{
	chomp;
	my @a=split(/[\t|\s]+/,$_,-1);
    ##读取区间
	if(!$ck{$a[-2]}{$q})
		{$ck{$a[-2]}{$q}="$a[1]\-$a[2]";}
	else
		{$ck{$a[-2]}{$q}.=",$a[1]\-$a[2]";}
	}
close C;

if(!-s "$opt_o/mummer_out/$q\_vs_$ref.res")
	{
open SNP,"$opt_o/mummer_out/$q\_vs_$ref.snps";
open RE,">$opt_o/mummer_out/$q\_vs_$ref.res";
my $n;
while(<SNP>)
	{
    chomp;
	$n++;
	next if($n<6);
	my @a=split;
	next if($a[1] eq "." || $a[2] eq ".");
    next if(length($a[1]) >1 || length($a[2]) >1);
	print RE "$a[-2]|$a[0]|$a[1]\t$a[2]\n";
   }
close SNP;
close RE;
   }
  }##队列@q结束

if(!-s "$opt_o/merge_out")
{
$c="$RealBin/bin/merge.pl -o $opt_o/merge_out -d NA -m chr -n -u e $opt_o/mummer_out/*.res \n";
system($c);
print LOG $c;
}

}##mummer结束
elsif($opt_m eq "b")#blat程序
{
`mkdir -p $opt_o/blat_shell`;
open SH,">$opt_o/blat_shell/multi.sh";
for(@q)
{
chomp;
my ($q,$path)=($1,$2)if(/^(.+?)\t(.+?)$/);
next if($q eq "$r");
my $c;
$c="$RealBin/bin/blat2nuc.pl -i $path -f $opt_r -n all -c 0.05 -s 0.05 -v -o $opt_o/blat_out -p $q\_vs_$ref \n";
print SH "$c" if(!-s "$opt_o/blat_out/$q\_vs_$ref.blastseq");
}
close SH;           

###执行多线程

&msh("$opt_o/blat_shell/",$opt_n,"$opt_o/blat_shell/multi.sh");

##

for(@q)
{
my ($q,$path)=($1,$2)if(/^(.+?)\t(.+?)$/);

##记录比对结果
if(!-s "$opt_o/blat_out/$q\_vs_$ref.res")
	{
##读取比对结果
open BS,"$opt_o/blat_out/$q\_vs_$ref.blastseq";
my %base;
while(<BS>)
	{
    chomp;
	next if(/^$/);
	my ($s,$qs,$ss);
	chomp($s=<BS>);
	my($n,$f,$e)=($1,$2,$3) if($s=~/^(.+?)\:(\d+?)\-(\d+?)$/);
    chomp($qs=<BS>);
	chomp($ss=<BS>);
	next if(!$f ||!$e);
	$qs=uc($qs);
	$ss=uc($ss);
	my @qs=split(//,$qs);
	my @ss=split(//,$ss);
    my $m1=0;
	if($f<$e)
	{
	for(0..$#ss)
		{
        next if($ss[$_] eq "-");
		my $m=$f+$m1;
##设定-v时 只保留SNP位置
         if($opt_a)
			{
            $base{"$n\|$m\|$ss[$_]"}=$qs[$_];
		    }
		 elsif($ss[$_] ne "$qs[$_]")
			{
			$base{"$n\|$m\|$ss[$_]"}=$qs[$_];
			}
		$m1++;
		}
	}
	else
	{
		for(0..$#ss)
		{
        next if($ss[$_] eq "-");
		my $m=$f-$m1;
##设定-v时 只保留SNP位置
         if($opt_a)
			{
            $base{"$n\|$m\|$ss[$_]"}=$qs[$_];
		    }
		 elsif($ss[$_] ne "$qs[$_]" && $ss[$_]!~/^[X|N]$/ && $qs[$_]!~/^[X|N]$/)
			{
			$base{"$n\|$m\|$ss[$_]"}=$qs[$_];
			}
		$m1++;
		}
	}
	<BS>;
    }
close BS;

#记录比对结果
open RE,">$opt_o/blat_out/$q\_vs_$ref.res";
for(sort keys %base)
	{
	print RE "$_\t$base{$_}\n";
    }
close RE;
	}

#读取比对区间
open BL,"$opt_o/blat_out/$q\_vs_$ref.blastout";
<BL>;
while(<BL>)
	{
	chomp;
	my @a=split(/\t/,$_,-1);
    ##读取区间
	if(!$ck{$a[1]}{$q})
		{$ck{$a[1]}{$q}=$a[-1];}
	else
		{$ck{$a[1]}{$q}.=",$a[-1]";}
	##计算cov和ID
	$stat{$q}{$ref}{id}+=$a[2]*$a[4];
    $stat{$q}{$ref}{qal}+=$a[4];
	}
close BL;


}#@in结束

##整合比对结果
if(!-s "$opt_o/merge_out")
{
$c="$RealBin/bin/merge.pl -o $opt_o/merge_out -d NA -t n -m chr -n -a -u e $opt_o/blat_out/*.res\n";
system($c);
print LOG $c;
}
} #blat结束

##记录core位点
`mkdir -p $opt_o/core_site 2>/dev/null`;
my %core;
open STAT,">$opt_o/core.stat";
print STAT "RefID\tAlnRatio\tAverageCoreRatio\tCutOff\tCoreSite\n";
for(sort keys %ck)
	{
	my $k1=$_;
	my @t;
	my $site;
	my @cq=sort keys %{$ck{$k1}};
    for(@cq)
		{
		my $k2=$_;
		my $t;
		$t=$ck{$k1}{$k2};
	    push(@t,$t);
		}
	my $tq=@t/@q;
	my $fq=0;
	if($tq <$opt_q)
		{
	    $core{$k1}="0-0";
	    }
	else
		{
###core的比例0-1
	@t=coreseq((@t,$opt_q));
	$t[0]="0-0" if($t[0] eq "-");
	$core{$k1}=$t[0];
	$fq=$t[1];
	$site=$t[2];
		}
    print STAT "$k1\t$tq\t$fq\t$opt_q\t$core{$k1}\n";
##输出同源位点
	open SITE,">$opt_o/core_site/$k1.coresite";
    map{print SITE "\t$_"}@cq;
	print SITE "\n$site\n";
	close SITE;
	}
close STAT;

####ANI值
for(sort keys %stat)
{
my $id=$stat{$_}{$ref}{id}/$stat{$_}{$ref}{qal};
print "$_\t$id\n";
}


if(!-s "$opt_o/CoreBase.aln")
{
open M,"$opt_o/merge_out";
open O,">$opt_o/CoreBase.matrix";
my %seq;
my $h;
chomp($h=<M>);
my @h=split(/\t/,$h,-1);
$h[0]=$ref;
print O "$h\n";
while(<M>)
{
chomp;
my @a=split("\t",$_);
my @ref=split(/\|/,$a[0]);
next if(@ref<3);
$ref[0]=join("\|",@ref[0..$#ref-2]) if(@ref>3);

###检查位点是否为Core位点
next if(!$core{$ref[0]});

my $c=&check(($core{$ref[0]},$ref[-2],$ref[-1]));
next if($c eq "-");

$seq{$h[0]}.=$ref[-1];

##检查整合位置的比对信息
for(1..$#a)
	{
	if($a[$_] eq "NA")
		{
		$a[$_]=&check(($ck{$ref[0]}{$h[$_]},$ref[-2],$ref[-1]));
        }
    $seq{$h[$_]}.=$a[$_];
	$snp{$h[$_]}++ if($a[$_]=~/^[A|T|C|G]$/);
	}


print O join("\t",@a);
print O "\n";
}
close M;

$snp{$h[0]}=length($seq{$h[0]});

open ALN,">$opt_o/CoreBase.aln";
for(0..$#h)
	{
       if($snp{$h[$_]} >= $opt_s)
		{
		print ALN ">$h[$_]\n$seq{$h[$_]}\n";
		}
		else
		{
		print "$h[$_] COV不足或亲缘关系过远，与ref的同源比对碱基数目$snp{$h[$_]}低于设定最低个数$opt_s！\n";
		}
	}
close ALN;
}


$c="";
if($opt_t)
{
$c.="$h{fasttree} -nt -gtr < $opt_o/CoreBase.aln > $opt_o/CoreBase.tree \n" if(!-s "$opt_o/CoreBase.tree");
$c.="$h{Rscript} $RealBin/bin/ggtree_v3.2.R -i $opt_o/CoreBase.tree -g $ref \n" if(!-s "$opt_o/CoreBase.tree.pdf");
}
print LOG $c;
system($c) if($c ne "");




sub check
{
my ($nc,$n,$d)= @_[0..2];
my $check=0;
my @s=split("\,",$nc);
for(@s)
	{
	my ($f,$e);
	($f,$e)=($1,$2) if(/^(\d+)\-(\d+)$/);
    next if(!$f || !$e);
	if($n>= $f && $n <= $e)
		{
		$check=1;
		last;
		}
	}
if($check==1)
	{return $d;}
else
	{return "-";}
}

##多线程比对
sub thread_sh
{
my @sp=@_;
my @tsp;
for(@sp)
	{
	push @tsp,threads->create(\&sh,$_);
	}
for(@tsp)
	{
	$_->join();
	}
system(sleep 3);
print  "多线程任务结束！\n";
}
sub sh
{
system("sh $_ 2>/dev/null");
}

##输出路径 线程数 记录SH脚本队列的文本
sub msh
{
my ($o,$n,$m)=@_;
if(!-s "$m")
	{
	return "Error!\n";
    }
my @muti=`less $m 2>/dev/null`;
$n=abs(int($n));
$n =1 if($n<1);
my $l=length $n;
`rm $o/split_*.sh 2>/dev/null`;
while(@muti)
	{
    chomp;
	for(1..$n)
		{
		last if(@muti==0);
		$_=sprintf("%0".$l."d",$_);
		open SP,">>$o/split_$_.sh";
		my $sp=pop @muti;
	    print SP $sp;
		close SP;
		}
	}
##
#执行多线程
my @sp;
@sp=`ls -d $o/split_*.sh 2>/dev/null`;
map{chomp;}@sp;
if(@sp>0)
	{
&thread_sh(@sp);
	}
##
}

sub uniname
{
my @in=@_;
map{chomp}@in;
map{$_=(split("\/",$_))[-1];$_}@in;
my ($f,$e)=("","");
my $n;
for(@in)
	{
	next if(/^$/);
	$n++;
	my @f;
	my @e;
	my @str=split(//,$_);
	if($n ==1)
		{
		@f=@str;
		@e=@str;
		}
	else
		{
		@f=split(//,$f);
		@e=split(//,$e);
		}
	$f="";
	$e="";
	for(0..$#f)
		{
		if($f[$_] eq $str[$_])
			{$f=$f."$f[$_]";}
		else
			{last;}
	    }
	for((1..$#e+1))
		{
		if($e[-$_] eq $str[-$_])
			{$e="$e[-$_]".$e;}
		else
			{last;}
	    }
    }
my @out;
if(@in==1)
	{
	$f="";
	$e="";
    }
@out=map{$_ =$1 if(/^$f(.+?)$e$/)}@in;
@out=(@out,$f,$e);
return @out;
}

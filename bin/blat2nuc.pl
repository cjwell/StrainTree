#!/workdir/staff/chenjing/Biosoft/bin/miniconda/bin/perl

############################################################
#	Copyright (C) 2020 matridx-hangzhou-researchgroup
#	Written by chenjing (cjwell@163.com)
#	Version: 2.0 (April 23th, 2020)
############################################################

#use strict;
use Getopt::Std;
use FindBin qw($Bin $RealBin $Script);
use List::Util qw/max min uniq/;
use List::MoreUtils qw(firstidx);
use Cwd 'abs_path';


sub usage{
	print STDOUT "usage: $0  -i 输入文件 -d 1 -c 10 -s 50 -q 99\n";
	print STDOUT "当前路径下需要有blat程序，Options \n";
	print STDOUT "！输入参数\n";
	print STDOUT "\t-i : 输入比对用的核酸fasta格式文件\n";
	print STDOUT "\t-m : -i 输入的类型 dna or rna 默认dna\n";
	print STDOUT "\t-d : 比对预设数据库，如果需要自定义，请使用-f参数
	预设card.rRNA，card.cds，vfdb.setA，vfdb.setB，fungi,SILVA等核酸数据库
	1=card.rRNA	耐药数据库，rRNA序列
	2=card.cds 耐药数据库，cds核酸序列
	3=vfdb.setA 毒力数据库，核心序列
	4=vfdb.setB 毒力数据库，全部序列
	5=fungi_18s 真菌物种标签的数据库，默认-s=50
	6=SILVA_SSU 细菌物种标签的rRNA小亚基数据库，默认-s=50
	7=SILVA_LSU 细菌物种标签的rRNA大亚基数据库，默认-s=50
	8=AD Illumina 标准接头数据库
	9=metaphlan 对应版本为mpa_v30_CHOCOPhlAn_201901 物种标签的maker基因数据库
	10=mOTU 对应版本为CAMI_Sep_2015_mOTU 物种标签的maker基因数据库
	11=influenza 流感数据库
	\n";
	print STDOUT "\t-f : 指定参考核酸fasta文件作为临时数据库\n";
	print STDOUT "！可选参数\n";
	print STDOUT "\t-v : 设定此参数，将生成snps文件\n";	
    print STDOUT "\t-p : 输出文件名前缀\n";
	print STDOUT "\t-o : 输出文件夹，默认当前路径 blat_output\n";
    print STDOUT "\t-s : subject的coverage阈值,0-100默认10\n";
    print STDOUT "\t-c : query的coverage阈值，0-100默认0\n";
	print STDOUT "\t-q : Identity的阈值，默认90\n";
	print STDOUT "\t-n : 输出匹配数，可选1~n，设定为all表示全部阳性，默认1\n";
	print STDOUT "\t-a : 设定该参数将输出匹配成功的结果文件\n";
	print STDOUT "\t-h : show the help message\n";
	exit(1);
}

our($opt_i,$opt_o,$opt_l,$opt_c,$opt_d,$opt_b,$opt_g,$opt_q,$opt_s,$opt_t,$opt_f,$opt_k,$opt_p,$opt_r,$opt_n,$opt_v,$opt_a,$opt_m,$opt_u,$opt_h);
getopts('i:o:l:d:c:b:g:s:t:f:p:k:r:n:q:v,a,u:m:h');



###检查输入
&usage if ($opt_h);
my $str;
unless((defined $opt_i)==1)
{
print "please check -i parmeter!\n\n";
exit;
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
####生成脚本日志
open LOG,">$opt_o/run_log_$nowtime";

###设定默认参数
###生成出输出文件夹
$opt_o ||="blat_output";
`mkdir -p $opt_o`;
$opt_o = abs_path ($opt_o);
$opt_n ||= 1;

if($opt_d ==1 ||$opt_d eq "card.rRNA")
{
$opt_f ="$Bin/../database/CARD/card.rRNA.fa";
$opt_v ="y";
}
elsif($opt_d ==2||$opt_d eq "card.cds")
{
$opt_f ="$Bin/../database/CARD/card.cds.fa";
}
elsif($opt_d ==3 ||$opt_d eq "vfdb.setA")
{
$opt_f ="$Bin/../database/VFDB/VFDB_setA_nt_core.fas";
}
elsif($opt_d ==4 ||$opt_d eq "vfdb.setB")
{
$opt_f ="$Bin/../database/VFDB/VFDB_setB_nt_all.fas";
}
elsif($opt_d ==5||$opt_d eq "fungi_18s")
{
$opt_f ="$Bin/../database/Fungi/UNITE_Fungi";
$opt_s ||=50; 
}
elsif($opt_d ==6||$opt_d eq "SSU")
{
$opt_f ="$Bin/../database/SILVA_SSU/SILVA_SSURef_NR99.fa";
$opt_s ||=50; 
}
elsif($opt_d ==7||$opt_d eq "LSU")
{
$opt_f ="$Bin/../database/SILVA_LSU/SILVA_LSURef_NR99.fa";
$opt_s ||=50; 
}
elsif($opt_d ==8||$opt_d eq "AD")
{
$opt_f ="$Bin/../database/adapter/illumina_ad.fa";
$opt_s ||=50; 
}
elsif($opt_d ==9||$opt_d eq "metaphlan")
{
$opt_f ="$Bin/../database/metaphlan/metaphlan.fa";
$opt_s ||=50; 
}
elsif($opt_d ==10||$opt_d eq "mOTU")
{
$opt_f ="$Bin/../database/mOTU/mOTU.v2b.fa";
$opt_s ||=50; 
}
elsif($opt_d ==11||$opt_d eq "inf")
{
$opt_f ="$Bin/../database/influenza/influenza.fa";
$opt_s ||=50; 
$opt_c ||=50;
}


my %db;
my %aln;
if(-s "$Bin/../database/taxa.length")
{
open DB,"$Bin/../database/taxa.length" or die;
while(<DB>)
	{
chomp;
my @a=split(/\t/,$_);
$db{"$a[0]\t$a[1]"}{n}=$a[2];
$db{"$a[0]\t$a[1]"}{l}=$a[3];
	}
}

$opt_m ||="dna";
$opt_q||=90;
$opt_t ||= "e";

my $command;

####开始运行
#print LOG "###input seq is $opt_i !\n";



###-f 输入 blat比对的cds文件
if ((defined $opt_f)==1)
{
$opt_f=abs_path ($opt_f);
$str=(split("\/",$opt_f))[-1];
#$str=$1 if($str=~/^(.+)\./);
}
elsif(!$opt_d)
{
print "-f -d参数必须选择一个！！\n";
exit;
}


my @in;
if(-f "$opt_i")
{
$in[0]=$opt_i;
}
elsif(-d "$opt_i")
{
@in =`ls -d $opt_i/* 2>/dev/null`;
@in =`ls -d $opt_i/*/ 2>/dev/null` if(@in==0);
@in =grep{chomp;/[fasta|fa|fna|ffn|faa|cds]$/}@in;
print "$opt_i 路径下没有 fasta|fa|fna|ffn|faa|cds 后缀名的fasta文件\n" and exit if(@in==0);
}
elsif(!-s "$opt_i")
{
print "$opt_i 不存在，尝试匹配，如果输入中存在*匹配符，请用-i \"*\" 格式输入\n";
@in =`ls -d $opt_i 2>/dev/null`;
@in =grep{chomp;-f "$_"}@in;
print "运行ls -d $opt_i 没有成功匹配到文件\n" and exit if(@in==0);
}

###输入名字修剪
my @in1=@in;
@in1=&uniname(@in);
if($opt_t eq "f")
{
map{$in1[$_].=$in1[-1];}(0..$#in);
}
elsif($opt_t eq "e")
{
map{$in1[$_]=$in1[-2].$in1[$_];}(0..$#in);
}
elsif($opt_t eq "all")
{
map{$in1[$_]=$in1[-2].$in1[$_].$in[-1];}(0..$#in);
}


if(@in==1)
{
$in1[0]=$1 if($in1=~/^(.+)\./);
}

for(0..$#in)
{
chomp;
my $in=$in[$_];
if($#in>0)
	{
$opt_p = $in1[$_]."_blat2_".$str;
	}
else
	{
$opt_p ||= $in1[$_]."_blat2_".$str;
    }
#####生成blat脚本
my $c="";
if(!-s "$opt_o/$opt_p\.blast")
{
if($opt_m eq "dna")
	{
$c="$Bin/blat  $opt_f  $in  $opt_o/$opt_p\.blast -out=blast\n";
	}
	elsif ($opt_m eq "rna")
	{
$c="$Bin/blat  $opt_f  $in  $opt_o/$opt_p\.blast -out=blast -q=rna\n";	
	}
}



if($c ne "")
{
print LOG "#生成$opt_o/$opt_p\.blast文件！\n$c\n";
system("$c");
$c="";
}


###执行blat信息过滤
open O,">$opt_o/$opt_p.blastout";
print O "Query\tSubject\tIdentity\tQuery_Length\tQAL\tQuery_Coverage\tSubject_Length\tSAL\tSubject_Coverage\tDistance\tQuery_site\tSuject_site\n";
open ALN,">$opt_o/$opt_p.blastseq";

if($opt_v)
{
open SNP,">$opt_o/$opt_p.snpsite";
print SNP "SubjectRef\tQueryBase\n";
open V,">$opt_o/$opt_p.snps";
print V "Query\tSubject\tSNP\n";
my %snp1;
my %snp2;
}

my $inblast=`less $opt_o/$opt_p\.blast 2>/dev/null`;
my @a=split(/\nReference/,$inblast);

for(@a)
{
	##一个query的结果
	my @b=split(/\n>/,$_);
    my($query,$subject,$ql,$sl,$aln);
	my (%st,%qt);
	if($b[0]=~/Query=\s*(.+?)\n\s*\((\d+) letters\)\n/)
	{
	##Query
	$query=$1;
	$ql=$2;
	}
	###输出匹配数
	my $n;
	if($opt_n eq "all" || $opt_n>$#b)
	{
	$n = $#b;
	}
	else
	{
	$n=$opt_n;
	}
    my %core;##片段重复检查
	for(1..$n)
	{
	###最佳匹配$b[1]
    if($b[$_]=~/^(.+?)\s*\n\s*Length = (\d+?)\n(.+?)$/s)
	{
	##最佳subject
	$subject=$1;
	$sl=$2;

	###记录aln
	$aln=$3;


    ###最佳匹配的所有片段
    my($sal,$qal,$nm);
	my @c=split(/\n\sScore/,$aln);
	next if($#c==0);
    
    ###每次循环一个片段
	for(1..$#c)
		{
	    my @d=split(/\n/,$c[$_]);
        my($q,$s,$first,$end,$qf,$qe);

	    for(@d)
		  {
	       next if(/^$/);
	       if(/^Query:\s+(\d+?)\s+(.+?)\s(\d+)$/)
			{
	        $qf=$1 if($qf eq "");
			$q.=$2;
			$qe=$3;
	        }
	       elsif(/^Sbjct:\s+(\d+?)\s+(.+?)\s(\d+)$/)
			{
			$first=$1 if($first eq "");
			$s.=$2;
			$end=$3;
	        }
		   elsif(/Identities = (\d+)\/(\d+) \(/)
			{
			$qal+=$1;
			$sal+=$2;
			}
		  }
       next if($sal<50);
        ###反向匹配时反转互补
        if($first > $end)
		  {
           $s=~ tr/atcguATCGU/tagcaTAGCA/;
		   $s= scalar reverse $s;
		   $q=~ tr/atcguATCGU/tagcaTAGCA/;
		   $q= scalar reverse $q;
		   my $tmp=$first;
		   $first=$end;
		   $end=$tmp;
           $tmp=$qf;
		   $qf=$qe;
		   $qe=$tmp;
		  }
         ####检查重复匹配
		 my @tpan;
         if(!$core{$subject})
			{$core{$subject}="$first-$end";}
		 else
			{
			 @tpan=&panseq($core{$subject},"$first-$end");
		     my $tp=join("\,",@tpan);
			 if($tp eq $core{$subject})
				{
				 next;##跳过完全重叠的匹配
			    }
             else
				{$core{$subject}=$tp;}
			}

		 ##记录比对
		 $qt{$query}{$subject}.=",$qf-$qe";
		 $st{$query}{$subject}.=",$first-$end";
         $aln{$query}{$subject}.="$query:$qf-$qe\n$subject:$first-$end\n$q\n$s\n\n";

         ###最佳匹配的突变输出
         if($opt_v)
         {
         $q=uc($q);
		 $s=uc($s);
		 my @q=split(//,$q);
         my @s=split(//,$s);
         my $n=$first;
			 for(0..$#s)
			 {
			 ###输出突变
			 if($s[$_] ne $q[$_])
				 {
                  if($q[$_] !~/^[X|N]$/ && $s[$_] !~/^[X|N]$/)
					 {
					   ##格式1
					 $snp1{"$query$subject"}.="$subject\|$n\|$s[$_]\t$q[$_]\n" if($s[$_] ne "-" && $q[$_] ne "-"); #忽略插入缺失
					 if($s[$_] eq "-")
					   {$s[$_]="ins";}
				     elsif($q[$_] eq "-")
					   {$q[$_] = "del";}
				     my $snp="$s[$_]$n$q[$_]";
					  ##格式2
					 $snp2{"$query$subject"}.="$query\t$subject\t$snp\n";
					}
					else
					{
					$nm++;
					}
				 }
				 $n++ if($s[$_] ne "ins");
			 }
		 }
		}
my (@qout,@sout);
@qout=split("\,",$qt{$query}{$subject});
@sout=split("\,",$st{$query}{$subject});
@qout=&panseq(@qout);
@sout=&panseq(@sout);


##检查错误
        if($ql==0 ||$sl==0||$sal ==0)
		{
print LOG "$query 和 $subject 的比对中存在错误！\n";
next;
        }

#####输出比对结果

$sal=$sal-$nm;
my $id=100*$qal/$sal;

$qal=0;
map{$qal+=$2-$1+1 if(/^(\d+)-(\d+)$/)}@qout;
$sal=0;
map{$sal+=$2-$1+1 if(/^(\d+)-(\d+)$/)}@sout;
my $dis;
my $qc=100*$qal/$ql;
my $sc=100*$sal/$sl;

for((@qout,@sout))
		{
        my @sp=split("\-",$_);
		$dis+=abs($sp[1]-$sp[0])+1;
		};
$dis=100*$dis/($sl+$ql);

if($sc>$opt_s && $qc >$opt_c && $id> $opt_q)
		{
print  O "$query\t$subject\t$id\t$ql\t$qal\t$qc\t$sl\t$sal\t$sc\t$dis";
print  O "\t".join("\,",@qout);
print  O "\t".join("\,",@sout);
print  O "\n";
print ALN "$aln{$query}{$subject}\n";
		}
	}
	}#for结束

if($opt_v)
   	{
     my @snp;
	 @snp=split(/\n/,$snp1{"$query$subject"});
	 @snp=grep{next if($_ eq "");++$c{$_}==1}@snp;
	 print  SNP join("\n",sort @snp);
	 print SNP "\n";
	 @snp=split(/\n/,$snp2{"$query$subject"});
	 @snp=grep{next if($_ eq "");++$c{$_}==1}@snp;
	 print V join("\n",sort @snp);
	 print V "\n";
	}
}


close O;
if($opt_v)
	{
close SNP;
close V;
    }


###提取文件
if($opt_a)
{
 $c.="$Bin/pick_seq.pl  -i $in -f 0 -g $opt_o/$opt_p.blastout  -s >$opt_o/$opt_p\.query.fa\n";
 #$c.="$Bin/pick_seq.pl  -i $opt_f  -f 1 -g $opt_o/$opt_p.blastout  -s >$opt_o/$opt_p\.subject.fa\n";
 if($c ne "")
  {
   print LOG "#生成$opt_o/$opt_p\.query.fa文件！\n$c\n";
   system("$c") if(!-s "$opt_o/$opt_p.query.fa");
   $c="";
  }
close ALN;
}
close LOG;

##@in结束

if($opt_d ==7||$opt_d eq "LSU"||$opt_d ==9||$opt_d==10)
	{
my @out=&stat_BlatTaxa("$opt_o/$opt_p.blastout");
open O,">$opt_o/$opt_p.taxa";
print O join("\n",@out);
close O;
    }

}


##############子程序

sub stat_BlatTaxa
{
open TMP,"$_[0]" or die;
<TMP>;
my %sum;
my %q;
my %n;
my %dbl;
my $t;
if($opt_d ==9)
	{
$t="Metaphaln";
	}
elsif($opt_d ==10)
	{
$t="mOTU";
	}
elsif($opt_d ==7)
	{
$t="LSU";
	}
my @pick;
while(<TMP>)
	{
	chomp;
	next if(/^$/);
	my @a=split(/\t/,$_);
	my @b;
	@b=split(/\|/,$a[1]);
	my $g;
	$g=$1 if($b[0]=~/^(.+?)\;/);
	if(@b==1)
		{
	@b=split(/\;/,$a[1]);
	$g=$1 if($b[0]=~/^(.+?)\_/);
	print "数据库格式应为LSU\n" if($t ne "LSU");
		}
	$b[-1]=$1 if($b[-1]=~/^\s*(.+?)\s*$/);
	$sum{$b[-1]}+=$a[4];
	$q{$b[-1]}+=$a[2]*$a[4];
	$n{$b[-1]}++;
    $dbl{$b[-1]}{$g}+=$a[4] if(!$dbl{$b[-1]}{$g} | $a[4] > $dbl{$b[-1]}{$g});
	}
close TMP;

if(!%n)
	{
return "$t\tNA\tNA\tNA\n";
    }
for(sort keys %n)
	{
	my $k=$_;
	next if(/^:/);
	$k=$1 if($k=~/^\s*(.+?)\s*$/);
	next if($sum{$k}==0);
    my $q;
	$q=$q{$k}/$sum{$k};
	$q{$k}=$q;
 	my $db="$t\t$k"; 
	next unless($db{$db}{n} && $db{$db}{l});
	my($dn,$dl)=($db{$db}{n},$db{$db}{l});
	$out="$t\t$k\t$q{$k}\t$sum{$k}\t$n{$k}\t$dl\t$dn";
	my $dbn=keys %{$dbl{$k}};
    my $dbl;
	map{$dbl+=$dbl{$k}{$_}}keys %{$dbl{$k}};
	$dbn=100*$dbn/$dn;
    $dbl=100*$dbl/$dl;
	$out.="\t$dbl\t$dbn";
	push(@out,$out);
	}

	return @out;
}

sub uniname
{
my @in=@_;
chomp @in;
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
@out=map{chomp;$_ =$1 if(/^$f(.+?)$e$/);$_}@in;
@out=(@out,$f,$e);
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
@t1=uniq(@t1);
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
for(@t)
	{
    next if(/^$/);
	$n++;
	my @s=split(/\,/,$_);
	my (@f,@e);
	@s=panseq(@s);
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
undef @t;
for(sort {$a<=>$b}keys %c)
	{
my $k1=$_;
for(keys %{$c{$k1}})
		{
	    push(@t,"$k1\-$_")unless($c{$k1}{$_}<$cut);
        }
    }
@t=panseq(@t);


return @t;
}
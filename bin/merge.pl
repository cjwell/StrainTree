#!/usr/bin/env perl

############################################################
#	Copyright (C) 2020 matridx-hangzhou-researchgroup
#	Written by chenjing (cjwell@163.com)
#	Version: 2.1 (April 11th, 2022)
############################################################

#use strict;
use Getopt::Std;
use FindBin;
use List::Util qw/max min uniq/;
use Cwd 'abs_path';
use lib "$FindBin::Bin/lib";

sub usage{
	print STDOUT "usage: $0  -i samplelist\nOR\n\t$0  *inputfiles\n";
	print STDOUT "Options\n";
	print STDOUT "\t-i : 输入需要整合的文件list 路径，没有输入 -i参数时，默认输入的是整合文件名，整合参数为：-o meger_out -t n -r n -c n\n";
	print STDOUT "\t-s : 输出文件进行排序，可选c,r,n，c对行排序，r 对列排序，n不排序 \"r,c\" 表示对行列排序  default r\n";
	print STDOUT "\t-t : 文件是否包含样本抬头，y代表有抬头为样本名，n代表无样本名抬头，无抬头文件将以\"文件名_列数\"命名sample，default n\n";
    print STDOUT "\t-n : 设定该参数将使用文件名作为抬头的前缀\n";
	print STDOUT "\t-d : 当内容缺省时的default值，默认0.0\n";
	print STDOUT "\t-r : 检查feature是否有重复名字，y:重复名字将加入后缀_dup进行区分,n:不检查重复，此时相同feature将根据-m参数合并计算，默认n\n";
	print STDOUT "\t-c : 检查sample是否有重复名字，y:重复名字将加入后缀_dup进行区分,n:不检查重复，相同sample将根据-m参数进行合并计算，默认n\n";
	print STDOUT "\t-v : 选取特定列作为value值，从0起始，逗号\",\"分割,如 \"1，2\";默认设定为all，除了-f以外的所有列；\n";
	print STDOUT "\t-f : 选取特定的列作为feature，从0起始，逗号\",\"分割,如 \"0，1\"的特定复数列作为整合文件的feature名,默认第1列（0）\n";;
	print STDOUT "\t-p : 设定此参数时，相同文件中复数的feature和value将合并为1列\n";
	print STDOUT "\t-m : 文件内容为 num数字（相同feature/sample进行加和） 还是 chr字符串（相同feature/sample进行合并，\"；\"号链接） default: num\n";
	print STDOUT "\t-x : 过滤含向前匹配特殊符号的行，比如\"#\"抬头，用\n";
	print STDOUT "\t-g : 设定此参数时，将忽略以逗号\",\"分割,如 \"1,2\"的特定行数,行数从1起始\n";
    print STDOUT "\t-o : 输出文件夹，默认当前路径merge_out\n";
    print STDOUT "\t-u : merge_out,抬头去掉前后共同字符（c），去掉前端（f），去掉后端（e），不做任何处理（n）；默认：n\n";
    print STDOUT "\t-a : 设定此参数，将保留空文件\n";
	print STDOUT "\t-h : show the help message\n";
	exit(1);
}

our($opt_i,$opt_o,$opt_d,$opt_s,$opt_t,$opt_n,$opt_f,$opt_c,$opt_r,$opt_m,$opt_x,$opt_p,$opt_g,$opt_v,$opt_u,$opt_h);
getopts('i:o:d:s:t:f:c:r:m:x:p,g:v:u:n,a,h');

&usage if ($opt_h);
&usage if (@ARGV ==0 && (defined $opt_i)==0);


##设置默认参数
my $tmp_list=0;
if ((defined $opt_i)==0)
{
	$tmp_list=rand(10000);
    open LIST,">$tmp_list.tmplist";
    my @tmp_list=@ARGV;
    for(@tmp_list)
	{
	 chomp;
	 next if(/^$/);
	 next if(!-s "$_" && !$opt_a);
	 print LIST "$_\n";
    }
$opt_i ="$tmp_list.tmplist";
close LIST;
}
$opt_i= abs_path ($opt_i);
$opt_o ||="./merge_out";
$opt_o = abs_path ($opt_o);

$opt_f||="0";
my @f1;
@f1 =(split(/\,/,$opt_f));
my %gv;
map{$gv{$_}=1}@f1;
my @v =(split(/\,/,$opt_v));

my %gn;
if($opt_g)
{
for((split(/\,/,$opt_g)))
	{
	$gn{$_}=$_;
    }
}
$opt_t ||= "n" ;
$opt_c ||= "n";
$opt_r||="n";
$opt_s||="r";
$opt_m ||= "chr";
$opt_u||="n";
if($opt_m eq "num")
{
$opt_d ||= "0.0" ;
}
elsif($opt_m eq "chr")
{
$opt_d ||="";
}
else
{
print "-m 参数只能是num或者chr！\n";
exit;
}
##排序
my @s=split(/\,/,$opt_s);

open IN,"$opt_i" || die $!;
my $head;
my %m;
my @list;
my (@H,@F,@r,@c);
my (%V,%v);
my $nid;
my $hn=-1;
while(<IN>)
{
chomp;
next if($_ eq "");
my $id=(split(/\//,$_))[-1];

###检查文件是否存在
if(!-s "$_")
	{
    if($opt_a)
		{
	$c[$#c+1]=@c+1;
	my $ftn=$f1[0]+1;
	$H[$#H+1]="$id\_ValueColumn\_$ftn";
		}
		else
		{next;}
    }

open IN1,"$_" or die "$!";
my $r=0;
my (@h,@f);

##忽略的行
my $gn=0;

 while(<IN1>)
   {
   chomp;
   ##忽略特定行
   if($opt_g)
	   {
        $gn++;
        next if($gn{$gn});
	   }
   next if($_ eq "");

##忽略特殊符号
   next if($opt_x && $_=~/$opt_x/);
   
   my @a=split(/\t/,$_,-1);
##缺失值替换   
   map{$_ = $opt_d if($_ eq "")}@a;

##整合特定的feature列
   $f=join ("\t",@a[@f1]);
###是否过滤特定的列
   if(@v==0 |$v[0] eq "all")
	   {
	   @v=grep{!$gv{$_}}(0..$#a);
	   }
##重新排列
	   @a=($f,@a[@v]);

  $r++;

###如果有抬头
  if($r==1 && $opt_t eq "y")
	   {
		   @h=@a;
		   $h=$a[0];
		   ##储存第一行样本名
           for(1..$#h)
			   {
			   if($opt_n)
				   {$H[$#H+1]="$id\_$h[$_]";}
			   else
				   {$H[$#H+1]="$h[$_]";}
			   $c[$#c+1]=@c+1;
			   }
       next;
	   }
	   


   ###如果没有抬头
   if($r == 1 && $opt_t eq "n")
	    {
	     @h=@a;
		 $h=join("\t",map{"Column_".$_}@f1) if(@f1>0);
		   ##储存第一行样本名
           for(1..$#h)
			   {
			   $c[$#c+1]=@c+1;
			   if($opt_n)
				   {$H[$#H+1]="$id\_ValueColumn\_$v[$_-1]";}
			   else
				   {$H[$#H+1]="ValueColumn\_$v[$_-1]";}
	   
			   }
		}

   next if($a[0] eq "");

	if($#a != $#h)
		 {
		 print "$id\t$a[0]\t单元格数目不对，已忽略!!\n\n";
		 next;
	     }


    ##储存第一列特征名

     push(@F,$a[0]);

     $r[$#r+1]=@r+1;
	##哈希记录行列values
    for(1..$#h)
	    {
        $v{$r[-1]}{$c[-$#h+$_-1]}=$a[$_];
		}

   }

}


if($opt_c eq "y")
{
@H=&check_dup(@H);
}

if($opt_r  eq "y")
{
@F=&check_dup(@F);
}

###如果为数字，重复名加和
open LOG,">$opt_o/merge_log";
for(0..$#r)
{
  my $k1=$_;
  for(0..$#c)
	{
	  my $k2=$_;
	  next unless(exists $v{$r[$k1]}{$c[$k2]});
	  ####相同行列抬头加和处理
	  if($opt_m eq "num")
	  {
		  if(!$V{$F[$k1]}{$H[$k2]})
			{
			$V{$F[$k1]}{$H[$k2]} = $v{$r[$k1]}{$c[$k2]};
			}
			else
			{
            print LOG "$H[$k2] 文件中feature $F[$k1] 存在重复！将对Value $V{$F[$k1]}{$H[$k2]} 和 $v{$r[$k1]}{$c[$k2]} 进行去重复合并！\n";
			
			$V{$F[$k1]}{$H[$k2]} += $v{$r[$k1]}{$c[$k2]};
			}
	  }
	  elsif($opt_m eq "chr")
		{
	    if(!$V{$F[$k1]}{$H[$k2]})
			{
			$V{$F[$k1]}{$H[$k2]}="$v{$r[$k1]}{$c[$k2]}";
			}
			else
			{
		    print LOG "$H[$k2] 文件中feature $F[$k1] 存在重复！将对Value $V{$F[$k1]}{$H[$k2]} 和 $v{$r[$k1]}{$c[$k2]} 进行去重复合并！\n";
			$V{$F[$k1]}{$H[$k2]}.=";$v{$r[$k1]}{$c[$k2]}";
			my @vtmp=split("\;",$V{$F[$k1]}{$H[$k2]});
			@vtmp=uniq(@vtmp);
			$V{$F[$k1]}{$H[$k2]}=join("\;",@vtmp);
			}
		}
	}
}
close LOG;
###删除临时文件
system("rm $tmp_list.tmplist")if($opt_i =~ /$tmp_list.tmplist/);

###输出
@H=uniq(@H);
@F=uniq(@F);
@H=sort @H if(grep($_ =~/^c$/i,@s)>0);
@F=sort @F if(grep($_ =~/^r$/i,@s)>0);
@F=reverse @F if(grep($_ =~/^v$/i,@s)>0);


open O,">$opt_o";
print O "$h";
my @h;
if($opt_u eq "n")
{@h=@H;}
elsif($opt_u eq "c")
{@h=&uniname(@H);}
elsif($opt_u eq "e")
{
	@h=&uniname(@H);
	@h=map{$_=$h[-2].$_;$_}@h;
}
elsif($opt_u eq "f")
{
	@h=&uniname(@H);
	@h=map{$_=$_.$h[-1];$_}@h;
}

for(0..$#H)
{
print O "\t$h[$_]";
}
print O "\n";

for(@F)
{
my $tf=$_;
my @out;
 for(@H)
	{
     if(exists $V{$tf}{$_})
		{
		 push(@out,$V{$tf}{$_});
		 }
	 else
		{
		 push(@out,$opt_d);
		 }
    }


if($opt_m eq "num")
	{
my	@out1=grep{$_==0}@out;
next if($#out1==$#out);
	}
	elsif($opt_m eq "chr")
	{
my	@out2=grep{$_ eq "$opt_d"}@out;	
next if($#out2==$#out);
	}

print O $tf;
for(@out)
	{
    print O "\t$_";
    }


	print O "\n";
}


close O;



sub check_dup
{
my %num;
my @tmp;
for(@_)
	{
	my $t=$_;
	$t=$1 if(/^(.+?)_dup(\d+)$/);
	$num{$t}++;
	my $num=$num{$t}-1;
	if(grep($_ eq $t,@tmp))
	{
	print "$t have dup!\n";
	$t="$t\_dup$num";
	}
	push(@tmp,$t);
    }
    return @tmp;
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



##############################################################
$date=localtime;
print "file merge has been completed! 
Finish Date:$date\n";

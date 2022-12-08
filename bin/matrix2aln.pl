#!/workdir/staff/chenjing/Biosoft/bin/miniconda/bin/perl
use Getopt::Std;
use FindBin qw($Bin $RealBin $Script);
use List::Util qw/max min/;
use List::MoreUtils qw(firstidx);
use Cwd 'abs_path';
use lib "$FindBin::Bin/lib";

sub usage{
	print STDOUT "usage: $0  -i matrix文件路径  -o ./\n";
	print STDOUT "Options\n";
	print STDOUT "\t-i : 输入matrix文件\n";
    print STDOUT "\t-o : 输出文件夹路径，将生成两份输出文件\n";
    print STDOUT "\t-p : 输出文件的前缀名\n";
    print STDOUT "\t-f : 过滤 - 的比率，不能同时和-c使用\n";
    print STDOUT "\t-c : 过滤 - 的数目，设置该参数-f将失效\n";
	print STDOUT "\t-h : show the help message\n";
	exit(1);
}

our($opt_i,$opt_o,$opt_f,$opt_d,$opt_c,$opt_g,$opt_m,$opt_r,$opt_p,$opt_h);
getopts('i:o:d:f:c:r:p:g:m:h');
&usage if($opt_h);
$opt_p ||="SNPs";
$opt_o||="./";
$opt_r ||= "n";
$opt_d ||= "-";
my $out;
open M,"$opt_i";
chomp(my $h=<M>);
my @h=split(/\t/,$h);
$opt_c ||= $opt_f*$#h;
@h=map{$_=$1 if(/^(.+?)\.(ffn|fa|fna|fasta|cds|faa|pro|res)$/);$_=$1 if(/^(.+?)\_vs_/);$_}@h;
push(@out,join("\t",@h));

while(<M>)
{
chomp;
my @a=split(/\t/,$_,-1);
my @ref=split(/\|/,$a[0]);
my @tf=grep{$_ eq $opt_d }@a;
next if(@tf>$opt_c);
push(@out,join("\t",@a));
$s{'ref'}.="$ref[-1]";
for(1..$#h)
	{
    $s{$h[$_]}.="$a[$_]";
    }
}
open O,">$opt_o/$opt_p.matridx_filter";
print O join("\n",@out);
close O;

open ALN,">$opt_o/$opt_p.aln";
for(keys %s)
{
next if($_ eq 'ref' && $opt_r eq "n");
print ALN ">$_\n$s{$_}\n";
}
close ALN;

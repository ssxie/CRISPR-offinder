# The purpose of the CRISPR-offinder is to find CRISPR targets in sequence
# provided by a user that have identified as having no off-targets in
# a specfied genome through an algorithm that searches the genome for similar sequences.

# More information on the CRISPR-offinder can be found in the article

#   Shengsong Xie#, Wubin Qu#, Changzhi Zhao, Guanglei Li, Xiangdong Liu, Qianzhi Shao, 
#   Haiwei Li and Shuhong Zhao (2016) CRISPR-offinder: a CRISPR guide RNA design and 
#   off-target searching tool for user-defined protospacer adjacent motif. Submitting.

# and on the web site www.crispr-offinder.com or http://crispr.igenetech.com. These web sites
# also contains full documentation on how to install and run the system.

# CRISPR-offinder requires an OpenCL-enabled device (CPU, GPU, or etc..) and corresponding 
# runtime driver pre-installed to run properly.
# Please download and install the latest driver below:
# 1.Intel (Download 'OpenCL runtime' in the middle of the page): 
#         http://software.intel.com/en-us/vcsource/tools/opencl-sdk
# 2.NVidia: http://www.nvidia.com/Download/index.aspx
# 3.AMD: http://support.amd.com/en-us/download
# Before installing CRISPR-offinder, please check whether your device is an OpenCL-supported one.

# Supported OS
# 1.Linux (with proprietary drivers installed)
# 2.Max OS X (Snow leopard or higher)
# 3.Windows 7 or higher (XP or below is not supported)

# Whole genome of target organism is needed (in FASTA format). You can find one in one of the below links:
# 1.UCSC genome sequences library
# 2.Ensembl sequence library
# Extract all FASTA files in a directory. Remember the full path of the FASTA files directory.

# COPYRIGHT NOTICE
# 
# GNU GENERAL PUBLIC LICENSE
#
# Copyright (c) 2016 Shengsong Xie and Wubin Qu 
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
# The CRISPR-offinder is being distributed under the GPL license
# Please cite the article listed above if you use, copy, or employ code
# from the CRISPR-offinder.

#!/usr/bin/perl -w
use Getopt::Long;
use File::Basename;
use threads;
use FindBin qw($Bin);

BEGIN{
	push(@INC,$Bin);
};

my $usage="Welcome to CRISPR-offinder
	--a CRISPR guide RNA design and off-target searching tool for 
	user-defined protospacer adjacent motif
Usage:
	perl $0 [Options]
Version:
	1.0
Options:
	-i            [s] Input file <required>
	-o            [s] Output dir <default: ./>
	-pam          [s] protospacer adjacent motif(PAM) <required>
	-pamtype      [s] motif type f/r f:forword,5' r:reverse,3' <deafult: f>
	-length       [i] Length of protospacer <default: 20>
	-gc_min       [i] The minimum value of GC content <default: 20>
	-gc_max       [i] The maximum value of GC content <default: 80>
	-mismatches   [i] Number of mismatches[0-9] <default: 5>
	-strand       [s] Searching CRISPR target sites using DNA strands based option(s/a/b) <default: s>
	-cga          [s] (C: using CPUs, G: using GPUs, A: using accelerators) <default: C>
	-gd           [s] genome dir <default: $Bin/genome>
	-system       [s] run system (Linux32/Linux64/Mac) <default: Linux64>
	-offset_start [i] The minimum value of sgRNA offset <default: -2>
	-offset_end   [i] The maximum value of sgRNA offset <default: 32>

PAM requirement:
    NGG - SpCas9 from Streptococcus pyogenes - direction: 3¡¯
    NRG - SpCas9 from Streptococcus pyogenes - direction: 3¡¯
    NNAGAAW - StCas9 from Streptococcus thermophilus - direction: 3¡¯
    NNNNGMTT - NmCas9 from Neisseria meningitidis - direction: 3¡¯
    NNGRRT - SaCas9 from Staphylococcus aureus - direction: 3¡¯
    NNNRRT - SaCas9 KKH variant - direction: 3¡¯
    NGG(reduced NAG binding) - SpCas9 D1135E variant - direction: 3¡¯
    NGCG - SpCas9 VRER variant - direction: 3¡¯
    NGAG - SpCas9 EQR variant - direction: 3¡¯
    NGAN-NGNG - SpCas9 VQR variant - direction: 3¡¯
    NGG - FnCas9 from Francisella novicida - direction: 3¡¯
    YG - FnCas9 RHA variant - direction: 3¡¯
    TTTN - AsCpf1 from Acidaminococcus, LbCpf1 from Lachnospiraceae - direction: 5¡¯
    TTN - FnCpf1 from Francisella novicida strain U112 - direction: 5¡¯
    CTA - FnCpf1 from Francisella novicida strain U112 - direction: 5¡¯
    TTN-CTA - FnCpf1 from Francisella novicida strain U112 - direction: 5¡¯
    TTN - C2c1 from four major taxa: Bacilli, Verrucomicrobia, a-proteobacteria, and d-proteobacteria - direction: 5¡¯
    Custom, Enter your PAM - enter user defined-PAM

    Note that CRISPR-offinder allows mixed bases to account for the degeneracy in PAM sequences
    Code  Base     Code  Base
    A 	Adenine 	K 	G or T
    C 	Cytosine 	M 	A or C
    G 	Guanine 	B 	C or G or T
    T 	Thymine 	D 	A or G or T
    R 	A or G  	H 	A or C or T
    Y 	C or T  	V 	A or C or G
    S 	G or C  	N 	any base
    W 	A or T 	  	 
Authors:
	Shengsong Xie, Wubin Qu
Email:
	ssxieinfo\@gmail.com, wubin.qu\@igenetech.com 
Homepage:
	http://www.crispr-offinder.com, http://crispr.igenetech.com
Copyright:
	Free for research and educational users.
";

unless(@ARGV>0){
	print "$usage";
	exit;
}

GetOptions(
	"i=s" => \$Inputfile_Fasta,             #Input file
	"o=s" => \$outdir,                      #Output dir
	"pam=s" => \$PAM,                       #protospacer adjacent motif(PAM)
	"pamtype=s" => \$PAMtype,               #motif type f/r f:5' r:3' [f]
	"length=i" => \$truncat,                #Length of protospacer[20]
	"gc_min=i" => \$GC_l,                   #The minimum value of GC content [20]
	"gc_max=i" => \$GC_m,                   #The maximum value of GC content [80]
	"mismatches=i" => \$Mismatches_num,     #Number of mismatches[0-9]
	"strand=s" => \$Option,                 #Searching CRISPR target sites using DNA strands based option(s/a/b)
	"system=s" => \$System,                 #run system (Linux32/Linux64/Mac) <default: Linux64>
	"cga=s" => \$CGA,                       #(C: using CPUs, G: using GPUs, A: using accelerators)
	"gd=s" => \$genomedir,                  #genome dir
	"offset_start=i" => \$offset_s,         #The minimum value of sgRNA offset [-2]
	"offset_end=i" => \$offset_e,           #The maximum value of sgRNA offset [32]
);

# default #
$truncat ||= "20";  
$GC_l ||= "20";
$GC_m ||= "80";
$CGA ||= "C";
$offset_s ||="-2";
$offset_e ||="32";
$outdir||="./";
$Option||="s";
$genomedir||="$Bin/genome";
$System||="Linux64";
$Mismatches_num||=5;
$PAMtype||="f";
unless($System eq "Linux32" || $System eq "Linux64" || $System eq "Mac"){print "Error: system type $System is not support!\n";print "$usage";exit;}
# code #
my %par = (
	'K' => "[GT]",
	'M' => "[AC]",
	'R' => "[AG]",
	'S' => "[GC]",
	'W' => "[AT]",
	'Y' => "[CT]",
	'B' => "[CGT]",
	'D' => "[AGT]",
	'H' => "[ACT]",
	'V' => "[ACG]",
	'N' => "[ATCG]",
);
my %sgRNA = (
	'NRG' => "NGG",
);
my %offtarget = (
	'TTN-CTA' => "YTN",
);
# create #
unless(-e "$outdir/CRISPR_offinder.report_$truncat.$Option.$PAM/input_file"){system "mkdir -p $outdir/CRISPR_offinder.report_$truncat.$Option.$PAM/input_file"};
unless(-e "$outdir/CRISPR_offinder.report_$truncat.$Option.$PAM/off-target"){system "mkdir -p $outdir/CRISPR_offinder.report_$truncat.$Option.$PAM/off-target"};
my $local_time = localtime;
print "$local_time begin ...\n";
#system "dos2unix $Inputfile_Fasta";
# step1 #
open IN,"$Inputfile_Fasta" or die $!;
open FASTA,">$outdir/CRISPR_offinder.report_$truncat.$Option.$PAM/TargetSeq.fa" or die $!;
my ($id,$seq_s,$seq_a);
while(<IN>){
	$_=~s/(\r\n|\n|\r)/\n/g;
	chomp;
	if(m/^\>/){$id=$_;next;}
	$seq_s = $_;
	if($Option eq "b" || $Option eq "s"){
		print FASTA "${id}_S\n$seq_s\n";
	}
	if($Option eq "b" || $Option eq "a"){
		$seq_s =~ tr/atucgACGUT/TAAGCTGCAA/;
		$seq_a = reverse($seq_s); 
		print FASTA "${id}_A\n$seq_a\n";
	}
}
close IN;
close FASTA;
# step2 single #
open IN,"$outdir/CRISPR_offinder.report_$truncat.$Option.$PAM/TargetSeq.fa" or die $!;
open OUT,">$outdir/CRISPR_offinder.report_$truncat.$Option.$PAM/report_protospacer_single.txt" or die $!;
open OUTFA,">$outdir/CRISPR_offinder.report_$truncat.$Option.$PAM/CRISPR.targets_single.fa" or die $!;
open OUTS,">$outdir/CRISPR_offinder.report_$truncat.$Option.$PAM/CRISPR.targets_S.txt" or die $!;
open OUTA,">$outdir/CRISPR_offinder.report_$truncat.$Option.$PAM/CRISPR.targets_A.txt" or die $!;
print OUT "sgRID\t"."Start\t"."End\t"."CRISPR_target_sequence(5'-3')\t"."Length(nt)\t"."GC%\n";
$id="";
my @filelist;
while(<IN>){
	chomp;
	if(m/^\>/){$id=$_;$id=~s/^\>//;next;}
	$_ = uc $_;
	$_=~tr/U/T/;
	my ($num,$match);
	if(exists $sgRNA{$PAM}){$match = ana_pam($sgRNA{$PAM});}
	else{$match = ana_pam($PAM);}
	if($PAMtype eq "f"){
		while($_=~/($match)[ATCG]{$truncat}/mg){
			my $sgR = $&;
			my $before = $`;
			my $after = $';
			pos($_)= (length $before) + 1;
			my $sgRGC = ana_gc($sgR);
			unless($sgRGC>=$GC_l && $sgRGC<=$GC_m){next;}
			my $offpam = ana_off($PAM);
			my $motif1 = join("","$offpam","N" x $truncat);
			my $motif2 = join("","N" x length $1,substr($sgR,length $1));$num++;
			my ($sgRidx_s,$sgRidx_e,$sgRidx_A_s,$sgRidx_A_e) = ana_sel($before,$after,$sgR,$_);
			if($id=~/._S$/){
				print OUTS "$id\_$num\t"."$sgRidx_s\t"."$sgRidx_e\t"."$sgR\t".(length $sgR)."\t$sgRGC%\n";
				print OUT "$id\_$num\t"."$sgRidx_s\t"."$sgRidx_e\t"."$sgR\t".(length $sgR)."\t$sgRGC%\n";
			}
			elsif($id=~ /._A$/){
				print OUTA "$id\_$num\t"."$sgRidx_A_s\t"."$sgRidx_A_e\t"."$sgR\t".(length $sgR)."\t$sgRGC%\n";
				print OUT "$id\_$num\t"."$sgRidx_A_s\t"."$sgRidx_A_e\t"."$sgR\t".(length $sgR)."\t$sgRGC%\n";
			}
			print OUTFA ">$id\_$num\n"."$sgR\n";
			open TEMP,">$outdir/CRISPR_offinder.report_$truncat.$Option.$PAM/input_file/$id\_$num.txt";
			print TEMP "$genomedir\n"."$motif1\n"."$motif2\t$Mismatches_num\n";
			push(@filelist,"$outdir/CRISPR_offinder.report_$truncat.$Option.$PAM/input_file/$id\_$num.txt");
			close TEMP;
		}
	}elsif($PAMtype eq "r"){
		while($_=~/[ATCG]{$truncat}($match)/mg){
			my $sgR = $&;
			my $before = $`;
			my $after = $';
			pos($_)= (length $before) + 1;
			my $sgRGC = ana_gc($sgR);
			unless($sgRGC>=$GC_l && $sgRGC<=$GC_m){next;}
			my $offpam = ana_off($PAM);;
			my $motif1 = join("","N" x $truncat,"$offpam");
			my $motif2 = join("",substr($sgR,0,$truncat),"N" x length $1);$num++;
			my ($sgRidx_s,$sgRidx_e,$sgRidx_A_s,$sgRidx_A_e) = ana_sel($before,$after,$sgR,$_);
			if($id=~/._S$/){
				print OUTS "$id\_$num\t"."$sgRidx_s\t"."$sgRidx_e\t"."$sgR\t".(length $sgR)."\t$sgRGC%\n";
				print OUT "$id\_$num\t"."$sgRidx_s\t"."$sgRidx_e\t"."$sgR\t".(length $sgR)."\t$sgRGC%\n";
			}
			elsif($id=~ /._A$/){
				print OUTA "$id\_$num\t"."$sgRidx_A_s\t"."$sgRidx_A_e\t"."$sgR\t".(length $sgR)."\t$sgRGC%\n";
				print OUT "$id\_$num\t"."$sgRidx_A_s\t"."$sgRidx_A_e\t"."$sgR\t".(length $sgR)."\t$sgRGC%\n";
			}
			print OUTFA ">$id\_$num\n"."$sgR\n";
			open TEMP,">$outdir/CRISPR_offinder.report_$truncat.$Option.$PAM/input_file/$id\_$num.txt";
			print TEMP "$genomedir\n"."$motif1\n"."$motif2\t$Mismatches_num\n";
			push(@filelist,"$outdir/CRISPR_offinder.report_$truncat.$Option.$PAM/input_file/$id\_$num.txt");
			close TEMP;
		}
	}
}
close IN;
close OUT;
close OUTFA;
close OUTS;
close OUTA;
# step2 pair #
open PA,"$outdir/CRISPR_offinder.report_$truncat.$Option.$PAM/CRISPR.targets_A.txt" or die $!;
open PS,"$outdir/CRISPR_offinder.report_$truncat.$Option.$PAM/CRISPR.targets_S.txt" or die $!;
open PAIRS, ">$outdir/CRISPR_offinder.report_$truncat.$Option.$PAM/report_protospacer_pairs.xls" or die $!;
print PAIRS "sgRID_S\ttarget_seq_S\tStart_S\tEnd_S\tGC%_S\t<->\tsgRID_A\ttarget_seq_A\tStart_A\tEnd_A\tGC%_A\tsgRNA_offset(bp)\n";
my %hash;
while(<PS>){
	chomp;
	$hash{$_}=1;
}
while(<PA>){
	chomp;
	foreach my $ps (keys %hash){
		my ($sgRID_A,$Start_A,$End_A,$target_seq_A,$Pattern_A,$GC_A)=(split /\t/,$_)[0,1,2,3,4,5];
		my ($sgRID_S,$Start_S,$End_S,$target_seq_S,$Pattern_S,$GC_S)=(split /\t/,$ps)[0,1,2,3,4,5];
		$ID_A = $sgRID_A;$ID_A=~ s/_A_(\d+)//m;
		$ID_S = $sgRID_S;$ID_S=~ s/_S_(\d+)//m;
		my $offset_value = $Start_S -$End_A;
		if (($ID_A eq $ID_S) and ($offset_value >= $offset_s and $offset_value <= $offset_e)) {      # -2 to 32 bp or 5 to 35 bp
			print PAIRS "$sgRID_A"."\t"."$target_seq_A"."\t"."$Start_A"."\t"."$End_A"."\t"."$GC_A"."\t<->\t";
			print PAIRS "$sgRID_S"."\t"."$target_seq_S"."\t"."$Start_S"."\t"."$End_S"."\t"."$GC_S"."\t"."$offset_value"."\n";
		}
	}
}
close PA;
close PS;
close PAIRS;
# step3 cas-offinder #
my @fileresult;my @thread;
my $sth=0;my $eth=0;
for(my $i=0;$i<@filelist;$i++){
	my $sample = basename $filelist[$i];
	push(@fileresult,"$outdir/CRISPR_offinder.report_$truncat.$Option.$PAM/off-target/$sample");
	$thread[$i] = threads->new(\&offinder,$filelist[$i],$sample);
	$eth=$i;
	if($i%5==0){
		for(my $m=$sth;$m<=$eth;$m++){$thread[$m]->join;}
		$sth=$eth+1;
	}
}
for(my $m=$sth;$m<=$eth;$m++){$thread[$m]->join;}
# count of target
open OUT,">$outdir/CRISPR_offinder.report_$truncat.$Option.$PAM/count_offtarget.txt" or die $!;
#print OUT "#ID\t0M\t1M\t2M\t3M\t4M\t5M\t6M\t7M\t8M\t9M\tTotal\n";
print OUT "0M\t1M\t2M\t3M\t4M\t5M\t6M\t7M\t8M\t9M\tTotal\n";
foreach my $file (@fileresult){
	open IN,"$file" or die $!;
	open TEMP,">$file.temp" or die $!;
	my $name = basename $file;
	my $total=0;
	my %hashmis;
	$name =~ s/.txt$//;
	my $offmatch = ana_pam($PAM);
	while(<IN>){
		chomp;
		my ($crRNA, $Chr, $Position, $DNA, $Direction, $Mismatches2) =split(/\s+/, $_);
		if($PAMtype eq "f"){
			unless($DNA=~/^$offmatch/mg){print "message: drop improper result $DNA\n";next;}
		}
		elsif($PAMtype eq "r"){
			unless($DNA=~/$offmatch$/mg){print "message: drop improper result $DNA\n";next;}
		}
		$hashmis{$Mismatches2}++;
		$total++;
		print TEMP "$_\n";
	}
	close IN;
	close TEMP;
	system "mv $file.temp $file";
	#print OUT "$name";
	for(my $i=0;$i<10;$i++){
		if(exists $hashmis{$i}){print OUT "$hashmis{$i}\t";}
		else{print OUT "0\t";}
	}
	if(exists $hashmis{0}){
		$total--;
	}
	print OUT "$total\n";
}
close OUT;
system "paste $outdir/CRISPR_offinder.report_$truncat.$Option.$PAM/report_protospacer_single.txt $outdir/CRISPR_offinder.report_$truncat.$Option.$PAM/count_offtarget.txt >$outdir/CRISPR_offinder.report_$truncat.$Option.$PAM/CRISPR_offinder_report.xls";
$local_time = localtime;
print "$local_time end ...\n";

#delete Temporary files.
unlink ("$outdir/CRISPR_offinder.report_$truncat.$Option.$PAM/count_offtarget.txt")||die "Can't delete count_offtarget.txt file";
unlink ("$outdir/CRISPR_offinder.report_$truncat.$Option.$PAM/CRISPR.targets_A.txt")||die "Can't delete CRISPR.targets_A.txt file";
unlink ("$outdir/CRISPR_offinder.report_$truncat.$Option.$PAM/CRISPR.targets_S.txt")||die "Can't delete CRISPR.targets_S.txt";
unlink ("$outdir/CRISPR_offinder.report_$truncat.$Option.$PAM/CRISPR.targets_single.fa")||die "Can't delete CRISPR.targets_single.fa file";
unlink ("$outdir/CRISPR_offinder.report_$truncat.$Option.$PAM/report_protospacer_single.txt")||die "Can't delete report_protospacer_single.txt file";

# hanshu #
sub ana_pam
{
	my $p = shift;
	my @a = split(//,$p);
	foreach my $v (@a){
		if(exists $par{$v}){$v=$par{$v};}
		if($v eq "-"){$v="|";}
	}
	return join("",@a);
}

sub ana_gc
{
	my $seq = shift;
	my $len = length $seq;
	my @all = split (//,$seq);
	my $gc;
	foreach my $a (@all){
		if($a=~/[GCgc]/){$gc++;}
	}
	my $gp = sprintf("%.2f",$gc*100/$len);
	return $gp;
}

sub ana_sel
{
	my ($s,$e,$l,$all) = @_;
	my $sgRidx_s = (length $s) + 1;
	my $sgRidx_e = $sgRidx_s + (length $l) - 1;
	my $sgRidx_A_s = (length $all)-$sgRidx_e+1;
	my $sgRidx_A_e = (length $all)-$sgRidx_s+1;
	return ($sgRidx_s,$sgRidx_e,$sgRidx_A_s,$sgRidx_A_e);
}

sub ana_real
{
	my ($p,$m,$v) = @_;
	unless($p=~/-/){return $p;}
	my @a = split (/-/,$p);
	my @b = split (/\|/,$m);
	for(my $i=0;$i<@b;$i++){
		if($v=~/$b[$i]/){return $a[$i];}
	}
}

sub ana_off
{
	my $p = shift;
	if(exists $offtarget{$p}){return $offtarget{$p};}
	else{return $p;}
}

sub offinder
{
	my ($file,$sample) = @_;
	system "$Bin/bin/cas-offinder_${System} $file $CGA $outdir/CRISPR_offinder.report_$truncat.$Option.$PAM/off-target/$sample";
}
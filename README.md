ABOUT
CRISPR-offinder 
- A program used for designing CRISPR guide RNA and
searching off-targets for user-defined protospacer adjacent motif(PAM)

CRISPR-offinder requires an OpenCL-enabled device (CPU, GPU, or etc..) and corresponding 
runtime driver pre-installed to run properly.

Please download and install the latest driver below:
1.Intel (Download 'OpenCL runtime' in the middle of the page): 
http://software.intel.com/en-us/vcsource/tools/opencl-sdk
2.NVidia: http://www.nvidia.com/Download/index.aspx
3.AMD: http://support.amd.com/en-us/download

Before installing CRISPR-offinder, please check whether your device is an OpenCL-supported one and make sure that Cas-OFFinder program work well. Bae, S. et al. (2014) Cas-OFFinder: a fast and versatile algorithm that searches for potential 
off-target sites of Cas9 RNA-guided endonucleases. Bioinformatics., 30, 1473-1475.

Supported OS
1.Linux (with proprietary drivers installed)
2.Max OS X (Snow leopard or higher)
3.Windows 7 or higher (XP or below is not supported)

VERSION
         Version: 1.0
         Feb 24, 2016

SYNOPSIS
         $ perl CRISPR-offinder.pl <option>
         For help information, type 'perl CRISPR-offinder.pl'

DESCRIPTION
         This script takes a FASTA file of the target sites and queries the
         reference database by CRISPR-offinder tolerant designated mismatches, and gives
         the human readable and computer mineable target/off target sites, 
	 also the summary file.

AUTHOR
         Shengsong Xie and Wubin Qu
         E-mail: ssxie@mail.hzau.edu.cn or wubin.qu\@igenetech.com

Homepage
         http://www.crispr-offinder.com, http://crispr.igenetech.com

LICENSE
         Free for research and educational users.
		 
CRISPR-offinder.pl version 1.0
==============================

Usage:

  perl CRISPR-offinder.pl <option>
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

Whole genome of target organism is needed (in FASTA format). You can find one in one of the below links:
1.UCSC genome sequences library
2.Ensembl sequence library
Extract all FASTA files in a directory. Remember the full path of the FASTA files directory(for option -gd).

For help information: perl CRISPR-offinder.pl

PAM requirement:
    NGG - SpCas9 from Streptococcus pyogenes - direction: 3’
    NRG - SpCas9 from Streptococcus pyogenes - direction: 3’
    NNAGAAW - StCas9 from Streptococcus thermophilus - direction: 3’
    NNNNGMTT - NmCas9 from Neisseria meningitidis - direction: 3’
    NNGRRT - SaCas9 from Staphylococcus aureus - direction: 3’
    NNNRRT - SaCas9 KKH variant - direction: 3’
    NGG(reduced NAG binding) - SpCas9 D1135E variant - direction: 3’
    NGCG - SpCas9 VRER variant - direction: 3’
    NGAG - SpCas9 EQR variant - direction: 3’
    NGAN-NGNG - SpCas9 VQR variant - direction: 3’
    NGG - FnCas9 from Francisella novicida - direction: 3’
    YG - FnCas9 RHA variant - direction: 3’
    TTTN - AsCpf1 from Acidaminococcus, LbCpf1 from Lachnospiraceae - direction: 5’
    TTN - FnCpf1 from Francisella novicida strain U112 - direction: 5’
    CTA - FnCpf1 from Francisella novicida strain U112 - direction: 5’
    TTN-CTA - FnCpf1 from Francisella novicida strain U112 - direction: 5’
    TTN - C2c1 from four major taxa: Bacilli, Verrucomicrobia, a-proteobacteria, and d-proteobacteria - direction: 5’
    Custom, Enter your PAM - enter user defined-PAM

    Code  Base     Code  Base
    A 	Adenine 	K 	G or T
    C 	Cytosine 	M 	A or C
    G 	Guanine 	B 	C or G or T
    T 	Thymine 	D 	A or G or T
    R 	A or G  	H 	A or C or T
    Y 	C or T  	V 	A or C or G
    S 	G or C  	N 	any base
    W 	A or T 	 
    
    Note that CRISPR-offinder allows mixed bases to account for the degeneracy in PAM sequences

Install of CRISPR-offinder:

unzip CRISPR-offinder_1.0.zip

Change the access permissions to Cas-OFFinder in Linux system, type:
    cd bin
    chmod u+x cas-offinder_Linux32 cas-offinder_Linux64 cas-offinder_Mac

For citation:

    Shengsong Xie#, Wubin Qu#, Changzhi Zhao, Guanglei Li, Xiangdong Liu, Qianzhi Shao, 
    Haiwei Li and Shuhong Zhao (2016) CRISPR-offinder: a CRISPR guide RNA design and 
    off-target searching tool for user-defined protospacer adjacent motif.  
    Submitting.	  


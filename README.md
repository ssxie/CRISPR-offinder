ABOUT

CRISPR-offinder - A program used for designing CRISPR guide RNA and
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
         E-mail: ssxie\@mail.hzau.edu.cn or wubin.qu\@igenetech.com

Homepage
         http://www.crispr-offinder.com, http://crispr.igenetech.com

LICENSE
         Free for research and educational users.
		 
CRISPR-offinder.pl version 1.0
==============================

Usage:

  perl CRISPR-offinder.pl <option>

Whole genome of target organism is needed (in FASTA format). You can find one in one of the below links:
1.UCSC genome sequences library
2.Ensembl sequence library
Extract all FASTA files in a directory. Remember the full path of the FASTA files directory(for option -gd).

For help information: perl CRISPR-offinder.pl

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


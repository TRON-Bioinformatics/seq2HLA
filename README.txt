##########################################################################################################
#Title:
#seq2HLA - HLA typing from RNA-Seq sequence reads
#
#Release: 2.3
#
#Author:
#Sebastian Boegel, 2012 - 2014 (c)
#TRON - Translational Oncology at the University Medical Center Mainz, 55131 Mainz, Germany
#University Medical Center of the Johannes Gutenberg-University Mainz, III.  Medical Department, Mainz, Germany
#
#Contact:
#boegels@uni-mainz.de
#
#Synopsis:
#We developed an in-silico method "Seq2HLA", written in python and R, which takes standard RNA-Seq sequence reads in fastq format
#as input, uses a bowtie index comprising all HLA alleles and outputs the most likely HLA class I and class II genotypes (in 4 digit resolution),
#a p-value for each call, and the expression of each class
#
#Usage:
#python seq2HLA.py -1 <readfile1> -2 <readfile2> -r "<runname>" [-p <int>]* [-3 <int>]**
#*optional: number of parallel search threads for bowtie optional (Default:6)
#**optional: trim int bases from the low-quality end of each read
#readfile can be uncompressed or gzipped fastq file
#runname should contain path information, e.g. "folder/subfolder/..../run", in order to store all resulting into to folder and all filenames will have the suffix run-

#Output:
#The results are output to stdout and to textfiles. Most important are:
#i) <prefix>-ClassI.HLAgenotype2digits => 2 digit result of Class I
#ii) <prefix>-ClassII.HLAgenotype2digits => 2 digit result of Class II
#iii) <prefix>-ClassI.HLAgenotype4digits => 4 digit result of Class I
#iv) <prefix>-ClassII.HLAgenotype4digits => 4 digit result of Class II
#v) <prefix>.ambiguity => reports typing ambuigities (more than one solution for an allele possible)
#vi) <prefix>-ClassI.expression => expression of Class I alleles
#vii) <prefix>-ClassII.expression => expression of Class II alleles

#Dependencies:
#0.) seq2HLA is a python script, developed with Python 2.6.8
#1.) bowtie must be reachable by the command "bowtie". seq2HLA was developed and tested with bowtie version 0.12.7 (64-bit). The call to bowtie is invoked with 6 CPUs. You can change that by paramter -p.
#2.) R must be installed, seq2HLA.py was developed and tested with R version 2.12.2 (2011-02-25)
#3.) Input must be paired-end reads in fastq-format
#4.) Index files must be located in the folder "references".
#5.) Packages: biopython (developed with V1.58), numpy (1.3.0)

#Version history:
#2.3: typing of HLA II loci DRA, DPA1 and DPB1, typing of non-classical HLA I alleles (e.g. HLA-G...), cleaning up after execution (deletion of intermediate files) (August 2017)
#2.2: improved performance, automatic detection of read length (option -l no longer required), user can choose number of parralel search threads (-p), seq2HLA now works with automatic path recognition, so it can be invoked from every path (April 2014)
#2.1: supports gzipped fastq files as input
#2.0: 4-digit typing
#1.0: 2-digit typing

#References:
#Boegel, Sebastian; Loewer, Martin; Schaefer, Michael; Bukur, Thomas; Graaf, Jos de; Boisguerin, Valesca et al. (2013): HLA typing from RNA-Seq sequence reads. In: Genome Med 4 (12), S. 102. DOI: 10.1186/gm403.
#Sebastian Boegel, Jelle Scholtalbers, Martin Löwer, Ugur Sahin, John C Castle (2015): In Silico HLA Typing Using Standard RNA-Seq Sequence Reads. In: Molecular Typing of Blood Cell Antigens
#Sebastian Boegel, Martin Löwer, Thomas Bukur, Ugur Sahin, John C Castle (2014): A catalog of HLA type, HLA expression, and neo-epitope candidates in human cancer cell lines. In: OncoImmunology
#Jelle Scholtalbers, Sebastian Boegel, Thomas Bukur, Marius Byl, Sebastian Goerges, Patrick Sorn, Martin Loewer, Ugur Sahin, John C Castle (2015): TCLP: an online cancer cell line catalogue integrating HLA type, predicted neo-epitopes, virus and gene expression. In: Genome Med

#License:
#The MIT License (MIT)
#Copyright (c) 2012 Sebastian Boegel
#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:
#
#The above copyright notice and this permission notice shall be included in
#all copies or substantial portions of the Software.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
#THE SOFTWARE.

###########################################################################################################

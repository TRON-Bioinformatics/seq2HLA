# seq2HLA - HLA typing from RNA-Seq sequence reads

We developed an in-silico method "Seq2HLA", written in python and R, which takes standard RNA-Seq sequence reads in fastq format
as input, uses a bowtie index comprising all HLA alleles and outputs the most likely HLA class I and class II genotypes (in 4 digit resolution), a p-value for each call, and the expression of each class.


## Dependencies

- Python (3.7.16)
  - biopython (>= 1.79)
  - numpy (>= 1.21.5)
  - scipy (>= 1.7.3)
- bowtie (>= 1.3.1)


## Installation

```
git clone git@github.com:TRON-Bioinformatics/seq2HLA.git
conda env create -f environment.yml
conda activate seq2HLA
```


## Usage

```
usage: seq2HLA.py [-h] -1 FQ1 -2 FQ2 [-r RUN_NAME] [-p THREADS] [-3 TRIM3]
                  [-V]

Predict HLA type of paired-end FASTQs

optional arguments:
  -h, --help            show this help message and exit
  -1 FQ1, --fq1 FQ1     Specify first FASTQ file. Gzipped input supported.
  -2 FQ2, --fq2 FQ2     Specify second FASTQ file. Gzipped input supported.
  -r RUN_NAME, --run_name RUN_NAME
                        Specify run name. This will contain path information
                        as well as the output file prefix
  -p THREADS, --threads THREADS
                        Specify number of threads to be used by bowtie
  -3 TRIM3, --trim3 TRIM3
                        Specify trimming cutoff for the low quality end.
                        Bowtie option: -3 <int> trims <int> bases from the low
                        quality 3' end of each read. Default: 0
  -V, --version         show the version number and exit
```


## Output

The results are output to stdout and to textfiles. Most important are:
- `hla_1_*.HLAgenotype2digits` - 2 digit result of Class I
- `hla_2_classical.HLAgenotype2digits` - 2 digit result of Class II
- `hla_1_*.HLAgenotype4digits` - 4 digit result of Class I
- `hla_2_classical.HLAgenotype4digits` - 4 digit result of Class II
- `hla.ambiguity` - reports typing ambuigities (more than one solution for an allele possible)
- `hla_1_*.expression` - expression of Class I alleles
- `hla_2_classical.expression` - expression of Class II alleles


## Version history

- 2.4: Python 3 supported, updated dependencies and installation, cleaned code
- 2.3: typing of HLA II loci DRA, DPA1 and DPB1, typing of non-classical HLA I alleles (e.g. HLA-G...), cleaning up after execution (deletion of intermediate files) (August 2017)
- 2.2: improved performance, automatic detection of read length (option -l no longer required), user can choose number of parralel search threads (-p), seq2HLA now works with automatic path recognition, so it can be invoked from every path (April 2014)
- 2.1: supports gzipped fastq files as input
- 2.0: 4-digit typing
- 1.0: 2-digit typing

## References

- Boegel, Sebastian; Loewer, Martin; Schaefer, Michael; Bukur, Thomas; Graaf, Jos de; Boisguerin, Valesca et al. (2013): HLA typing from RNA-Seq sequence reads. In: Genome Med 4 (12), S. 102. DOI: 10.1186/gm403.
- Sebastian Boegel, Jelle Scholtalbers, Martin Löwer, Ugur Sahin, John C Castle (2015): In Silico HLA Typing Using Standard RNA-Seq Sequence Reads. In: Molecular Typing of Blood Cell Antigens
- Sebastian Boegel, Martin Löwer, Thomas Bukur, Ugur Sahin, John C Castle (2014): A catalog of HLA type, HLA expression, and neo-epitope candidates in human cancer cell lines. In: OncoImmunology
- Jelle Scholtalbers, Sebastian Boegel, Thomas Bukur, Marius Byl, Sebastian Goerges, Patrick Sorn, Martin Loewer, Ugur Sahin, John C Castle (2015): TCLP: an online cancer cell line catalogue integrating HLA type, predicted neo-epitopes, virus and gene expression. In: Genome Med

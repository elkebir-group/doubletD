# doubletD

![Overview of doubletD](doubletD_overview.png)
(a) The first step of most single-cell sequencing technologies involves cell capture where the goal is to encapsulate single cells into droplets, known as *singlets*.
However, errors in this process can lead to three kind of doublets -- *neotypic* doublets, *nested* doublets and *selflets*.
(b) The cells in each isolated droplet *i* undergo whole-genome amplification and sequencing independently.
These processes introduce errors such as allelic dropouts and imbalance in amplification.
(c) The resulting aligned reads are used for variant calling yielding alternate v_{i,j} and total c_{i,j} read counts at each locus of interest *j*.
(d) doubletD uses the observed variant allele frequencies v_{i,j}/c_{i,j} as the key signal, while accounting for sequencing and amplification errors to detect doublets in the sample.

## Contents

  1. [Pre-requisites](#pre-requisites)
  2. [Usage instcructions](#usage)
     * [I/O formats](#io)
     * [DoubletD](#doubletD)
     * [simulation](#simulation)

<a name="pre-requisites"></a>
## Pre-requisites
+ python3 (>=3.6)
+ [numpy](https://numpy.org/doc/)
+ [pandas](https://pandas.pydata.org/pandas-docs/stable/index.html)
<!-- + (optional for simulation pipeline) [snakemake (>=5.2.0)](https://snakemake.readthedocs.io) --?

<a name="usage"></a>
## Usage instructions

<a name="io"></a>
### I/O formats
The input for Jumper is a text based with two input comma-separated dataframes -- one containing the total read counts and another containing the alternate read counts.
For both the files, each row is a different droplet and each column is a loci.
See `data/sample_DP.tsv` and `data/sample_AD.tsv` for an example for both files.
The output is also a datafram with each row for a different droplet and columns, from left to right, posterior probability that the dorplet is a singlet, posterior probability that the droplet is a doublet and prediction for the droplet to be either 'singlet' or 'doublet'.
See `data/sample_prediction.tsv` for an example.

### Arguments

Parameters with default value `None` are estimated from data

    usage: doubletD.py [-h] [--inputTotal INPUTTOTAL]
                       [--inputAlternate INPUTALTERNATE] [--delta DELTA]
                       [--beta BETA] [--mu_hetero MU_HETERO] [--mu_homo MU_HOMO]
                       [--alpha_fp ALPHA_FP] [--alpha_fn ALPHA_FN] [-o OUTPUTFILE]
                       [--noverbose] [--binomial] [--prec PREC] [--missing]

    optional arguments:
      -h, --help            show this help message and exit
      --inputTotal INPUTTOTAL
                            csv file with a table of total read counts for each
                            position in each cell
      --inputAlternate INPUTALTERNATE
                            csv file with a table of alternate read counts for
                            each position in each cell
      --delta DELTA         doublet rate [0.1]
      --beta BETA           Allelic dropout (ADO) rate [0.05]
      --mu_hetero MU_HETERO
                            heterozygous mutation rate [None]
      --mu_homo MU_HOMO     homozygous mutation rate [None]
      --alpha_fp ALPHA_FP   copy false positive error rate [None]
      --alpha_fn ALPHA_FN   copy flase negative error rate [None]
      -o OUTPUTFILE, --outputfile OUTPUTFILE
                            output file name
      --noverbose           do not output statements from internal solvers
                            [default is false]
      --binomial            use cellcoal doublet model [default is false]
      --prec PREC           Precision for Beta Distribution [None]
      --missing             use missing data in the model? [No]



### Example

Here we will show an example of how to run `doubletD`.
The input files are located in the `example` directory.
We run `doubletD` with a prior doublet probabiltity of 0.2 and ADO rate of 0.5 without using missing data in our model.

    $ python scripts/doubletD.py --inputAlternate example/AD.csv --inputTotal example/DP.csv --delta 0.2 --beta 0.05 -o example/prediction.tsv 

This command generates output file `prediction.tsv` in directory `example`.

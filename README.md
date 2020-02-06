[![Release](https://img.shields.io/github/v/release/bcgsc/RNA-Bloom?include_prereleases)](https://github.com/bcgsc/RNA-Bloom/releases)
[![Downloads](https://img.shields.io/github/downloads/bcgsc/RNA-Bloom/total?logo=github)](https://github.com/bcgsc/RNA-Bloom/releases)
[![Conda](https://img.shields.io/conda/dn/bioconda/rnabloom?label=Conda)](https://anaconda.org/bioconda/rnabloom)

<p align="center">
  <img src="rnabloom_logo.png" alt="RNA-Bloom's logo"/>
</p>

**RNA-Bloom** is a fast and memory-efficient *de novo* transcript sequence assembler. It is designed for the following sequencing data types:
* paired-end bulk RNA-seq (strand-specific/agnostic)
* paired-end single-cell RNA-seq (strand-specific/agnostic)
* nanopore RNA-seq (PCR cDNA/direct cDNA/direct RNA)

Written by [Ka Ming Nip](mailto:kmnip@bcgsc.ca) :email:

:copyright: 2018-2020 Canada's Michael Smith Genome Sciences Centre, BC Cancer

--------------------------------------------------------------------------------

## Dependency :pushpin:

* [Java SE Runtime Environment (JRE) 8](http://www.oracle.com/technetwork/java/javase/downloads/jre8-downloads-2133155.html)

* External software used:

| software                                            | short reads            | long reads             |
| --------------------------------------------------- | ---------------------- | ---------------------- |
| [minimap2](https://github.com/lh3/minimap2)         | required               | required               |
| [Racon](https://github.com/lbcb-sci/racon)          | not used               | required               |
| [ntCard](https://github.com/bcgsc/ntCard) >=1.2.1   | required for `-ntcard` | required for `-ntcard` |

:warning: Their executables must be accessible from your `PATH`!

## Installation :wrench:

You may either install RNA-Bloom with `conda` run:
```
conda install -c bioconda rnabloom
```
OR download from GitHub:
 
1. Download the binary tarball `rnabloom_vX.X.X.tar.gz` from the [releases](https://github.com/bcgsc/RNA-Bloom/releases) section.
2. Extract the downloaded tarball with the command:
```
tar -zxf rnabloom_vX.X.X.tar.gz
```
3. RNA-Bloom is ready to use, ie. `java -jar /path/to/RNA-Bloom.jar ...`



## Quick Start for Short Reads :running:

:warning: Input reads must be in either FASTQ or FASTA format and may be compressed with GZIP. 

### (A) assemble bulk RNA-seq data:
```
java -jar RNA-Bloom.jar -left LEFT.fastq -right RIGHT.fastq -revcomp-right -ntcard -t THREADS -outdir OUTDIR
```

### (B) assemble single-cell RNA-seq data:
```
java -jar RNA-Bloom.jar -pool READSLIST.txt -revcomp-right -ntcard -t THREADS -outdir OUTDIR
```
RNA-Bloom was tested on Smart-seq2 and SMARTer datasets.

#### file format for the `-pool` option:

This text file is expected to have 3 columns, ie.

| column 1 | column 2           | column 3            |
| -------- | ------------------ | ------------------- |
| cell ID  | path of left reads | path of right reads |

Columns are separated by space/tab characters, eg.
```
cell1 /path/to/cell1/left.fastq /path/to/cell1/right.fastq
cell2 /path/to/cell2/left.fastq /path/to/cell2/right.fastq
cell3 /path/to/cell3/left.fastq /path/to/cell3/right.fastq
```

### (C) strand-specific assembly:
```
java -jar RNA-Bloom.jar -stranded ...
```
The `-stranded` option indicates that input reads are strand-specific.

Strand-specific reads are typically in the F2R1 orientation, where `/2` denotes *left* reads in *forward* orientation and `/1` denotes *right* reads in *reverse* orientation.

Configure the read file paths accordingly for bulk RNA-seq data:

`-left /path/to/reads_2.fastq -right /path/to/reads_1.fastq`

and for scRNA-seq data:
```
cell1 /path/to/cell1/reads_2.fastq /path/to/cell1/reads_1.fastq
```

### (D) reference-guided assembly:
```
java -jar RNA-Bloom.jar -ref TRANSCRIPTS.fasta ...
```
The `-ref` option specifies the reference transcriptome FASTA file for guiding short-read assembly.



## Quick Start for Nanopore Reads :running:

### (A) assemble nanopore PCR cDNA sequencing data:
```
java -jar RNA-Bloom.jar -long READS.fasta -ntcard -t THREADS -outdir OUTDIR
```
Input reads are expected to be in a mix of both forward and reverse orientations.

### (B) assemble nanopore direct cDNA sequencing data:
```
java -jar RNA-Bloom.jar -long READS.fasta -stranded -revcomp-long -ntcard -t THREADS -outdir OUTDIR
```
Input reads are expected to be only in the reverse orientation.

### (C) assemble nanopore direct RNA sequencing data:
```
java -jar RNA-Bloom.jar -long READS.fasta -stranded -ntcard -t THREADS -outdir OUTDIR
```
Input reads are expected to be only in the forward orientation.

By default, uracil (`U`) is written as `T`. Use the `-uracil` option to write `U` instead of `T` in the output assembly.

ntCard v1.2.1 supports uracil in reads.


## General Settings :gear:

### (A) set Bloom filter sizes automatically:
```
java -jar RNA-Bloom.jar -fpr 0.01 -nk 28077715 ...
```
This sets the size of Bloom filters automatically to accommodate 28,077,715 unique k-mers for a max false positive rate (FPR) of 1%.

Instead of specifying the exact number of k-mers, you may use ntCard to count k-mers:
```
java -jar RNA-Bloom.jar -fpr 0.01 -ntcard ...
```
To use the `-ntcard` option, `ntcard` must be found in your `PATH`.

As a rule of thumb, a lower FPR may result in a better assembly but requires more memory for a larger Bloom filter.

### (B) set the total size of Bloom filters:
```
java -jar RNA-Bloom.jar -mem 10 ...
```
This sets the total size to 10 GB. If neither `-nk`, `-ntcard`, or `-mem` are used, then the total size is configured based on the size of input read files.

### (C) stop at an intermediate stage:
```
java -jar RNA-Bloom.jar -stage N ...
```
| N   | short reads          | long reads           |
| --- | -------------------- | -------------------- |
| 1   | construct graph      | construct graph      |
| 2   | assemble fragments   | correct reads        |
| 3   | assemble transcripts | cluster reads        |
| 4   | N/A                  | assemble transcripts |

This is a very useful option if you only want to assemble fragments or correct long reads (ie. with `-stage 2`)!

### (D) list all available options in RNA-Bloom:
```
java -jar RNA-Bloom.jar -help
```

### (E) limit the size of Java heap:
```
java -Xmx2g -jar RNA-Bloom.jar ...
```
This limits the maximum Java heap to 2 GB with the `-Xmx` option. Note that `java` options has no effect on Bloom filter sizes.

See documentation for other [JVM options](https://docs.oracle.com/cd/E37116_01/install.111210/e23737/configuring_jvm.htm#OUDIG00071).


## Implementation :pencil:

RNA-Bloom is written in Java with Apache NetBeans IDE. It uses the [Apache Commons CLI library](https://commons.apache.org/proper/commons-cli/) and [JGraphT core library](https://jgrapht.org/).


## Citing RNA-Bloom :scroll:

If you use RNA-Bloom in your work, please cite us:

> Ka Ming Nip, Readman Chiu, Chen Yang, Justin Chu, Hamid Mohamadi, Rene L Warren, Inanc Birol. (2019) RNA-Bloom provides lightweight reference-free transcriptome assembly for single cells. bioRxiv 701607. doi: https://doi.org/10.1101/701607

--------------------------------------------------------------------------------

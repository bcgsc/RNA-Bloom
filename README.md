[![Release](https://img.shields.io/github/v/release/bcgsc/RNA-Bloom?include_prereleases)](https://github.com/bcgsc/RNA-Bloom/releases)
[![Downloads](https://img.shields.io/github/downloads/bcgsc/RNA-Bloom/total?logo=github)](https://github.com/bcgsc/RNA-Bloom/releases)
[![Conda](https://img.shields.io/conda/dn/bioconda/rnabloom?label=Conda)](https://anaconda.org/bioconda/rnabloom)

<p align="center">
  <img src="rnabloom_logo.png" alt="RNA-Bloom's logo"/>
</p>

**RNA-Bloom** is a fast and memory-efficient *de novo* transcript sequence assembler. It is designed for the following sequencing data types:
* single-end/paired-end bulk RNA-seq (strand-specific/agnostic)
* paired-end single-cell RNA-seq (strand-specific/agnostic)
* nanopore RNA-seq (cDNA/direct RNA)

Written by [Ka Ming Nip](mailto:kmnip@bcgsc.ca) :email:

:copyright: 2018-present Canada's Michael Smith Genome Sciences Centre, BC Cancer

--------------------------------------------------------------------------------

## Dependency :pushpin:

* [Java SE Development Kit (JDK) 11](https://www.oracle.com/java/technologies/downloads/#java11)

* External software used:

| software                                            | short reads            | long reads             |
| --------------------------------------------------- | ---------------------- | ---------------------- |
| [minimap2](https://github.com/lh3/minimap2)         | required               | required               |
| [Racon](https://github.com/lbcb-sci/racon)          | not used               | required               |
| [ntCard](https://github.com/bcgsc/ntCard) >=1.2.1   | required for `-ntcard` | required for `-ntcard` |

:warning: Their executables must be accessible from your `PATH`!



## Installation :wrench:

RNA-Bloom can be installed in two ways:

### (A) install with `conda`:
```
conda install -c bioconda rnabloom
```
All dependent software (listed above) will be installed. RNA-Bloom can be run as `rnabloom ...`

### (B) download from GitHub: 
1. Download the binary tarball `rnabloom_vX.X.X.tar.gz` from the [releases](https://github.com/bcgsc/RNA-Bloom/releases) section.
2. Extract the downloaded tarball with the command:
```
tar -zxf rnabloom_vX.X.X.tar.gz
```
RNA-Bloom can be run as `java -jar /path/to/RNA-Bloom.jar ...`



## Quick Start for Short Reads :running:

:warning: Input reads must be in either FASTQ or FASTA format and may be compressed with GZIP. 

### (A) assemble bulk RNA-seq data:

* paired-end reads only
  * when `left` reads are sense and `right` reads are antisense, use `-revcomp-right` to reverse-complement `right` reads
  * when `left` reads are antisense and `right` reads are sense, use `-revcomp-left` to reverse-complement `left` reads
  * for non-stranded data, use either `-revcomp-right` or `-revcomp-left`
```
java -jar RNA-Bloom.jar -left LEFT.fastq -right RIGHT.fastq -revcomp-right -ntcard -t THREADS -outdir OUTDIR
```

* single-end reads only
  * use `-sef` for forward reads and `-ser` for reverse reads
```
java -jar RNA-Bloom.jar -sef SE.fastq -ntcard -t THREADS -outdir OUTDIR
```

* paired-end and single-end reads
```
java -jar RNA-Bloom.jar -left LEFT.fastq -right RIGHT.fastq -revcomp-right -sef SE.fastq -ntcard -t THREADS -outdir OUTDIR
```

### (B) assemble multi-sample RNA-seq data with pooled assembly mode:
```
java -jar RNA-Bloom.jar -pool READSLIST.txt -revcomp-right -ntcard -t THREADS -outdir OUTDIR
```
This is especially useful for single-cell datasets. RNA-Bloom was tested on Smart-seq2 and SMARTer datasets.

#### file format for the `-pool` option:

This is a tabular file that describes the read file paths for all cells/samples to be used pooled assembly.
- Column header is on the first line, leading with `#`
- Columns are separated by space/tab characters
- Lines sharing the same `name` will be grouped together as the same sample during assembly

##### (i) paired-end reads only:
Only `name`, `left`, and `right` columns are specified for a total of 3 columns. The legacy header-less tri-column format is still supported.
```
#name left right
cell1 /path/to/cell1/left.fastq /path/to/cell1/right.fastq
cell2 /path/to/cell2/left.fastq /path/to/cell2/right.fastq
cell3 /path/to/cell3/left.fastq /path/to/cell3/right.fastq
```

##### (ii) paired and unpaired reads:
In addition to `name`, `left`, and `right` columns, either `sef`, `ser` or both are specified for a total of 4~5 columns. 
```
#name left right sef ser
cell1 /path/to/cell1/left.fastq /path/to/cell1/right.fastq /path/to/cell1/sef.fastq /path/to/cell1/ser.fastq
cell2 /path/to/cell2/left.fastq /path/to/cell2/right.fastq /path/to/cell2/sef.fastq /path/to/cell2/ser.fastq
cell3 /path/to/cell3/left.fastq /path/to/cell3/right.fastq /path/to/cell3/sef.fastq /path/to/cell3/ser.fastq
```

### (C) strand-specific assembly:
```
java -jar RNA-Bloom.jar -stranded ...
```
The `-stranded` option indicates that input reads are strand-specific.

Strand-specific reads are typically in the F2R1 orientation, where `/2` denotes *left* reads in *forward* orientation and `/1` denotes *right* reads in *reverse* orientation.

Configure the read file paths accordingly for bulk RNA-seq data and indicate read orientation:

`-stranded -left /path/to/reads_2.fastq -right /path/to/reads_1.fastq -revcomp-right`

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

### (A) assemble nanopore cDNA sequencing data:
```
java -jar RNA-Bloom.jar -long LONG.fastq -ntcard -t THREADS -outdir OUTDIR
```
Input reads are expected to be in a mix of both forward and reverse orientations.

### (B) assemble nanopore direct RNA sequencing data:
```
java -jar RNA-Bloom.jar -long LONG.fastq -stranded -ntcard -t THREADS -outdir OUTDIR
```
Input reads are expected to be only in the forward orientation.

By default, uracil (`U`) is written as `T`. Use the `-uracil` option to write `U` instead of `T` in the output assembly.

ntCard v1.2.1 supports uracil in reads.

### (C) assemble nanopore sequencing data with short-read polishing:
```
java -jar RNA-Bloom.jar -long LONG.fastq -sef SHORT.fastq -ntcard -t THREADS -outdir OUTDIR
```

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
| 3   | assemble transcripts | assemble transcripts |

This is a very useful option if you only want to assemble fragments or correct long reads (ie. with `-stage 2`)!

### (D) list all available options in RNA-Bloom:
```
java -jar RNA-Bloom.jar -help
```

### (E) limit the size of Java heap:
```
java -Xmx2g -jar RNA-Bloom.jar ...
```
or if you installed with `conda`:
```
export JAVA_TOOL_OPTIONS="-Xmx2g"
rnabloom ...
```
This limits the maximum Java heap to 2 GB with the `-Xmx` option. Note that `java` options has no effect on Bloom filter sizes.

See documentation for other [JVM options](https://docs.oracle.com/cd/E37116_01/install.111210/e23737/configuring_jvm.htm#OUDIG00071).


## Implementation :pencil:

RNA-Bloom is written in Java with Apache NetBeans IDE. It uses the [Apache Commons CLI library](https://commons.apache.org/proper/commons-cli/) and [JGraphT core library](https://jgrapht.org/).


## Citing RNA-Bloom :scroll:

If you use RNA-Bloom in your work, please cite [our publication](https://genome.cshlp.org/content/30/8/1191.full):

> Ka Ming Nip, Readman Chiu, Chen Yang, Justin Chu, Hamid Mohamadi, René L. Warren, and Inanc Birol. RNA-Bloom enables reference-free and reference-guided sequence assembly for single-cell transcriptomes. Genome Research. August 2020 30: 1191-1200; Published in Advance August 17, 2020, doi:10.1101/gr.260174.119

--------------------------------------------------------------------------------

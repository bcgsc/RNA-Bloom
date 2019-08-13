<p align="center">
  <img src="rnabloom_logo.png" alt="RNA-Bloom's logo"/>
</p>

**RNA-Bloom** is a fast and memory-efficient *de novo* transcript sequence assembler. It is designed for the following sequencing data types:
* paired-end bulk RNA-seq (strand-specific/agnostic)
* paired-end single-cell RNA-seq (strand-specific/agnostic)
* nanopore RNA-seq (cDNA/direct RNA)

See [Quick Start](#quick-start-running) for example usage.

Written by [Ka Ming Nip](mailto:kmnip@bcgsc.ca) :email:

:copyright: 2018 Canada's Michael Smith Genome Sciences Centre, BC Cancer

--------------------------------------------------------------------------------

## Dependency :pushpin:

* [Java SE Runtime Environment (JRE) 8](http://www.oracle.com/technetwork/java/javase/downloads/jre8-downloads-2133155.html)

Required for nanopore RNA-seq assembly:

* [minimap2](https://github.com/lh3/minimap2)
* [miniasm](https://github.com/lh3/miniasm)
* [Racon](https://github.com/isovic/racon)

Optional:

* [ntCard](https://github.com/bcgsc/ntCard)

Check your Java version:
```
java -version
```

Example:
```
java version "1.8.0_101"
Java(TM) SE Runtime Environment (build 1.8.0_101-b13)
Java HotSpot(TM) 64-Bit Server VM (build 25.101-b13, mixed mode)
```

## Installation :wrench:

1. Download the binary tarball `rnabloom_vX.X.X.tar.gz` from the [releases](https://github.com/bcgsc/RNA-Bloom/releases) section.
2. Extract the downloaded tarball with the command:
```
tar -zxf rnabloom_vX.X.X.tar.gz
```
3. RNA-Bloom is ready to use, ie. `java -jar /path/to/RNA-Bloom.jar ...`

**There is nothing to compile/configure/build!** :thumbsup:

## Quick Start :running:

:warning: RNA-Bloom only supports paired-end short-read RNA-seq data (eg. Illumina, BGISEQ) and nanopore RNA-seq data (Oxford Nanopore Technologies). Input reads must be in either FASTQ or FASTA format and may be compressed with GZIP. 

### assemble bulk RNA-seq data:
```
java -jar RNA-Bloom.jar -left LEFT.fastq.gz -right RIGHT.fastq.gz -revcomp-right -t THREADS -outdir OUTDIR
```

### assemble strand-specific bulk RNA-seq data:
```
java -jar RNA-Bloom.jar -stranded -left LEFT.fastq.gz -right RIGHT.fastq.gz -revcomp-right -t THREADS -outdir OUTDIR
```
Note that dUTP protocols produce reads in the F2R1 orientation, where `/2` denotes left reads in forward orientation and `/1` denotes right reads in reverse orientation. In this case, please specify your reads paths as `-left reads_2.fastq -right reads_1.fastq`.

### assemble single-cell RNA-seq data (Smart-seq2):
```
java -jar RNA-Bloom.jar -pool READSLIST.txt -revcomp-right -t THREADS -outdir OUTDIR
```

#### file format for the `-pool` option:

This text file is expected to have 3 columns, ie.

| column 1 | column 2           | column 3            |
| -------- | ------------------ | ------------------- |
| cell ID  | path of left reads | path of right reads |

Columns are separated by space/tab characters, eg.

```
cell1 /path/to/cell1/left.fastq.gz /path/to/cell1/right.fastq.gz
cell2 /path/to/cell2/left.fastq.gz /path/to/cell2/right.fastq.gz
cell3 /path/to/cell3/left.fastq.gz /path/to/cell3/right.fastq.gz
```

### assemble nanopore cDNA sequencing data:
```
java -jar RNA-Bloom.jar -ntcard -long READS.fa -t THREADS -outdir OUTDIR
```

### assemble nanopore direct RNA sequencing data:
```
java -jar RNA-Bloom.jar -stranded -ntcard -long READS.fa -t THREADS -outdir OUTDIR
```

### set the Bloom filter sizes based on the maximum allowable false positive rate and the expected number of unique k-mers:
```
java -jar RNA-Bloom.jar -fpr 0.05 -nk 28077715 ...
```
The number of unique k-mers in your dataset can be estimated efficiently with [ntCard](https://github.com/bcgsc/ntCard).

As an alternative to `-nk`, you can use the `-ntcard` option in RNA-Bloom if `ntcard` is already in your `PATH`, eg.

```
java -jar RNA-Bloom.jar -fpr 0.05 -ntcard ...
```
As a rule of thumb, a lower false positive rate may result in a better assembly but requires more memory for a larger Bloom filter.

### set the total size of Bloom filters to 3GB:
```
java -jar RNA-Bloom.jar -mem 3 ...
```
Otherwise, it is adjusted automatically based on the size of input read files.

### stop at an intermediate stage:
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

### list all available options in RNA-Bloom:
```
java -jar RNA-Bloom.jar -help
```

### limit the size of Java heap to 1GB:
```
java -Xmx1g -jar RNA-Bloom.jar ...
```
This option does not need to be set larger than the total Bloom filter size.

[Other JVM options](https://docs.oracle.com/cd/E37116_01/install.111210/e23737/configuring_jvm.htm#OUDIG00071) may also be used.


## Citing RNA-Bloom :scroll:

If you use RNA-Bloom in your work, please cite us:

> Ka Ming Nip, Readman Chiu, Chen Yang, Justin Chu, Hamid Mohamadi, Rene L Warren, Inanc Birol. (2019) RNA-Bloom provides lightweight reference-free transcriptome assembly for single cells. bioRxiv 701607. doi: https://doi.org/10.1101/701607

--------------------------------------------------------------------------------

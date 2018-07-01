<p align="center">
  <img src="rnabloom_logo.png" alt="RNA-Bloom's logo"/>
</p>

**RNA-Bloom** is a fast and memory-efficient *de novo* transcript sequence assembler for bulk and single cell paired-end RNA-seq data.

Written by [Ka Ming Nip](mailto:kmnip@bcgsc.ca)

Copyright 2018 Canada's Michael Smith Genome Sciences Centre, BC Cancer

--------------------------------------------------------------------------------

## Dependency

* [Java SE Runtime Environment 8](http://www.oracle.com/technetwork/java/javase/downloads/jre8-downloads-2133155.html)

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

## Installation

1. Download the binary tarball `rnabloom_vX.X.X.tar.gz` from the [releases](https://github.com/bcgsc/RNA-Bloom/releases) section.
2. Extract the downloaded tarball with the command:
```
tar -zxf rnabloom_vX.X.X.tar.gz
```
3. RNA-Bloom is ready to use, ie. `java -jar /path/to/RNA-Bloom.jar ...`
**There is nothing to compile/configure/build!**

## Quick Start

### assemble bulk RNA-seq data:
```
java -jar RNA-Bloom.jar -left LEFT.fastq.gz -right RIGHT.fastq.gz -revcomp-right -t THREADS -outdir OUTDIR
```

### assemble strand-specific RNA-seq data:
```
java -jar RNA-Bloom.jar -stranded -left LEFT.fastq.gz -right RIGHT.fastq.gz -revcomp-right -t THREADS -outdir OUTDIR
```
Note that dUTP protocols produce reads in the F2R1 orientation, where `/2` denotes left reads in forward orientation and `/1` denotes right reads in reverse orientation. In this case, please specify your reads paths as `-left reads_2.fastq -right reads_1.fastq`.

### assemble single cell RNA-seq data:
```
java -jar RNA-Bloom.jar -pool READSLIST.txt -revcomp-right -t THREADS -outdir OUTDIR
```

### example READSLIST.txt:
```
cell1 /path/to/cell1/left.fastq.gz /path/to/cell1/right.fastq.gz
cell2 /path/to/cell2/left.fastq.gz /path/to/cell2/right.fastq.gz
... ... ...
```
Columns are separated by space/tab characters.

This file consists of 3 columns, ie.
1. cell id
2. path of left reads
3. path of right reads

### limit the total size of Bloom filters to 3GB:
```
java -jar RNA-Bloom.jar -mem 3 ...
```
Otherwise, it is adjusted automatically based on the size of input read files.

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

--------------------------------------------------------------------------------
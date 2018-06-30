# RNA-Bloom

RNA-Bloom is a *de novo* transcript sequence assembler for bulk and single cell paired-end RNA-seq data.

Written by [Ka Ming Nip](mailto:kmnip@bcgsc.ca)

Copyright 2018 Canada's Michael Smith Genome Sciences Centre, BC Cancer

--------------------------------------------------------------------------------

## Quick Start

### bulk RNA-seq data:
```
java -jar RNA-Bloom.jar -left LEFT.fastq.gz -right RIGHT.fastq.gz -revcomp-right -t THREADS -outdir OUTDIR
```

### strand-specific RNA-seq data:
```
java -jar RNA-Bloom.jar -stranded -left LEFT.fastq.gz -right RIGHT.fastq.gz -revcomp-right -t THREADS -outdir OUTDIR
```

### single cell RNA-seq data:
```
java -jar RNA-Bloom.jar -pool READSLIST.txt -revcomp-right -t THREADS -outdir OUTDIR
```

### example READSLIST.txt:
```
cell1 /path/to/cell1/left.fastq.gz /path/to/cell1/right.fastq.gz
cell2 /path/to/cell2/left.fastq.gz /path/to/cell2/right.fastq.gz
```

### for a list of all available options:
```
java -jar RNA-Bloom.jar -help
```

--------------------------------------------------------------------------------
# RNA-Bloom
## *De novo* transcript sequence reconstruction from bulk and single cell RNA-seq data

[Ka Ming Nip](mailto:kmnip@bcgsc.ca)

Copyright 2018 Canada's Michael Smith Genome Sciences Centre, BC Cancer

--------------------------------------------------------------------------------

Sample command for bulk RNA-seq data:
```
java -jar RNA-Bloom.jar -left LEFT.fastq.gz -right RIGHT.fastq.gz -revcomp-right -t THREADS -outdir OUTDIR
```

Sample command for strand-specific RNA-seq data:
```
java -jar RNA-Bloom.jar -stranded -left LEFT.fastq.gz -right RIGHT.fastq.gz -revcomp-right -t THREADS -outdir OUTDIR
```

Sample command for single cell RNA-seq data:
```
java -jar RNA-Bloom.jar -pool READSLIST.txt -revcomp-right -t THREADS -outdir OUTDIR
```

Example READSLIST.txt:
```
cell1 /path/to/cell1/left.fastq.gz /path/to/cell1/right.fastq.gz
cell2 /path/to/cell2/left.fastq.gz /path/to/cell2/right.fastq.gz
```

For a list of other options:
```
java -jar RNA-Bloom.jar -help
```

--------------------------------------------------------------------------------
# From FASTQ to VCF file
This is a very short practical demo from the fastq files to VCF files.
Some commands continue on the next line. These commands have backslashes \ in the end of the line to indicate that the newline character (not shown) should be ignored. You should be able to copy (Ctrl+c) the full length command and then paste it in the console window. An alternative is to copy one line at time, leaving the backslashes out. In a script file, newline characters have to escaped with backslashes.

# What is FastQC
Modern high throughput sequencers can generate tens of millions of sequences in a single run. Before analysing this sequence to draw biological conclusions you should always perform some simple quality control checks to ensure that the raw data looks good and there are no problems or biases in your data which may affect how you can usefully use it.

Manual [FastQC](https://dnacore.missouri.edu/PDF/FastQC_Manual.pdf).</br>
Results [Good](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/good_sequence_short_fastqc.html) and [Bad](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/bad_sequence_fastqc.html).</br>

### Create a directory structure to work
```
mkdir dados
mkdir dados/fastq
mkdir dados/bwa
mkdir dados/fastqc
mkdir dados/freebayes
mkdir dados/gatk
```

## Get the FASTQ data
```
time cp /bioinfo/data/fastq/003.fastq.gz dados/fastq/
```


## Running the FASTQC analyses 
```
time fastqc -o dados/fastqc dados/fastq/003.fastq.gz
```

## Lets go to your home directory

```
cd ~/
name=sample_name;
library=library_name;
pl=miseq;

time bwa mem -M -R "@RG\tID:CAP\tSM:$name\tLB:$library\tPL:$pl" \
/bioinfo/referencia/hg19/chr1_13_17.fa \
dados/fastq/003.fastq.gz >dados/bwa/AMOSTRA01_S1.sam
```

## Converting SAM to BAM with samtools
```
time samtools fixmate dados/bwa/AMOSTRA01_S1.sam dados/bwa/AMOSTRA01_S1.bam
time samtools sort -O bam -o dados/bwa/AMOSTRA01_S1_sorted.bam dados/bwa/AMOSTRA01_S1.bam
time samtools index dados/bwa/AMOSTRA01_S1_sorted.bam
```


## Call variants assuming a diploid sample and require at least 0.3 fraction of observations supporting an alternate allele and also 15 supporting observations to consider a variant;

```
# -F --min-alternate-fraction N
#      Require at least this fraction of observations supporting
#      an alternate allele within a single individual in the
#      in order to evaluate the position.  default: 0.05
# -C --min-alternate-count N
#      Require at least this count of observations supporting
#      an alternate allele within a single individual in order
#      to evaluate the position.  default: 2

time /bioinfo/app/freebayes/bin/freebayes -f /bioinfo/referencia/hg19/chr1_13_17.fa \
-F 0.3 -C 15 \
--pooled-continuous dados/bwa/AMOSTRA01_S1_sorted.bam \
>dados/freebayes/AMOSTRA01_S1_sorted.vcf
```

## Variant-only calling on DNAseq
The HaplotypeCaller is capable of calling SNPs and indels simultaneously via local de-novo assembly of haplotypes in an active region. In other words, whenever the program encounters a region showing signs of variation, it discards the existing mapping information and completely reassembles the reads in that region. This allows the HaplotypeCaller to be more accurate when calling regions that are traditionally difficult to call, for example when they contain different types of variants close to each other.

```
time /bioinfo/app/gatk/gatk-4.1.2.0/gatk HaplotypeCaller -R /bioinfo/referencia/hg19/chr1_13_17.fa \
-I dados/bwa/AMOSTRA01_S1_sorted.bam \
-O dados/gatk/AMOSTRA01_S1_sorted.vcf
```

## References
https://github.com/lh3/bwa
https://github.com/broadinstitute/gatk
https://github.com/ekg/freebayes


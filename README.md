# From FASTQ to VCF
This is a very short practical demo from the fastq files to VCF files.
Some commands continue on the next line. These commands have backslashes \ in the end of the line to indicate that the newline character (not shown) should be ignored. You should be able to copy (Ctrl+c) the full length command and then paste it in the console window. An alternative is to copy one line at time, leaving the backslashes out. In a script file, newline characters have to escaped with backslashes.

# Análise de Qualidade das Sequências com o FastQC
Manual do [FastQC](https://dnacore.missouri.edu/PDF/FastQC_Manual.pdf).</br>
Exemplo de resultado [BOM](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/good_sequence_short_fastqc.html) e [RUIM](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/bad_sequence_fastqc.html).</br>

### Criar a estrutura de diretórios para trabalhar;
```
mkdir dados
mkdir dados/fastq
mkdir dados/bwa
mkdir dados/picard
mkdir dados/fastqc
mkdir dados/bedtools
mkdir dados/annovar
mkdir dados/freebayes
mkdir dados/gatk
```


## Copiar os FASTQ para sua pasta de análise;
```
time cp /bioinfo/data/fastq/003.fastq.gz dados/fastq/
```

## Listar os arquivos copiados;
```
ls -lh dados/fastq/*
```

## Executar o FASTQC para avaliar a qualidade das sequencias produzidas;
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

## Chamada de variantes com o GATK;
```
time /bioinfo/app/gatk/gatk-4.1.2.0/gatk HaplotypeCaller -R /bioinfo/referencia/hg19/chr1_13_17.fa \
-I dados/bwa/AMOSTRA01_S1_sorted.bam \
-O dados/gatk/AMOSTRA01_S1_sorted.vcf

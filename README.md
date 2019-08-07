# From FASTQ to VCF
This is a very short practical demo from the fastq files to VCF files.
Some commands continue on the next line. These commands have backslashes \ in the end of the line to indicate that the newline character (not shown) should be ignored. You should be able to copy (Ctrl+c) the full length command and then paste it in the console window. An alternative is to copy one line at time, leaving the backslashes out. In a script file, newline characters have to escaped with backslashes.
contra o hg19;
```
# Voltar para o HOME
cd ~/
NOME=NOME;
Biblioteca=Biblioteca;
Plataforma=Plataforma;

time bwa mem -M -R "@RG\tID:CAP\tSM:$NOME\tLB:$Biblioteca\tPL:$Plataforma" \
/bioinfo/referencia/hg19/chr1_13_17.fa \
dados/fastq/003.fastq.gz >dados/bwa/AMOSTRA01_S1.sam
```

## Utilizar o samtools: fixmate, sort e index
```
time samtools fixmate dados/bwa/AMOSTRA01_S1.sam dados/bwa/AMOSTRA01_S1.bam
time samtools sort -O bam -o dados/bwa/AMOSTRA01_S1_sorted.bam dados/bwa/AMOSTRA01_S1.bam
time samtools index dados/bwa/AMOSTRA01_S1_sorted.bam
```


## Chamada de variantes com o Freebayes;

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

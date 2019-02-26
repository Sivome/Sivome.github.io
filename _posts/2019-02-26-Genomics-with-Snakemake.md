

```python
# Author: Viswanadham Sridhara
# Date: 25 February 2019
# Title: Getting Started with Snakemake and Genomics

```


```python
As the name says i.e., Getting Started with Snakemake and Genomics, the snakemake rules described here are terse, but explicitly written to convey the information via *reading*

```


```python
Here, I started with a fastq file from NCBI SRA and with the fasta file, 1. built Bowtie2 index, 2. mapped fastq file to fasta to generate read alignments, 3. processed alignments
```


```python
To run the below pipeline, you need the NGS raw data (i.e., reads in fastq format)  and the fasta file of species of interest.
```


```python
# These SRRs are from the project description from NCBI portal
# https://www.ncbi.nlm.nih.gov/sra/SRX1981073[accn]
# Use sra-tools fasterq-dump to dump the SRRs used in the project
```


```python
# The fasta used in this analyses is from (name it to k12_profile.fasta)
# https://www.ncbi.nlm.nih.gov/genome?term=NC_000913.3&cmd=DetailsSearch
```


```python
# Rule 1: Build rule of all rules?
ecoli_samples = ['SRR5468393']

rule all:
  input: expand("{sample}.bed", sample=ecoli_samples)
```


```python
# Rule 2: build the index file from fasta file
rule bowtie_index_build:
    input:
        "k12_profile.fasta"
    shell:
        "bowtie2-build {input} k12Index/k12_profile"
```


```python
# Rule 3: Map fastq to fasta index
rule bowtie_map:
    input:
        expand("{sample}.fastq", sample=ecoli_samples)
    output:
        expand("{sample}.sam", sample=ecoli_samples)
    shell:
        "bowtie2 -x k12Index/k12_profile -p 4 --very-sensitive --no-unal -U {input} -S {output} -k 1"
```


```python
# Rule 4: Use samtools with q=30 for good alignments i.e., bam with a good cut-off
rule samtools_view:
    input:
        expand("{sample}.sam", sample=ecoli_samples)
    output:
        expand("{sample}.bam", sample=ecoli_samples)
    shell:
        "samtools view -q 30 -bS {input} > {output}"
```


```python
# Rule 5: Sort the bam
rule samtools_sort:
    input:
        expand("{sample}.bam", sample=ecoli_samples)
    output:
        expand("{sample}.Sorted", sample=ecoli_samples)
    shell:
        "samtools sort {input} > {output}"
```


```python
# Rule 6: convert bam to bed format (example)
rule bedtools_BAMtoBED:
    input:
        expand("{sample}.Sorted", sample=ecoli_samples)
    output:
        expand("{sample}.bed", sample=ecoli_samples)
    shell:
        "bedtools bamtobed -i {input} > {output}"
```


```python
Building DAG of jobs...
Using shell: /bin/bash
Provided cores: 1
Rules claiming more threads will be scaled down.
Job counts:
	count	jobs
	1	all
	1	bedtools_BAMtoBED
	1	bowtie_map
	1	samtools_sort
	1	samtools_view
	5

```


```python

[Mon Feb 25 21:05:14 2019]
rule bowtie_map:
    input: SRR5468393.fastq
    output: SRR5468393.sam
    jobid: 4

[Mon Feb 25 21:05:28 2019]
Finished job 4.
1 of 5 steps (20%) done

[Mon Feb 25 21:05:28 2019]
rule samtools_view:
    input: SRR5468393.sam
    output: SRR5468393.bam
    jobid: 3

[Mon Feb 25 21:05:32 2019]
Finished job 3.
2 of 5 steps (40%) done

[Mon Feb 25 21:05:32 2019]
rule samtools_sort:
    input: SRR5468393.bam
    output: SRR5468393.Sorted
    jobid: 2

[Mon Feb 25 21:05:35 2019]
Finished job 2.
3 of 5 steps (60%) done

[Mon Feb 25 21:05:35 2019]
rule bedtools_BAMtoBED:
    input: SRR5468393.Sorted
    output: SRR5468393.bed
    jobid: 1

[Mon Feb 25 21:05:38 2019]
Finished job 1.
4 of 5 steps (80%) done

[Mon Feb 25 21:05:38 2019]
localrule all:
    input: SRR5468393.bed
    jobid: 0

```

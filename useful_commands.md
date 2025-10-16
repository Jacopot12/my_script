# Command-Line Snippets Collection

## Sorting and Counting

#### Count specific lines

```bash
sort | uniq -n
```

#### Generate random numbers in a certain range

```bash
shuf -i 0-1 -n 1
```

#### Sort a BED file with two sub-genomes (e.g., Teff)

```bash
sort -k1,1V -k2,2n input.bed > output_sorted.bed
```

## Filtering and Selecting Lines

#### Get all lines that do not start with "#"

```bash
grep -v '^#' file.txt
```

#### Select lines that do not contain the word region

```bash
awk '!/region/' input.txt > output.txt
```

#### Select lines that do not contain a comma in column 2

```bash
awk '$2 !~ /,/' file.txt
```

## Genome and FASTA Operations

#### Generate genome size file

```bash
cut -f1,2 reference.fasta.fai > genome_file.genome
```

#### Remove a specific chr from a multi-fasta file

```bash
sed '/^>B73-Chr10/,/^>/d' ZmB73_chr5.fa > output.fa
```

#### Check genome size

```bash
zcat genome.fa.gz | grep -v ">" | tr -d '\n' | wc -m
```

## GTF and Feature Analysis

#### Check the number of genes in a GTF file

```bash
zcat file.gtf.gz | grep -v "#" | awk '$3=="gene"' | wc -l
```

#### Count every feature in a GTF file

```bash
zcat file.gtf.gz | grep -v "#" | cut -f3 | sort | uniq -c
```

## awk Operations

#### If the value in column $2 is less than 540 print the whole line

```bash
awk '($2 < 540) {print $0}' file.txt
```

#### Print a header and at the end print "DONE"

```bash
awk 'BEGIN{printf "Col1\tCol2\tCol3\n"} {print $1"\t"$2"\t"$3} END{print "DONE"}' file.txt
```

#### The -v flag specifies variables, FS stands for Field Separator

```bash
awk -v FS=: '{print $1}' file.txt
```

#### Print var1 = hello using a custom variable

```bash
awk -v var1=hello 'BEGIN{printf "var1 = %s\n", var1}'
```

#### Search and print the pattern /tty/

```bash
awk '/tty/ {print}' file.txt
```

#### The variable c counts matches

```bash
awk '/test/{++c} END {print "Total matched: ", c}' file.txt
```

#### Print lines shorter than 100 characters

```bash
awk 'length($0) < 100 {print}' file.txt
```

#### Replace a certain pattern in text and print everything

```bash
awk '{sub(/text to replace/, "replacement text"); print}'
```

#### Print lines only if the sum of column 1 and 2 is less than 10

```bash
awk '($1 + $2) < 10 {print}' test.txt
```

#### awk command to get only lines containing B73-Chr1 or B73-Chr2

```bash
awk '/B73-Chr1/ || /B73-Chr2/ {print}' file.txt > file_sub.txt
```

### awk command to get only lines where column1 is exactly chr1 or chr2

```bash
awk '$1 == "chr1" || $1 == "chr2"' file1.txt > file2.txt
```

#### Compare a single-column file with another single-column file and print matching lines

```bash
awk 'NR==FNR {a[$0]; next} $0 in a' file1.txt file2.txt > file3.txt
```

#### Convert a .vcf file to .bed file

```bash
awk -F'\t' '{if($0 !~ /^#/) print $1 "\t" $2-1 "\t" $2}' snp_file.vcf > snp_file.bed
```

#### Compare SNPs between two cultivars in my SNPs table

```bash
awk '{print $3, $4}' SNPs_file_onlyBiallelic.txt | awk '$1==$2' | wc -l
```

## File Format Conversions

#### Convert a sam file to bam file with samtools

```bash
samtools view -S -b input.sam > output.bam
```

#### Count reads in a fastq

```bash
echo $(zcat file.fastq.gz | wc -l)/4 | bc
```

## Fixing Files

#### "Fix" corrupted bam (without EOF block)

```bash
samtools view -@ 20 -b -o fixed.bam corrupted.bam
```
## Docker

#### Run a container

```bash
docker run --rm -v "$(pwd):/data" -u $(id -u):$(id -g) quay.io/biocontainers/ngmerge:0.3--0 NGmerge --help
```

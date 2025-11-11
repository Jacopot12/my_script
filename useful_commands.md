# Command-Line Snippets Collection

#### Basic bash command line:

**display processes**

```bash
top
```

**username**

```bash
whoami
```

**which users are logged in**

```bash
who
```

**disk free**

```bash
df -h
```

**disk usage**

```bash
du -h
```

---

#### Transfer and download data

wget can handle HTTP and FTP links

```bash
wget https://github.com/jacopoM28/CompOmics_2022/archive/refs/heads/main.zip
```

**-O saves the file with its original filename**

```bash
curl –O link
```

**copy a file from remote host to local host**

```bash
scp username@ip:path/to/file/to/copy /where/to/paste/it
```

**copy a file from local host to remote host**

```bash
scp path/to/file/to/copy username@ip:/where/to/paste/it
```

**copy a directory**

```bash
scp -r username@ip:path/to/file/to/copy /where/to/paste/it
```

---

#### Compress and decompress data

**compress file file.gz**

```bash
gzip file
```

**slower than gzip but higher compression ratio**

```bash
bzip2 file
```

**keep also the not compressed file**

```bash
gzip –k file
```

**uncompress file**

```bash
gunzip file.gz
```

**less compressed file**

```bash
zless file.gz
```

**use grep in compressed file**

```bash
zgrep "word" file.gz
```

**c create archive; z gzip archive; f specify new output; v verbose**

```bash
tar -zcvf myfolder.tar.gz myfolder
```

**decompress archive**

```bash
tar xvfz ./nome_archivio.tgz
```

---

#### Merge and sort files

**merge multiple files in 1**

```bash
cat file1 file2 file3 …
```

**sort the file, careful to computational sorting of file**

```bash
sort file
```

**human numeric sort**

```bash
sort –h file
```

**sort by first column adn then numerically by second column**

```bash
sort -k1,1 -k2,2n file
```

**sort by first column adn then numerically by second column in reversed order**

```bash
sort -k1,1 -k2,2nr file
```

**as before but human sorted**

```bash
sort -k1,1V -k2,2n file
```

**join two files according to first column of each file**

```bash
join -1 1 -2 1 sorted_file1 sorted_file2
```

**keep also non joined rows**

```bash
join -1 1 -2 1 -a 1 sorted_file1 sorted_file2
```

**merge lines of files**

```bash
paste file1 file2
```

---

#### Compare two sorted files

**Compare files line by line and show side by side**

```bash
diff -y file1 file2
```

**compare two sorted files line by line**

```bash
comm file1 file2
```

**lines unique to file1**

```bash
comm -1 file1 file2
```

**lines unique to file2**

```bash
comm -2 file1 file2
```

**print only lines present in both file1 and file2**

```bash
comm -12 file1 file2
```

**print lines in file1 not in file2, and vice versa**

```bash
comm -3 file1 file2
```

---

#### Grep

**print all rows that contains "word"**

```bash
grep "word" file
```

**print all rows that contains exactly the pattern "word"**

```bash
grep -w "word" file
```

**inverted match, print all rows that not contain the patter "word"**

```bash
grep -V "word" file
```

**ignore case distinctions, grep both "word" and "WORD"**

```bash
grep -i "word" file
```

**count how many rows contain the patter "word"**

```bash
grep -c "word" file
```

**print rows containing pattern "word" and the 10 rows after**

```bash
grep –A10 "word" file
```

**print rows containing pattern "word" and the 10 rows before**

```bash
grep –B10 "word" file
```

**print rows containing pattern "word" and the 10 rows after and before**

```bash
grep –C10 "word" file
```

**print Locus101 and Locus102**

```bash
grep "Locus10[12]" file
```

**print Locus101 and Locus102**

```bash
greo -E "(Locus101|Locus102)" file
```

---

#### Awk

**fiels separator is ";"**

```bash
awk -F";" '{print $3,$4}' file
```

**if first column = second column, print all columns**

```bash
awk '$1==$2 {print}' file
```

**if first column contain "chr2" or "chr3", print all columns, and a column with the difference between 3 and 2**

```bash
awk '$1 ~ /chr2|chr3/ { print $0 "\t" $3 - $2 }' file
```

**if both the first column contain "chr1" AND 3-2>0 , print all columns**

```bash
awk '$1 ~ /chr1/ && $3 - $2 > 10 '{print}' file
```

**print length instead of sequence in fasta file**

```bash
awk '{if ($1~">") print $0; else print length}' "fasta_file
```

---

#### Sed

**for each line subtitute "Locus" with "Transcripts" at first occurrance**

```bash
sed 's/Locus/Transcript/' file
```

**for each line subtitute "Locus" with "Transcripts" at first occurrance**

```bash
sed 's/Locus/Transcript/g' file
```

**overwrite input with the output**

```bash
sed -i 's/Locus/Transcript/g' file
```

**delete any row containing "Locus"**

```bash
sed '/Locus/d' file
```

---

#### Concatenate commands and programs

| connects the standard output of one process to the standard input of another

```bash
grep "word" file1 | sed 's/ /\t/g' | program1 > file2
```

; concatenate different commands or programs sequentially

```bash
grep ">" file1.fasta >output1 ; grep "_" file2.fasta output2
```

& Concatenate two programs so that program2 run only if program1 completed successfully

```bash
program1 input.txt > intermediate-results.txt && program2 intermediate-results.txt > results.txt
```

---

#### Standard output and standard error

```bash
program1 file 2> program1.stderr > results.txt
```

---

#### Command substitution

Command substitution runs a command inline and returns the output as a string that can be used for command.

```bash
echo "$(cat file.fasta)"
```

**show the string "There are 416 entries in my FASTA file."**

```bash
echo "There are $(grep -c '^>' input.fasta) entries in my FASTA file."
```

---

#### For loop

```bash
for i in *.fasta; do echo $i; done
```

```bash
for i in *.fasta; do mv $i ${i::-5}”_for_trimming”; done
```

```bash
for i in *.fasta; do mv $i ${i:1:3}”.fasta” ; done
```

```bash
for i in *fasta; do sed ‘s/>Locus/>/’ > $i”_editname” ; done
```

```bash
for i in *fasta; do grep –c”>” $i ; done > counts
```

```bash
for i in *fasta; do program1 $i > “output_”$i; done
```

```bash
for i in */ ; do cd $i; cp *.fasta ../; cd ..; done 
```

---

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

#### Rename headers in fastq.gz (bioawk)

```bash
zcat file_2.fastq.gz | bioawk -c fastx '{print "@" "SRR19651167." NR "/2"; print $seq; print "+"; print $qual}' | gzip > renamed_file_2.fastq.gz
```


## Docker

#### Run a container

```bash
docker run --rm -v "$(pwd):/data" -u $(id -u):$(id -g) quay.io/biocontainers/ngmerge:0.3--0 NGmerge --help
```

## Download reads

#### Download with fastq-dump (SRA toolkit)

```bash
fastq-dump --gzip --origfmt --split-files --skip-technical --defline-seq '@$ac.$sn.$si/$ri' SRRXXXXXXX
```





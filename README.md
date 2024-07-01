# my_scripts
This is a repository of the script that I generate for my bioinformatic purposes.

The first script, randomSNP_V2.sh, aims to put SNPs random in a fasta file (single line or multi-line). You can set a percentage of variation, from 0% to 100% of SNPs insertion.
Nevertheless this script is not able to assess if a base was already changed. As consequence, the percentage of SNPs insertion is the maximum variation that could be insert, but not in all cases.

The second script, randomSNP_V3.sh, corrects the problem of the randomSNP_V2.sh. Thus each base of the sequence can be changed only once, and the insertion rate of SNPs is effective.

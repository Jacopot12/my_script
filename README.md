# my_scripts
This is a repository of the script that I generate for my bioinformatic purposes.

The first, randomSNP_V2.sh, aim to put SNPs random in a file fasta (single line or multi-line). In this case you can set a percentage of variation, from 0% to 100% of SNPs insertion.
But, the problem in this script is that it is not able to assess if a base was already changed. As consequence, the percentage of SNPs insertion is the maximum variation that could be insert, baut not in al cases.

The second script, randomSNP_V3.sh, corrects the problem of changing the bases that were already changed. So every bases in the sequence could chance only once and the percentage of SNPs insertion is actual. 

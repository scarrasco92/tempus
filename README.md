# Tempus challenge data

This python script uses a vcf file and an annotated impact table from VEP.

In case of multiple consequences per variant the script prioritizes and picks the most impactful type based on the VEP annotation from https://m.ensembl.org/info/genome/variation/prediction/predicted_data.html

For ties, the script will report all cases from the most impactful category.

In the case of multiple variants at a single position, the script will evaluate each variant independently and report them in multiple lines

Usage:
python VCF_parsing.py

Output:
output_vcf_03_02.txt

Output columns:
chromosome, position, reference allele, alternate allele, depth of sequence coverage, number of reads for alternate, alternate allele fraction, allele frequency in exac, variant consequence, impact, biotype of consequence

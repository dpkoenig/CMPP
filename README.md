# Genotyping script for CMP Population
Python script to call genotypes from low coverage Illumina sequencing datset. Input is a BAM file containing mapped reads from a single individual and al list of sites with the columns:
* Chromosome
* Position
* Refence allele
* Alternative allele

Caller meant for SNPs only and returns single individual genotypes that need to be merged with additional data to form th vcf

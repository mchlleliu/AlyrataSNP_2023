## done on SNPs called from no-family pool

module load vcftools

tabix mpileup.vcf.gz

## rename SNP ID (non model sp gets dafaulted to ".", so plink can't filter out the variants by ID)
bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' mpileup.vcf.gz | bgzip > mpileup_NAMED.vcf.gz

## remove INDELs, non-biallelic SNPs, sites with more than 50% of 22 samples missing, and SNPs with MAF < 0.20
vcftools --gzvcf mpileup_NAMED.vcf.gz --max-missing 0.50 --min-alleles 2 --max-alleles 2 \
--remove-indels --maf 0.20 --recode --recode-INFO-all --stdout | bgzip -c > noFam_noIND_MAF20_MM50.vcf.gz
tabix noFam_noIND_MAF20_MM50.vcf.gz

###### ----------------- LD PRUNING -----------------------

## use plink to perform pairwise LD analysis
## the initial pruning is done to reduce the number of SNPs to a number that can be analysed
# perform LD r2 calculation for each pair of SNPs within a 200 kbp region, shifting the window 10 bp each time. 
# Separate SNP IDs into prune.in (r2 < 0.1) and prune.out (r2 > 0.1) 
plink --vcf noFam_noIND_MAF20_MM50.vcf.gz --make-founders --indep-pairwise 200 10 0.01 --out MAF20_MM50 --allow-extra-chr --const-fid





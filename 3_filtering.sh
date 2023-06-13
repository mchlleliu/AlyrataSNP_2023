## done on SNPs called from no-family pool

module load vcftools

tabix mpileup.vcf.gz

## rename SNP ID (non model sp gets dafaulted to ".", so plink can't filter out the variants by ID)
bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' mpileup.vcf.gz | bgzip > mpileup_NAMED.vcf.gz

## remove INDELs, non-biallelic SNPs, sites with more than 50% of 22 samples missing, and SNPs with MAF < 0.20
vcftools --gzvcf mpileup.named.vcf.gz --max-missing-count 15 \
--remove-indels --maf 0.20 --min-alleles 2 --max-alleles 2 --max-non-ref-af 0.999 \
--recode --recode-INFO-all --stdout | bgzip -c > noFam_noIND_MAF20.vcf.gz
tabix noFam_noIND_MAF20.vcf.gz

## compute p-val for exact tests of HWE
vcftools --gzvcf noFam_noIND_MAF20.vcf.gz --hardy --out candidates
head -n5 candidates.hwe
# CHR	POS	OBS(HOM1/HET/HOM2)	E(HOM1/HET/HOM2)	ChiSq_HWE	P_HWE	P_HET_DEFICIT	P_HET_EXCESS
# scaffold_1	6286	8/7/0	8.82/5.37/0.82	1.389414e+00	5.279693e-01	1.000000e+00	4.045977e-01
# scaffold_1	93783	10/7/1	10.12/6.75/1.12	2.469136e-02	1.000000e+00	7.355792e-01	7.403782e-01
# scaffold_1	101927	11/8/0	11.84/6.32/0.84	1.351111e+00	5.384006e-01	1.000000e+00	3.956567e-01
# scaffold_1	102272	10/7/0	10.72/5.56/0.72	1.142661e+00	1.000000e+00	1.000000e+00	4.627364e-01

## compute allele freq
vcftools --gzvcf noFam_noIND_MAF20.vcf.gz --freq --out candidates
head -n5 candidates.frq
# CHROM	POS	N_ALLELES	N_CHR	REF	REF_FREQ	ALT	ALT_FREQ
# scaffold_1	6286	2	30	G	0.766667	T	0.233333
# scaffold_1	93783	2	36	G	0.75	T	0.25
# scaffold_1	101927	2	38	G	0.789474	A	0.210526
# scaffold_1	102272	2	34	A	0.794118	T	0.205882


#### move .hwe and .frq files to R for duplicates filtering (see 3.5_ExcessHet_Filtering.Rmd)

## additional filtering for SNPs with only 15 samples genotyped
## at least 11 samples need to be successfully genotyped with at least 2 read depth per sample
vcftools --gzvcf noFam_noIND_MAF20.vcf.gz --positions Candidate15_ALL.pos \
--minGQ 15 --minDP 2 --max-missing 0.50 --recode --recode-INFO-all --stdout | bgzip -c > minDP2_cand15.vcf.gz
# get list of SNPs. download this and reload to R
bcftools view minDP2_cand15.vcf.gz --no-header | cut -f1-3 > minDP2.pos

###### ----------------- LD PRUNING -----------------------

## use plink to perform pairwise LD analysis
## the initial pruning is done to reduce the number of SNPs to a number that can be analysed
# perform LD r2 calculation for each pair of SNPs within a 200 kbp region, shifting the window 10 bp each time. 
# Separate SNP IDs into prune.in (r2 < 0.1) and prune.out (r2 > 0.1) 
plink --vcf noFam_noIND_MAF20_MM50.vcf.gz --make-founders --indep-pairwise 200 10 0.01 --out MAF20_MM50 --allow-extra-chr --const-fid





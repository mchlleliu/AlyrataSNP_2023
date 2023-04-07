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



#### ----------------------------- get 1000 bp flanking region fasta files for each SNPs -----------------------------------

## given a list of snps, extract from the annotated vcf file, and make a new vcf file 
bcftools view --include ID==@MAF20_MM50.prune.in noFam_noIND_MAF20_MM50.vcf.gz -Oz -o noFam_MAF20_MM50_pruned.vcf.gz

# get just the 75% sample coverage
vcftools --gzvcf noFam_MAF20_MM50_pruned.vcf.gz --max-missing 0.75 --recode --recode-INFO-all --stdout | bgzip -c > noFam_noIND_MAF20_MM75_pruned.vcf.gz

# extract regions and SNP id. Create tab-delimited file with the following format:
# CHR	POS	ID
zcat noFam_MAF20_MM75_pruned.vcf.gz | grep -v "^#" | cut -f1-3 > noFam_MAF20_MM75.pos

# calculate flanking region (sub 1000 with flanking region size in bp)
awk 'BEGIN{FS=OFS="\t"}; {start = $2 - 1000; end = $2 + 1000; print $1,start,end,$3}' noFam_MAF20_MM75.pos > noFam_MAF20_MM75_flank.bed

# bedfile looks like:
head -n5 noFam_MAF20_MM75_flank.bed
#scaffold_1	946617	948617	scaffold_1_947617_A_C
#scaffold_1	29965657	29967657	scaffold_1_29966657_T_C
#scaffold_1	29977522	29979522	scaffold_1_29978522_C_A
#scaffold_2	12544023	12546023	scaffold_2_12545023_C_A
#scaffold_2	14361014	14363014	scaffold_2_14362014_T_G


# change formatting to run for the script below
sed 's/\t/:/' noFam_MAF20_MM75_flank.bed | sed 's/\t/-/' >> tmp && mv noFam_MAF20_MM75_flank.bed

# positions list looks like:
head -n5 noFam_MAF20_MM75_flank.bed
#scaffold_1:946617-948617	scaffold_1_947617_A_C
#scaffold_1:29965657-29967657	scaffold_1_29966657_T_C
#scaffold_1:29977522-29979522	scaffold_1_29978522_C_A
#scaffold_2:12544023-12546023	scaffold_2_12545023_C_A
#scaffold_2:14361014-14363014	scaffold_2_14362014_T_G

## edit the bed file if the start position is <= 0 to 1 otherwise faidx can't work


# ----------------- script to get fasta file from bed files --------------------------
#!/bin/bash
#SBATCH --job-name=make_fasta
#SBATCH --mail-type=ALL
#SBATCH --mail-user=michelleliutbcs@gmail.com
#SBATCH -A rrg-shaferab
#SBATCH --cpus-per-task 32
#SBATCH --mem-per-cpu=3G
#SBATCH -t 0-2:59:00
#SBATCH -o %x-%j.log

module load samtools/1.15.1
module load bcftools

outputdir=/home/miliu/projects/def-shaferab/miliu/arabidopsis/cleandata/aligned/assembled/fasta
regions_file=/home/miliu/projects/def-shaferab/miliu/arabidopsis/cleandata/aligned/SNP/noFam_MAF20_MM75_flank.bed
reference=/home/miliu/projects/def-shaferab/miliu/arabidopsis/reference_genomes/index/Lyrata_reference.fa
vcf_file=/home/miliu/projects/def-shaferab/miliu/arabidopsis/cleandata/aligned/per_fam/mpileup_NAMED.vcf.gz

# generate flanking fasta files for each sample
# use the reference genome fasta, but feed it with the alternative allele of the sample if a SNP is present
# each fasta file would contain the flanking regions around the SNP of interest for each sample
while read -r line;
do
	region=$(echo "$line" | cut -f1)
	snp_name=$(echo "$line" | cut -f2)
	for j in $(cat popmap2.txt)
	do
		sample_name=$(echo $j)
		samtools faidx $reference $region | bcftools consensus --haplotype A -s $sample_name $vcf_file >> $outputdir/$snp_name.fa
	done
done < $regions_file


#------------------- Rename headers to sample names, add reference genome into the flank fasta files, and index final files ----------------
#!/bin/bash
#SBATCH --job-name=ren_addRef_idx
#SBATCH --mail-type=ALL
#SBATCH --mail-user=michelleliutbcs@gmail.com
#SBATCH -A rrg-shaferab
#SBATCH --cpus-per-task 32
#SBATCH --mem-per-cpu=3G
#SBATCH -t 0-2:59:00
#SBATCH -o %x-%j.log

module load samtools/1.15.1
module load bcftools

outputdir=/home/miliu/projects/def-shaferab/miliu/arabidopsis/cleandata/aligned/assembled/fasta
regions_file=/home/miliu/projects/def-shaferab/miliu/arabidopsis/cleandata/aligned/SNP/noFam_MAF20_MM75_flank.bed
reference=/home/miliu/projects/def-shaferab/miliu/arabidopsis/reference_genomes/index/Lyrata_reference.fa

# first, rename each fasta header to sample name
for i in $(ls *.fa)
do
awk 'NR==FNR{names[NR]=$0; next} /^>/{$1=">"names[++c]}1' ../popmap2.txt $i >> tmp && mv tmp $i # change name to sample names
done

# add the reference fasta
while read -r line;
do
	region=$(echo "$line" | cut -f1)
	snp_name=$(echo "$line" | cut -f2)
	samtools faidx $reference $region >> $outputdir/$snp_name.fa
done < $regions_file

for i in $(ls *.fa)
do
samtools faidx $i # index fasta files
done



## ----------------------------- CALC ALT AF and COVERAGE for a given SNP ----------------------------------

## calc.sh script to check Alternative Allele frequency (since biallelelic, MAF = min(Ref,Alt)) and coverage in both allSamples and noFam pools.
## use this in a directory containing your noFam and allSamples raw vcf files.
# run with = 
# bash calc.sh $scaffold_no $SNP_pos
### check stats

#!/bin/bash
set -e
set -u
set -o pipefail

no=$1
snp_pos=$2
snp_start=$(expr $snp_pos - 1)
all_vcf=./all_mpileup.vcf.gz  #direct to VCF containing all sample data
noFam_vcf=./noFam_mpileup.vcf.gz  #direct to VCF containing trimmed data

# view snps
ID=$(bcftools query --regions scaffold_$no:$snp_start-$snp_pos -f '%CHROM %POS %REF %ALT\n' mpileup.vcf.gz)
echo "SNP id: $ID"

# calc allele freq
AF_all=$(bcftools query --regions scaffold_$no:$snp_start-$snp_pos -f '%CHROM %POS %AN %AC{0}\n' $all_vcf | awk '{printf "%s %s %f\n",$1,$2,$4/$3}')
echo "Alternate allele freq for all samples $MAF_all"

AF_no=$(bcftools query --regions scaffold_$no:$snp_start-$snp_pos -f '%CHROM %POS %AN %AC{0}\n' $noFam_vcf | awk '{printf "%s %s %f\n",$1,$2,$4/$3}')
echo "Alternate allele freq for no fam samples $MAF_no"

# calc missing data
all_MM=$(bcftools query --regions scaffold_$no:$snp_start-$snp_pos -f '[\t%GT]\n' $all_vcf | awk '{OFS=RS;$1=$1}1' | grep -v "\." | wc -l)
echo "present in $all_MM out of 48"

fam_MM=$(bcftools query --regions scaffold_$no:$snp_start-$snp_pos -f '[\t%GT]\n' $noFam_vcf | awk '{OFS=RS;$1=$1}1' | grep -v "\." | wc -l)
echo "present in $fam_MM out of 22"


# prepare VCF of selected candidate SNPs 

# get the list of snps that are in both the 48 individuals and 22 no-fam groups after random LD pruning
awk 'BEGIN{FS=OFS="\t"};NR==FNR{a[$0];next} $1 in a {print $0}' afterQC.prune.in ../../all_samples/LD_30/afterQC.prune.in > cons_SNP.list

## given a list of snps, extract from the annotated vcf file, 
# you can skip this step if you have a .bim file (formatting of columns might be a bit different, but we can shift them around with awk)
bcftools view --include ID==@snp_id.list annotated_filtered.vcf.gz -Oz -o catalog_family.vcf.gz




# ---------------- this chunk is for the SNPs with 75% coverage. Same procedure, but change files names as necessary --------------


## given a list of snps, extract from the annotated vcf file, and make a new vcf file 
bcftools view --include ID==@MAF20_MM50.prune.in noFam_noIND_MAF20_MM50.vcf.gz -Oz -o noFam_MAF20_MM50_pruned.vcf.gz

# prepare VCF with just the SNPs with 75% sample coverage from the 22 indv
vcftools --gzvcf noFam_MAF20_MM50_pruned.vcf.gz --max-missing 0.75 --recode --recode-INFO-all --stdout | bgzip -c > noFam_noIND_MAF20_MM75_pruned.vcf.gz


# ---------------------------------------------------------------------------------------------------------------------------------



# extract regions and SNP id
zcat snps.vcf.gz | grep -v "^#" | cut -f1-3 > snps.pos
# format is CHR, POS, ID

# calculate flanking region (sub 1000 with flanking region size in bp)
awk 'BEGIN{FS=OFS="\t"}; {start = $2 - 1000; end = $2 + 1000; print $1,start,end,$3}' snps.pos > snp_flank_reg.bed

cat snp_flank_reg.bed
# bedfile looks like:
# scaffold_2	14982722	14984722	scaffold_2_14983722_T_C
# scaffold_4	4457011	4459011	scaffold_4_4458011_A_T
# scaffold_6	12174445	12176445	scaffold_6_12175445_C_T
# scaffold_10	106835	108835	scaffold_10_107835_A_T
# scaffold_11	133237	135237	scaffold_11_134237_C_T


# change formatting to run for the script below
# pos can only start from min 1.
sed 's/\t/:/' snp_flank_reg.bed | sed 's/\t/-/' >> tmp && mv tmp snp_flank_reg.bed
# positions list looks like:
# scaffold_2:14982722-14984722	scaffold_2_14983722_T_C
# scaffold_4:4457011-4459011	scaffold_4_4458011_A_T
# scaffold_6:12174445-12176445	scaffold_6_12175445_C_T
# scaffold_10:106835-108835	scaffold_10_107835_A_T
# scaffold_11:133237-135237	scaffold_11_134237_C_T


# get fasta file from bed files
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

outputdir=/home/miliu/projects/def-shaferab/miliu/arabidopsis/cleandata/aligned/assembled/fasta/ADD_MISS75_MAF20
regions_file=/home/miliu/projects/def-shaferab/miliu/arabidopsis/cleandata/aligned/SNP/snp_ADDMAF20_MM75.bed 
reference=/home/miliu/projects/def-shaferab/miliu/arabidopsis/reference_genomes/index/Lyrata_reference.fa
vcf_file=/home/miliu/projects/def-shaferab/miliu/arabidopsis/cleandata/aligned/per_fam/mpileup.norm.named.vcf.gz
# i might need to change the rest to follow the all-sample mpileup data

while read -r line;
do
	region=$(echo "$line" | cut -f1)
	snp_name=$(echo "$line" | cut -f2)
	for j in $(cat popmap2.txt)
	do
		sample_name=$(echo $j)
		# the --haplotype A option is for alternative allele. Better to use I for IUPAC so that you can see which samples are heterozygous
		# if you want to get just read sites, use -a N -M N to assign absent and missing sites as N 
		samtools faidx $reference $region | bcftools consensus --haplotype A -s $sample_name $vcf_file >> $outputdir/$snp_name.fa
	done
done < $regions_file


#---------------------------

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

outputdir=/home/miliu/projects/def-shaferab/miliu/arabidopsis/cleandata/aligned/assembled/fasta/ADD_MISS30
regions_file=/home/miliu/projects/def-shaferab/miliu/arabidopsis/cleandata/aligned/SNP/snp_flank_3030_ADD.bed
reference=/home/miliu/projects/def-shaferab/miliu/arabidopsis/reference_genomes/index/Lyrata_reference.fa

# first, rename each fasta header to sample name
for i in $(ls *.fa)
do
awk 'NR==FNR{names[NR]=$0; next} /^>/{$1=">"names[++c]}1' ../popmap2.txt $i >> tmp && mv tmp $i # change name
samtools faidx $i # index fasta files
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

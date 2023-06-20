## Make FASTA files of flanking regions 1000 bp up- and downstream of each candidate SNP site

# prepare VCF of selected candidate SNPs. Use positions file saved from R from excessHet filtering
vcftools --gzvcf MAF20_noDupsnoFixed.vcf.gz --positions noMaxMAF/candidateList.txt --recode --recode-INFO-all --stdout | bgzip -c > candidateList.vcf.gz
# set SNP id names
bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' candidateList.vcf.gz | bgzip > candidateList_NAMED.vcf.gz
tabix candidateList_NAMED.vcf.gz

# make positions and SNP id file
bcftools view --regions-file noMaxMAF/candidateList.txt candidateList_NAMED.vcf.gz --no-header | cut -f1-3 > candidateList.pos
head -n 5 candidateList.pos
#scaffold_1	93783	scaffold_1_93783_G_T
#scaffold_1	328079	scaffold_1_328079_C_T
#scaffold_1	531905	scaffold_1_531905_G_A
#scaffold_1	917625	scaffold_1_917625_T_A
#scaffold_1	947896	scaffold_1_947896_A_T

# make bed file of flanking regions around the candidate SNP
awk 'BEGIN{FS=OFS="\t"}; {start = $2 - 1000; end = $2 + 1000; print $1,start,end,$3}' candidateList.pos > candidateList.bed
head -n5 candidateList.bed
#scaffold_1	92783	94783	scaffold_1_93783_G_T
#scaffold_1	327079	329079	scaffold_1_328079_C_T
#scaffold_1	530905	532905	scaffold_1_531905_G_A
#scaffold_1	916625	918625	scaffold_1_917625_T_A
#scaffold_1	946896	948896	scaffold_1_947896_A_T

sed 's/\t/:/' candidateList.bed | sed 's/\t/-/' >> tmp && mv tmp candidateList.bed
head -n5 candidateList.bed
#scaffold_1:92783-94783	scaffold_1_93783_G_T
#scaffold_1:327079-329079	scaffold_1_328079_C_T
#scaffold_1:530905-532905	scaffold_1_531905_G_A
#scaffold_1:916625-918625	scaffold_1_917625_T_A
#scaffold_1:946896-948896	scaffold_1_947896_A_T


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

outputdir=/home/miliu/scratch/arabidopsisL/noDups/fasta/cov15
regions_file=/home/miliu/scratch/arabidopsisL/noDups/candidateList.bed
reference=/home/miliu/projects/def-shaferab/miliu/arabidopsis/reference_genomes/index/Lyrata_reference.fa
vcf_file=/home/miliu/projects/def-shaferab/miliu/arabidopsis/cleandata/aligned/all_samples/mpileup_NAMED.vcf.gz
popmap=/home/miliu/projects/def-shaferab/miliu/arabidopsis/cleandata/aligned/assembled/popmap.txt

# generate flanking fasta files for each sample
# use the reference genome fasta, but feed it with the alternative allele of the sample if a SNP is present
# each fasta file would contain the flanking regions around the SNP of interest for each sample
while read -r line;
do
	region=$(echo "$line" | cut -f1)
	snp_name=$(echo "$line" | cut -f2)
	for j in $(cat "$popmap")
	do
		sample_name=$(echo $j)
		samtools faidx $reference $region | bcftools consensus --iupac-codes -M N -s $sample_name $vcf_file >> $outputdir/$snp_name.fa
	done
done < $regions_file


##### end of script

# submit job :)
sbatch get_fasta.sh 


# ----------------- script to adjust FASTA header names --------------------------
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

regions_file=/home/miliu/scratch/arabidopsisL/noDups/candidateList.bed
reference=/home/miliu/projects/def-shaferab/miliu/arabidopsis/reference_genomes/index/Lyrata_reference.fa
popmap=/home/miliu/projects/def-shaferab/miliu/arabidopsis/cleandata/aligned/assembled/popmap.txt

# first, rename each fasta header to sample name
for i in $(ls *.fa)
do
awk 'NR==FNR{names[NR]=$0; next} /^>/{$1=">"names[++c]}1' $popmap $i >> tmp && mv tmp $i # change name
done


# add the reference fasta
while read -r line;
do
	region=$(echo "$line" | cut -f1)
	snp_name=$(echo "$line" | cut -f2)
	samtools faidx $reference $region >> $snp_name.fa
done < $regions_file

for i in $(ls *.fa)
do
samtools faidx $i # index fasta files
done


# in your .bam file folders
ls *.bam > popmap.txt

### popmap.txt contains all samples 
### popmap.txt file looks like:
head -n 5 popmap.txt 
#10_1_S7_L001.bam
#10_3_S8_L001.bam
#10_4_S9_L001.bam
#10_5_S10_L001.bam
#10_6_S11_L001.bam


#!/bin/bash
#SBATCH --job-name=bcf_caller
#SBATCH --mail-type=ALL
#SBATCH --mail-user=michelleliutbcs@gmail.com
#SBATCH -A rrg-shaferab
#SBATCH --cpus-per-task 32
#SBATCH --mem-per-cpu=3G
#SBATCH -t 0-11:59:00
#SBATCH -o %x-%j.log

module load bcftools

reference=/home/miliu/projects/def-shaferab/miliu/arabidopsis/reference_genomes/index/Lyrata_reference.fa
regions_dir=/home/miliu/projects/def-shaferab/miliu/arabidopsis/reference_genomes
outputdir=/home/miliu/projects/def-shaferab/miliu/arabidopsis/cleandata/aligned

# call multivariant alleles for nuclear genomes
# popmap2.txt is the list of bam files used.
# -R specifies the region (nuclear scaffolds only. skip the chloroplast region)
bcftools mpileup --threads 32 -Q 20 -q 20 -Ou \
--annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR \
-f $reference -R $regions_dir/list_nuc_sc.txt --bam-list popmap.txt | \
bcftools call --threads 32 -m -v -Oz -o $outputdir/mpileup/mpileup.vcf.gz



# call consensus alleles for each chloroplast genome
for i in $(ls *.bam | sed 's/.bam//')
do
bcftools mpileup --threads 32 -Q 20 -q 20 -Ou -R $regionsdir/list_chlor_sc.txt -f $reference $i.bam | \
bcftools call --threads 32 --ploidy 1 -c -Oz -o $outputdir/chlor_vcf/$i.vcf.gz
done

## script (for nuclear genome, chloroplast is only run once) is run again using modified popmap2.txt samples list, 
# consisting only of individuals from different maternal plants


## -------------------- 1. FASTQ QUALITY CONTROL -------------------
# pwd
# /home/scratch/arabidopsis/rawdata
mkdir quality

## generate fastQC reports
#!/bin/bash
#SBATCH --job-name=fastQC
#SBATCH --mail-type=ALL 
#SBATCH --mail-user=     		# tells slurm to send updates to your email
#SBATCH -A rrg-shaferab
#SBATCH --cpus-per-task 32 	
#SBATCH --mem-per-cpu=4G 		
#SBATCH -t 0-02:59:00 			
#SBATCH -o %x-%j.log 				

module load fastqc
fastqc -t 32 ./*fastq.gz -o ./quality

# -----------------------------------------

# generate MultiQC report for all files
cd quality/
module use /cvmfs/soft.mugqic/CentOS6/modulefiles
module load mugqic/MultiQC/1.12
multiqc --interactive .





#### ---------------------- 2. TRIMMING READS -------------------

########## RUNNING TRIMMOMATIC ################
#!/bin/bash
#SBATCH --job-name=trimmomatic
#SBATCH --mail-type=ALL
#SBATCH --mail-user=michelleliutbcs@gmail.com
#SBATCH -A rrg-shaferab
#SBATCH --cpus-per-task 32
#SBATCH --mem-per-cpu=4G
#SBATCH -t 0-02:59:00
#SBATCH -o %x-%j.log

module load trimmomatic/0.39

# make new directory for trimmomatic output
mkdir /home/miliu/projects/def-shaferab/miliu/arabidopsis/cleandata
# variable for output directory
outputdir=/home/miliu/projects/def-shaferab/miliu/arabidopsis/cleandata

# run trimmomatic for every individual
# paired end data, Nextera adapter
for i in $(ls *.fastq.gz | sed 's/_L001_R1_001.fastq.gz//' | sed 's/_L001_R2_001.fastq.gz//' | uniq)
do
  java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -threads 32 $i\_L001_R1_001.fastq.gz \
  $i\_L001_R2_001.fastq.gz $outputdir/$i\_L001_R1_001_paired.gz $outputdir/$i\_L001_R1_001_unpaired.gz \
  $outputdir/$i\_L001_R2_001_paired.gz $outputdir/$i\_L001_R2_001_unpaired.gz \
  ILLUMINACLIP:$EBROOTTRIMMOMATIC/adapters/NexteraPE-PE.fa:2:30:10
done




#### -------------------- 3. SEQUENCE ALIGNMENT -----------------------

### --------------- A) adjust reference file --------------------------

# download reference fasta file:
# A. lyrata subsp. lyrata strain MN47 (IDs: 1085921 [UID] 111458 [GenBank] 4425008 [RefSeq])
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/004/255/GCF_000004255.2_v.1.0/GCF_000004255.2_v.1.0_genomic.fna.gz > Lyrata_reference.fa

############# RENAME REF FASTA FILES ##################
# rename:
# >NW_003301862.1 Arabidopsis lyrata subsp. lyrata unplaced genomic scaffold ARALYscaffold_1118, whole genome shotgun sequence
# >NC_034379.1 Arabidopsis lyrata subsp. lyrata chloroplast DNA, complete genome, strain: MN47
### new name(s) = scaffold_n , chloroplast_genome

sed -i 's/.*ARALYscaffold_/>scaffold_/' Lyrata_nuclear.fa
sed -i 's/, whole.*//' Lyrata_nuclear.fa 
sed -i 's/>NC.*/>chloroplast_genome' Lyrata_nuclear.fa

faidx Lyrata_nuclear.fa




### ------------------- B) BWA RUN -----------------------------------

#!/bin/bash
#SBATCH --job-name=bwa
#SBATCH --mail-type=ALL
#SBATCH --mail-user=michelleliutbcs@gmail.com
#SBATCH -A rrg-shaferab
#SBATCH --cpus-per-task 32
#SBATCH --mem-per-cpu=4G
#SBATCH -t 0-02:59:00
#SBATCH -o %x-%j.log

module load bwa/0.7.17

output=/home/miliu/projects/def-shaferab/miliu/arabidopsis/cleandata/aligned

for i in $(ls *.gz | sed 's/_R1_001_paired.gz//' | sed 's/_R1_001_unpaired.gz//' | sed 's/_R2_001_paired.gz//' | sed 's/_R2_001_unpaired.gz//' | uniq)
do
bwa mem -t 32 /home/miliu/projects/def-shaferab/miliu/arabidopsis/reference_genomes/index/Lyrata_reference.fa \
$i\_R1_001_paired.gz $i\_R2_001_paired.gz > $output\/$i\.paired.sam
bwa mem -t 32 /home/miliu/projects/def-shaferab/miliu/arabidopsis/reference_genomes/index/Lyrata_reference.fa \
$i\_R1_001_unpaired.gz > $output\/$i\.unpaired.1.sam
bwa mem -t 32 /home/miliu/projects/def-shaferab/miliu/arabidopsis/reference_genomes/index/Lyrata_reference.fa \
$i\_R2_001_unpaired.gz > $output\/$i\.unpaired.2.sam
done



### -------------- C) POST-ALIGNMENT PROCESSING ----------------------------

# ---------------------- CONVERT SAM TO BAM --------------------------
#!/bin/bash
#SBATCH --job-name=samtobam
#SBATCH --mail-type=ALL
#SBATCH --mail-user=michelleliutbcs@gmail.com
#SBATCH -A rrg-shaferab
#SBATCH --cpus-per-task 32
#SBATCH --mem-per-cpu=3G
#SBATCH -t 0-02:59:00
#SBATCH -o %x-%j.log

module load StdEnv/2020 gcc/9.3.0
module load samtools/1.13

ln -s /home/miliu/projects/def-shaferab/miliu/arabidopsis/cleandata/aligned/*sam .

for i in $(ls *.sam | sed 's/.sam//')
do 
samtools view -S -b $i\.sam -@ 32 \
-T /home/miliu/projects/def-shaferab/miliu/arabidopsis/reference_genomes/index/Lyrata_reference.fa | \
samtools sort -@ 32 > \
/home/miliu/projects/def-shaferab/miliu/arabidopsis/cleandata/aligned/bam/$i.bam
done



# -------------------- INDEX, MERGE, INDEX ---------------------------
## index all the bam files

#!/bin/bash
#SBATCH --job-name=samindsep
#SBATCH --mail-type=ALL
#SBATCH --mail-user=michelleliutbcs@gmail.com
#SBATCH -A rrg-shaferab
#SBATCH --cpus-per-task 32
#SBATCH --mem-per-cpu=3G
#SBATCH -t 0-02:59:00
#SBATCH -o %x-%j.log

module load samtools/1.15.1

ln -s /home/miliu/projects/def-shaferab/miliu/arabidopsis/cleandata/aligned2/bam/*.bam .

for i in $(ls *.bam | sed 's/.bam//')
do
samtools index $i\.bam -@ 32 \
/home/miliu/projects/def-shaferab/miliu/arabidopsis/cleandata/aligned2/bam/$i.bai
done

## ------------- merge bam files for each sample and index merged reads
#!/bin/bash
#SBATCH --job-name=mergeindex
#SBATCH --mail-type=ALL
#SBATCH --mail-user=michelleliutbcs@gmail.com
#SBATCH -A rrg-shaferab
#SBATCH --cpus-per-task 32
#SBATCH --mem-per-cpu=3G
#SBATCH -t 0-02:59:00
#SBATCH -o %x-%j.log

outputdir=/home/miliu/projects/def-shaferab/miliu/arabidopsis/cleandata/aligned2/bam/assembled

module load samtools/1.15.1

for i in $(ls *.bam | sed 's/.paired.bam//' | sed 's/.unpaired.1.bam//' | sed 's/.unpaired.2.bam//' | uniq)
do
samtools merge -@ 32 -f -o $outputdir/$i.bam $i.paired.bam $i.unpaired.1.bam $i.unpaired.2.bam
done
# option -f
# Force to overwrite the output file if present.

cd ./assembled/

ln -s $outputdir/*bam .

for i in $(ls *.bam | sed 's/.bam//')
do
samtools index $i\.bam -@ 32 \
$outputdir/$i.bai
done


#### ------------------- 4. QUALITY CONTROL :) -------------------------
#!/bin/bash
#SBATCH --job-name=stats
#SBATCH --mail-type=ALL
#SBATCH --mail-user=michelleliutbcs@gmail.com
#SBATCH -A rrg-shaferab
#SBATCH --cpus-per-task 32
#SBATCH --mem-per-cpu=3G
#SBATCH -t 0-02:59:00
#SBATCH -o %x-%j.log

module load StdEnv/2020 gcc/9.3.0
module load samtools/1.13
module load bcftools/1.13

mkdir coverage flagstats depth

for i in $(ls *.bam | sed 's/.bam//')
do
samtools coverage --min-BQ 20 --min-MQ 20 -o ./coverage/$i.tsv $i\.bam
samtools depth $i.bam | awk '{total += $3 } END { print total/NR }' > ./depth/$i.txt
samtools flagstat -@ 32 -O tsv $i\.bam > ./flagstats/$i.tsv
done

cd ./coverage/

for i in $(ls *.tsv | sed 's/.tsv//')
do
# print specified field (column) of the document, then output it to a new .tmp doc
awk '{print $5}' $i.tsv > $i\_covered_bases.tmp
awk '{print $7}' $i.tsv > $i\_depth.tmp
awk '{print $8}' $i.tsv > $i\_baseQ.tmp
awk '{print $6}' $i.tsv > $i\_coverage.tmp
awk '{print $9}' $i.tsv > $i\_mapQ.tmp
done

paste *_covered_bases.tmp > covered_bases.txt
paste *_depth.tmp > depth.txt
paste *_baseQ.tmp > baseQ.txt
paste *_coverage.tmp > coverage.txt
paste *_mapQ.tmp > mapQ.txt
rm -rf *tmp
# merge all output into a single .txt file then delete the tmp files


cd ../depth/
cat *txt > whole_depths.txt #combine all the depth files


cd ../flagstats/

for i in $(ls *.tsv | sed 's/.tsv//')
do
awk '{print $1}' $i.tsv > $i.tmp
done
paste *tmp > stats.txt
rm -rf *tmp

## check R_PLOTS on how to plot coverage files!

# ------------------------------------------------------

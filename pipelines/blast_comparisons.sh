#!/bin/sh
#SBATCH --error blast_comparison
#SBATCH --output blast_comparison
#SBATCH --job-name blast_comparison
#SBATCH --mail-user amy.fitzpatrick@teagasc.ie
#SBATCH --mail-type END,FAIL
#SBATCH -p Priority,Background
#SBATCH -N 1
#SBATCH --cpus-per-task=5

cd data/blastdb
# files prepared in R script prepare_blastdb.R
# compare each expected simulation to observed simulation for each pipeline

# move expected fasta files to folder
mkdir expected_fasta
mv expected*.fasta expected_fasta/

# move expected fasta files to folder
mkdir observed_fasta
mv *.fasta observed_fasta/

# make blast db
module load blast/2.8.1
for i in expected_fasta/*.fasta;
do
name=$(basename "$i" .fasta| cut -d '_' -f2)
makeblastdb -in $i -dbtype nucl -out expected_${name}
done

#################################################################################################################################################################################################################### 
# check quality of datasets
mkdir blastdb_output
blastn -db expected_1 -query observed_fasta/AmpliTaq_Gold.SSII_1.fasta -out blastdb_output/AmpliTaq_Gold.SSII_1_blastout.txt -num_threads 5 -outfmt "6 qseqid sseqid evalue qcovs pident bitscore qlen" -perc_identity 99 -qcov_hsp_perc 75 
blastn -db expected_1 -query observed_fasta/Kapa_HiFi.SSII_1.fasta -out blastdb_output/Kapa_HiFi.SSII_1_blastout.txt -num_threads 5 -outfmt "6 qseqid sseqid evalue qcovs pident bitscore qlen" -perc_identity 99 -qcov_hsp_perc 75 
blastn -db expected_1 -query observed_fasta/Kapa_Robust.SSII_1.fasta -out blastdb_output/Kapa_Robust.SSII_1_blastout.txt -num_threads 5 -outfmt "6 qseqid sseqid evalue qcovs pident bitscore qlen" -perc_identity 99 -qcov_hsp_perc 75 

blastn -db expected_2 -query observed_fasta/AmpliTaq_Gold.SSII_2.fasta -out blastdb_output/AmpliTaq_Gold.SSII_2_blastout.txt -num_threads 5 -outfmt "6 qseqid sseqid evalue qcovs pident bitscore qlen" -perc_identity 99 -qcov_hsp_perc 75 
blastn -db expected_2 -query observed_fasta/Kapa_HiFi.SSII_2.fasta -out blastdb_output/Kapa_HiFi.SSII_2_blastout.txt -num_threads 5 -outfmt "6 qseqid sseqid evalue qcovs pident bitscore qlen" -perc_identity 99 -qcov_hsp_perc 75 
blastn -db expected_2 -query observed_fasta/Kapa_Robust.SSII_2.fasta -out blastdb_output/Kapa_Robust.SSII_2_blastout.txt -num_threads 5 -outfmt "6 qseqid sseqid evalue qcovs pident bitscore qlen" -perc_identity 99 -qcov_hsp_perc 75 

blastn -db expected_LP6 -query observed_fasta/LP6_HTS.fasta -out blastdb_output/LP6_HTS_blastout.txt -num_threads 5 -outfmt "6 qseqid sseqid evalue qcovs pident bitscore qlen" -perc_identity 99 -qcov_hsp_perc 75 
module unload blast/2.8.1
 
# merge all blastout files so results are combined from each simulation
awk 'FNR>1 || NR==1 {print $0","FILENAME}' blastdb_output/*blastout.txt  > 99_blastout.txt

# merge all blast fasta db files so you can compare identified vs missing OTUs
#for f in observed_fasta/*.fasta; do sed -i "s/^>/>${f%.fasta}_/g" "${f}"; done

cat expected_fasta/*.fasta > blastdb_expected.fasta

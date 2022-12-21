#!/bin/sh
#SBATCH --error err-unifrac.txt
#SBATCH --job-name unifrac
#SBATCH --mail-user amy.fitzpatrick@teagasc.ie
#SBATCH --mail-type END,FAIL
#SBATCH -p Background
#SBATCH -N 1
#SBATCH --cpus-per-task=5

cd data/unifrac_input

module load qiime2/2021.2
source activate qiime2-2021.2

for i in `find . -name "*unifrac_biom.tsv" -type f`; do

simulation=$(echo $i | awk -F_ '{print $1}')
echo $simulation 

biom convert \
  -i $i \
  -o ${simulation}_long_expected_observed.biom \
  --to-hdf5 --table-type="OTU table" 

qiime tools import \
   --type 'FeatureData[Sequence]' \
   --input-path ${simulation}_expected_observed.fasta \
   --output-path ${simulation}_sequences.qza

qiime tools import \
 --type FeatureTable[RelativeFrequency] \
 --input-path ${simulation}_long_expected_observed.biom \
 --input-format BIOMV210Format\
 --output-path ${simulation}_long_biom.qza

qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences ${simulation}_sequences.qza \
  --o-alignment ${simulation}_aligned-rep-seqs.qza \
  --o-masked-alignment ${simulation}_masked-aligned-rep-seqs.qza \
  --o-tree ${simulation}_unrooted-tree.qza \
  --o-rooted-tree ${simulation}_long_rooted-tree.qza

qiime diversity-lib unweighted-unifrac \
 --i-table ${simulation}_long_biom.qza \
 --i-phylogeny ${simulation}_long_rooted-tree.qza \
 --p-no-bypass-tips \
 --o-distance-matrix ${simulation}_unweighted_unifrac_long

qiime diversity-lib weighted-unifrac \
 --i-table ${simulation}_long_biom.qza \
 --i-phylogeny ${simulation}_long_rooted-tree.qza \
 --p-no-bypass-tips \
 --o-distance-matrix ${simulation}_weighted_unifrac_long

qiime diversity adonis --i-distance-matrix ${simulation}_unweighted_unifrac_long.qza \
 --m-metadata-file ${simulation}_long_metadata.tsv  \
 --p-formula source\
 --p-permutations 9999 \
 --o-visualization ${simulation}_unweighted_unifrac_adonis.qzv 
                                                        
qiime diversity adonis --i-distance-matrix ${simulation}_weighted_unifrac_long.qza \
 --m-metadata-file ${simulation}_long_metadata.tsv \
 --p-formula source \
 --p-permutations 9999 \
 --o-visualization ${simulation}_weighted_unifrac_adonis.qzv 

qiime tools export \
  --input-path ${simulation}_unweighted_unifrac_long.qza \
  --output-path ${simulation}_unweighted_unifrac_long

qiime tools export \
  --input-path ${simulation}_weighted_unifrac_long.qza \
  --output-path ${simulation}_weighted_unifrac_long

done

module unload qiime2/2021.2
conda deactivate



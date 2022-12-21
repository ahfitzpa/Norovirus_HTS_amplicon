
result=${PWD##*/}
# add library name to metadata and expected_composition
awk -v var="$result" 'BEGIN{FS=OFS="\t"}NR>1{$1=$1 "." var}1' metadata1.tsv > metadata.tsv
awk -v var="$result" 'BEGIN{FS=OFS="\t"}NR>1{$1=$1 " " var}1' expected_composition1.tsv > expected_composition2.tsv

#column headers
module load R/4.0.2
Rscript append_suffix_columnnames.R  
module unload R/4.0.2
#rm -r expected_composition1.tsv
rm -r expected_composition2.tsv










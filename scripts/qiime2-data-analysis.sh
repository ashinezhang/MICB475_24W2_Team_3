ed_##Following from qiime2-data-merging.sh and working in the merged_data directory
cd /data/diabetes/denoise_test/merged_data

# Generate a tree for phylogenetic diversity analyses
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences merged_rep-seqs.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza 

# Alpha-rarefaction - need to adjusted depth
qiime diversity alpha-rarefaction \
  --i-table merged_table.qza \
  --i-phylogeny rooted-tree.qza \
  --p-max-depth 25000 \
  --m-metadata-file /data/diabetes/new_merged_diabetes_metadata.tsv \
  --o-visualization alpha-rarefaction.qzv

  

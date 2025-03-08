##Following from qiime2-data-merging.sh and working in the merged_data directory
cd /data/diabetes/denoise_test/merged_data

# Generate a tree for phylogenetic diversity analyses
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences merged_rep-seqs.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza 

# Alpha-rarefaction
qiime diversity alpha-rarefaction \
  --i-table table.qza \
  --i-phylogeny rooted-tree.qza \
  --p-max-depth 8000 \
  --m-metadata-file /mnt/datasets/project_1/moving_pictures/sample-metadata.tsv \
  --o-visualization alpha-rarefaction.qzv

  

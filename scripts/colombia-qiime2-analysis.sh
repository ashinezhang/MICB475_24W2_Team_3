# We already have colombia_demux_seqs.qza, colombia-rep-seqs.qza, colombia-table.qza, colombia-stats.qza and the respective .qzv files completed 
# ^ Code is in the /scripts/qiime2-data-merging.sh file
# Now that we abandoned the Mexico diabetes dataset due to large regional differences, we need to do the analysis on just the Colombia dataset

# Mar 11

# make a new working directory to keep it separate from the merged analyses
mkdir /data/diabetes/colombia
cd /data/diabetes/colombia

## Import new_colombia_diabetes_metadata.tsv from local computer Github to the server
scp new_colombia_diabetes_metadata.tsv root@10.19.139.163:/data/diabetes/colombia

## Taxonomy Analysis

# Code takes a while to load, so make a new screen
screen -S colombia-only

##Skip training classifier and use provided classifier to assign taxonomy to your reads since they are V4 region
qiime feature-classifier classify-sklearn \
  --i-classifier /mnt/datasets/classifiers/silva-138-99-515-806-nb-classifier.qza \
  --i-reads /data/diabetes/denoise_test/colombia-rep-seqs.qza \
  --o-classification colombia_taxonomy.qza

qiime metadata tabulate \
 --m-input-file colombia_taxonomy.qza \
 --o-visualization colombia_taxonomy.qzv

## Create taxa barplot
qiime taxa barplot \
 --i-table /data/diabetes/denoise_test/colombia-table.qza \
 --i-taxonomy colombia_taxonomy.qza \
 --m-metadata-file /data/diabetes/colombia/new_colombia_diabetes_metadata.tsv \  
 --o-visualization colombia_taxa-bar-plots.qzv


###Filtering
#Filter out mitochondria and chloroplast data
qiime taxa filter-table \
  --i-table /data/diabetes/denoise_test/colombia-table.qza \
  --i-taxonomy colombia_taxonomy.qza \
  --p-exclude mitochondria,chloroplast \
  --o-filtered-table table-no-mitochondria-no-chloroplast.qza

#Visualization of table
qiime feature-table summarize \
  --i-table table-no-mitochondria-no-chloroplast.qza \
  --o-visualization table-no-mitochondria-no-chloroplast.qzv \
  --m-sample-metadata-file /data/diabetes/colombia/new_colombia_diabetes_metadata.tsv


# Generate a tree for phylogenetic diversity analyses
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences /data/diabetes/denoise_test/colombia-rep-seqs.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza 

# Alpha-rarefaction - need to adjusted depth
qiime diversity alpha-rarefaction \
  --i-table table-no-mitochondria-no-chloroplast.qza \
  --i-phylogeny rooted-tree.qza \
  --p-max-depth 20000 \
  --m-metadata-file /data/diabetes/colombia/new_colombia_diabetes_metadata.tsv\
  --o-visualization alpha-rarefaction.qzv

  ##Based on table and rare-faction: Sampling depth of 5900 was chosen to prevent loss of rarer groups

#Get Core Metrics
  qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted-tree.qza \
  --i-table  table-no-mitochondria-no-chloroplast.qza \
  --p-sampling-depth 5900 \
  --m-metadata-file /data/diabetes/colombia/new_colombia_diabetes_metadata.tsv \
  --output-dir core-metrics-results

  #Shannon Diversity qzv
  qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/shannon_vector.qza \
  --m-metadata-file  /data/diabetes/colombia/new_colombia_diabetes_metadata.tsv \
  --o-visualization core-metrics-results/shannon_significance.qzv

  #Beta Diversity
  qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file  /data/diabetes/colombia/new_colombia_diabetes_metadata.tsv \
  --m-metadata-column Group And Sex \
  --o-visualization core-metrics-results/weighted-unifrac-significance.qzv \
  --p-pairwise

### PICRUSt analysis on qiime2
# Mar 22, 2025
# make a new working directory to keep it separate from the merged analyses
mkdir /data/diabetes/picrust-analysis
cd /data/diabetes/picrust-analysis

# Filter the table.qza to remove all the features with 5 or lower counts to make PICRUSt2 analysis faster.
qiime feature-table filter-features \
  --i-table /data/diabetes/denoise_test/colombia-table.qza \
  --p-min-frequency 5 \
  --o-filtered-table feature-frequency-filtered-table.qza

# Run picrust2 

qiime picrust2 full-pipeline \
  --i-table feature-frequency-filtered-table.qza \
  --i-seq /data/diabetes/denoise_test/colombia-rep-seqs.qza \
  --output-dir q2-picrust2_output \
  --p-placement-tool sepp \
  --p-hsp-method pic \
  --p-max-nsti 2 \
  --verbose

# Convert the output files to human readable files
qiime tools export \
  --input-path q2-picrust2_output/pathway_abundance.qza \
  --output-path pathabun_exported

biom convert \
   -i pathabun_exported/feature-table.biom \
   -o pathabun_exported/pathway_abundance.tsv \
   --to-tsv

# Copy files to github
scp -r root@10.19.139.163:~/data/diabetes/picrust-analysis/pathabun_exported .
scp -r root@10.19.139.163:~/data/diabetes/picrust-analysis/q2-picrust2_output .

# Next, the data will be analysed on R

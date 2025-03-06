# Processing and merging Colombia and Mexico Diabetes datasets using QIIME2 for downstream analysis

## MARCH 1
# Create a directory for the project and navigate to it - this will be the working directory for separately processing the datasets 
mkdir /data/diabetes
cd /data/diabetes

# Import the Colombia dataset using the manifest file and SingleEnd code
qiime tools import \
  --type "SampleData[SequencesWithQuality]" \
  --input-format SingleEndFastqManifestPhred33V2 \
  --input-path /mnt/datasets/project_2/colombia/colombia_manifest.txt \
  --output-path ./colombia_demux_seqs.qza

# Import the Mexico Diabetes dataset using the manifest file and PairedEnd code
qiime tools import \
  --type "SampleData[PairedEndSequencesWithQuality]" \
  --input-format PairedEndFastqManifestPhred33V2 \
  --input-path /mnt/datasets/project_2/diabetes/mexico_manifest.tsv \
  --output-path ./mexico_demux_seqs.qza

## MARCH 3
# Create visualization of demultiplexed samples
qiime demux summarize \
  --i-data mexico_demux_seqs.qza \
  --o-visualization mexico_demux_seqs.qzv

qiime demux summarize \
  --i-data colombia_demux_seqs.qza \
  --o-visualization colombia_demux_seqs.qzv

# Make directory in home directory on local machine and copy file
scp root@<IP>:/data/diabetes/mexico_demux_seqs.qzv .
scp root@<IP>:/data/diabetes/colombia_demux_seqs.qzv .
# use https://view.qiime2.org/ to visualize and chose truncation of 220

## MARCH 4, code took very long to load
# Denoise: Determine ASVs with DADA2
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs mexico_demux_seqs.qza \
  --p-trim-left-f 0 \
  --p-trunc-len-f 220 \
  --p-trim-left-r 0 \
  --p-trunc-len-r 220 \
  --o-representative-sequences mexico-rep-seqs.qza \
  --o-table mexico-table.qza \
  --o-denoising-stats mexico-stats.qza

qiime dada2 denoise-single \
  --i-demultiplexed-seqs colombia_demux_seqs.qza \
  --p-trim-left 0 \
  --p-trunc-len 220 \
  --o-representative-sequences colombia-rep-seqs.qza \
  --o-table colombia-table.qza \
  --o-denoising-stats colombia-stats.qza
  
# Visualize DADA2 stats
qiime metadata tabulate \
  --m-input-file mexico-stats.qza \
  --o-visualization mexico-stats.qzv

qiime metadata tabulate \
  --m-input-file colombia-stats.qza \
  --o-visualization colombia-stats.qzv

# Visualize ASVs stats
qiime feature-table summarize \
  --i-table mexico-table.qza \
  --o-visualization mexico-table.qzv \
  --m-sample-metadata-file /mnt/datasets/project_2/diabetes/mexico_metadata.tsv

qiime feature-table summarize \
  --i-table colombia-table.qza \
  --o-visualization colombia-table.qzv \
  --m-sample-metadata-file /mnt/datasets/project_2/colombia/colombia_metadata.txt

#Match 5th This code does not run. Error: (1/1) Invalid value for '--i-data': mexico-rep-seqs.qza is not a QIIME archive. 
#Some source said to update qiime

qiime feature-table tabulate-seqs \
  --i-data mexico-rep-seqs.qza \
  --o-visualization mexico-rep-seqs.qzv

qiime feature-table tabulate-seqs \
  --i-data colombia-rep-seqs.qza \
  --o-visualization colombia-rep-seqs.qzv

###David thinks the previous denoising code didn't complete so I am going to trymaking a new directory and denoising again
mkdir denoise_test
cd denoise_test
screen -S denoise_mexico 

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs /data/diabetes/mexico_demux_seqs.qza \
  --p-trim-left-f 0 \
  --p-trunc-len-f 220 \
  --p-trim-left-r 0 \
  --p-trunc-len-r 220 \
  --o-representative-sequences mexico-rep-seqs.qza \
  --o-table mexico-table.qza \
  --o-denoising-stats mexico-stats.qza

screen -S denoise_colombia

  qiime dada2 denoise-single \
  --i-demultiplexed-seqs colombia_demux_seqs.qza \
  --p-trim-left 0 \
  --p-trunc-len 220 \
  --o-representative-sequences colombia-rep-seqs.qza \
  --o-table colombia-table.qza \
  --o-denoising-stats colombia-stats.qza


  ##Repeat on new files
  # Visualize DADA2 stats
qiime metadata tabulate \
  --m-input-file mexico-stats.qza \
  --o-visualization mexico-stats.qzv

qiime metadata tabulate \
  --m-input-file colombia-stats.qza \
  --o-visualization colombia-stats.qzv

# Visualize ASVs stats
qiime feature-table summarize \
  --i-table mexico-table.qza \
  --o-visualization mexico-table.qzv \
  --m-sample-metadata-file /mnt/datasets/project_2/diabetes/mexico_metadata.tsv

qiime feature-table summarize \
  --i-table colombia-table.qza \
  --o-visualization colombia-table.qzv \
  --m-sample-metadata-file /mnt/datasets/project_2/colombia/colombia_metadata.txt

qiime feature-table tabulate-seqs \
  --i-data mexico-rep-seqs.qza \
  --o-visualization mexico-rep-seqs.qzv

qiime feature-table tabulate-seqs \
  --i-data colombia-rep-seqs.qza \
  --o-visualization colombia-rep-seqs.qzv


## Merge Datasets Together
# Make and change directory to /data/diabetes/merged_data, this will be the working directory for the merged datasets
mkdir merged_data
cd merged_data

# Merge tables
qiime feature-table merge \
 --i-tables /data/diabetes/denoise_test/mexico-table.qza \
 --i-tables /data/diabetes/denoise_test/colombia-table.qza \
 --o-merged-table merged_table.qza

# Merge rep-seqs
qiime feature-table merge-seqs \
 --i-data /data/diabetes/denoise_test/mexico-rep-seqs.qza \
 --i-data /data/diabetes/denoise_test/colombia-rep-seqs.qza \
 --o-merged-data merged_rep-seqs.qza


#######Skip training classifier and use exisiting one
## Train classifier on the merged datasets because there are only primer sequences for Colombia dataset and they are both V4 region 
  # CODE FROM CANVAS, NEED TO EDIT
  # Extract your amplicon of interest from the reference database
  # replace the ref-otus.qza with the representative sequence file on the server (e.g. /mnt/datasets/silva_ref_files/silva-138-99-seqs.qza)
  # replace primer sequences with your correct sequences
  # replace trunc-len with the one you defined in your denoising step
  qiime feature-classifier extract-reads \
    --i-sequences ref-otus.qza \ # merged_rep-seqs.qza?
    --p-f-primer GTGCCAGCMGCCGCGGTAA \
    --p-r-primer GGACTACHVGGGTWTCTAAT \
    --p-trunc-len 220 \
    --o-reads ref-seqs-trimmed.qza
  
  # Train classifier with your new ref-seq file
  # Replace ref-taxonomy.qza with the representative taxonomy file on the server (e.g. /mnt/datasets/silva_ref_files/silva-138-99-tax.qza)
  qiime feature-classifier fit-classifier-naive-bayes \
    --i-reference-reads ref-seqs-trimmed.qza \
    --i-reference-taxonomy ref-taxonomy.qza \ # i have no idea what this would be, don't we need the classifer to make the taxonomy?
    --o-classifier classifier.qza
  
  # Use the trained classifier to assign taxonomy to your reads (rep-seqs.qza)
  qiime feature-classifier classify-sklearn \
    --i-classifier classifier.qza \
    --i-reads merged_rep-seqs.qza \
    --o-classification merged_taxonomy.qza
   ######################


## Taxonomy Analysis
##Skip training classifier and use provided classifier:
qiime feature-classifier classify-sklearn \
  --i-classifier /mnt/datasets/classifiers/silva-138-99-515-806-nb-classifier.qza \
  --i-reads merged_rep-seqs.qza \
  --o-classification merged_taxonomy.qza

  qiime metadata tabulate \
  --m-input-file merged_taxonomy.qza \
  --o-visualization merged_taxonomy.qzv

  qiime taxa barplot \
  --i-table merged_table.qza \
  --i-taxonomy merged_taxonomy.qza \
  --m-metadata-file /data/metadata_files/team4_metadata_2.0.tsv \  ##Change for our metadata file, may need to import to server
  --o-visualization merged_taxa-bar-plots.qzv

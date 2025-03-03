# Processing and merging Colombia and Mexico Diabetes datasets using QIIME2 for downstream analysis

# MARCH 1
# Create a directory for the project and navigate to it - this will be the working directory
mkdir /data/diabetes
cd /data/diabetes

# Import the Colombia dataset using the manifest file
qiime tools import \
  --type "SampleData[SequencesWithQuality]" \
  --input-format SingleEndFastqManifestPhred33V2 \
  --input-path /mnt/datasets/project_2/colombia/colombia_manifest.txt \
  --output-path ./colombia_demux_seqs.qza

# Import the Mexico Diabetes dataset using the manifest file
qiime tools import \
  --type "SampleData[PairedEndSequencesWithQuality]" \
  --input-format PairedEndFastqManifestPhred33V2 \
  --input-path /mnt/datasets/project_2/diabetes/mexico_manifest.tsv \
  --output-path ./mexico_demux_seqs.qza

# MARCH 3
# Create visualization of demultiplexed samples
qiime demux summarize \
  --i-data mexico_demux_seqs.qza \
  --o-visualization mexico_demux_seqs.qzv

qiime demux summarize \
  --i-data colombia_demux_seqs.qza \
  --o-visualization colombia_demux_seqs.qzv
  
# Make directory in home directory on local machine
# Copy file
scp root@<IP>:/data/diabetes/mexico_demux_seqs.qzv . 
scp root@<IP>:/data/diabetes/colombia_demux_seqs.qzv . 

# chose truncation of 220

# MARCH ??
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
  
qiime feature-table tabulate-seqs \
  --i-data mexico-rep-seqs.qza \
  --o-visualization mexico-rep-seqs.qzv

qiime feature-table tabulate-seqs \
  --i-data colombia-rep-seqs.qza \
  --o-visualization colombia-rep-seqs.qzv

  # Maybe use this to split forward and reverse for Colombia
  qiime demux split-paired \
  --i-demux demux-paired.qza \
  --o-demux-forward only-forward-reads.qza \
  --o-demux-reverse only-reverse-reads.qza


mkdir colombia_export
cd colombia_export

#Copy and export metadata, OTU table, Taxonomy, and Phylogenetic tree into one folder
cp ../fixed_colombia_diabetes_metadata.tsv metadata.tsv

qiime tools export \
  --input-path ../colombia_taxonomy.qza \
  --output-path taxonomy_export
qiime tools export \
  --input-path ../rooted-tree.qza \
  --output-path rooted-tree.nwk

qiime tools export \
  --input-path ~/data/diabetes/colombia-table.qza \
  --output-path table_export

cd table_export
biom convert -i feature-table.biom --to-tsv -o feature-table.txt

#move to personal computer
scp -r root@10.19.139.163:~/data/diabetes/colombia/colombia_export .

#Clean up and move all files into one directory
move table_export\* .
move taxonomy_export\* .
move tree_export\* .
rmdir table_export
rmdir tree_export
rmdir taxonomy_export

# Directory for QIIME2 scripts

Script used so far

###import data in qiime and demultiplex
qiime tools import \
  --type "SampleData[SequencesWithQuality]" \
  --input-format SingleEndFastqManifestPhred33V2 \
  --input-path parkinsons_manifest.txt \
  --output-path ./demux.qza
###visualize demultiplexed data
qiime demux summarize \
  --i-data demux.qza \
  --o-visualization demux.qzv

### denoise
qiime quality-filter q-score \
 --i-demux demux.qza \
 --o-filtered-sequences demux-filtered.qza \
 --o-filter-stats demux-filter-stats.qza

qiime deblur denoise-16S \
  --i-demultiplexed-seqs demux-filtered.qza \
  --p-trim-length 251 \
  --o-representative-sequences rep-seqs.qza \
  --o-table table.qza \
  --p-sample-stats \
  --o-stats deblur-stats.qza

### Taxonomic analysis
qiime feature-classifier classify-sklearn \
  --i-classifier /mnt/datasets/classifiers/silva-138-99-515-806-nb-classifier.qza \
  --i-reads rep-seqs.qza \
  --o-classification taxonomy.qza
  
### fitler data to only healthy patients and no mito/chloro
qiime feature-table filter-samples \
  --i-table table.qza \
  --m-metadata-file parkinsons_metadata.txt \
  --p-where "[Disease]='Control'" \
  --o-filtered-table control-table.qza
  
qiime taxa filter-table \
  --i-table control-table.qza \
  --i-taxonomy taxonomy.qza \
  --p-exclude mitochondria,chloroplast \
  --o-filtered-table table-no-mitochondria-no-chloroplast.qza

### visualize data after denoising and filtering
qiime feature-table summarize \
  --i-table table-no-mitochondria-no-chloroplast.qza \
  --o-visualization table-no-mitochondria-no-chloroplast.qzv \
  --m-sample-metadata-file parkinsons_metadata.txt 

### Generate a tree for phylogenetic diversity analyses
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza 

### Alpha-rarefaction
qiime diversity alpha-rarefaction \
  --i-table table-no-mitochondria-no-chloroplast.qza \
  --i-phylogeny rooted-tree.qza \
  --p-max-depth 14000 \
  --m-metadata-file parkinsons_metadata.txt \
  --o-visualization alpha-rarefaction.qzv

qiime feature-table summarize \
  --i-table table.qza \
  --o-visualization table.qzv \
  --m-sample-metadata-file parkinsons_metadata.txt 

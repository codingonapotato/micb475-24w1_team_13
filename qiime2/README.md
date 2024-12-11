# Directory for QIIME2 scripts

Script used so far

### Import data in qiime and demultiplex
```
qiime tools import \
  --type "SampleData[SequencesWithQuality]" \
  --input-format SingleEndFastqManifestPhred33V2 \
  --input-path filtered_parkinsons_manifest.txt \
  --output-path ./demux.qza
```
### Visualize demultiplexed data
```
qiime demux summarize \
  --i-data demux.qza \
  --o-visualization demux.qzv
```
### Denoise
```
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
```
### Taxonomic analysis
```
qiime feature-classifier classify-sklearn \
  --i-classifier /mnt/datasets/classifiers/silva-138-99-515-806-nb-classifier.qza \
  --i-reads rep-seqs.qza \
  --o-classification taxonomy.qza
```
### Filter data to no mito/chloro
```
qiime taxa filter-table \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --p-exclude mitochondria,chloroplast \
  --o-filtered-table table-no-mitochondria-no-chloroplast.qza
```
### Visualize data after denoising and filtering
```
qiime feature-table summarize \
  --i-table table-no-mitochondria-no-chloroplast.qza \
  --o-visualization table-no-mitochondria-no-chloroplast.qzv \
  --m-sample-metadata-file pd_filtered_metadata.txt 
```
### Generate a tree for phylogenetic diversity analyses
```
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza 
```
### Alpha-rarefaction
```
qiime diversity alpha-rarefaction \
  --i-table table-no-mitochondria-no-chloroplast.qza \
  --i-phylogeny rooted-tree.qza \
  --p-max-depth 13000 \
  --m-metadata-file pd_filtered_metadata.txt \
  --o-visualization alpha-rarefaction.qzv

```
### Core metrics
```
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted-tree.qza \
  --i-table table-no-mitochondria-no-chloroplast.qza \
  --p-sampling-depth 4709 \
  --m-metadata-file pd_filtered_metadata.txt \
  --output-dir core-metrics-results

# Calculate alpha-group-significance
qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/faith_pd_vector.qza \
  --m-metadata-file pd_filtered_metadata.txt \
  --o-visualization core-metrics-results/faith-pd-group-significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/evenness_vector.qza \
  --m-metadata-file pd_filtered_metadata.txt \
  --o-visualization core-metrics-results/evenness-group-significance.qzv
  
# Calculate beta-group-significance
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file pd_filtered_metadata.txt \
  --m-metadata-column BDI_category_antidepressant_use \
  --o-visualization core-metrics-results/unweighted-unifrac-BA-significance.qzv \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file pd_filtered_metadata.txt \
  --m-metadata-column Antidepressant_use \
  --o-visualization core-metrics-results/unweighted-unifrac-A-significance.qzv \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file pd_filtered_metadata.txt \
  --m-metadata-column BDI_category \
  --o-visualization core-metrics-results/unweighted-unifrac-B-significance.qzv \
  --p-pairwise

```

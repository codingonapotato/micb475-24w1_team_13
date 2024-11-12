# Directory for QIIME2 scripts

Script used so far

### import data in qiime and demultiplex
qiime tools import \
  --type "SampleData[SequencesWithQuality]" \
  --input-format SingleEndFastqManifestPhred33V2 \
  --input-path parkinsons_manifest.txt \
  --output-path ./demux.qza
### visualize demultiplexed data
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
  --p-where "[Disease]='PD'" \
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


# A2/A3 Qiime2 Workflow
  > Prequisites: Run `categorize_metadata_by_bdi.R` to generate an updated metadata file with new column for BDI category ('low' or 'high')

  > Have to restart from filter step otherwise qiime feature-table will throw error:

```
  Plugin error from feature-table:

  The following IDs are not present in the metadata: 'F1195', 'F1256', 'F1262', 'F1282', 'F1290', 'F1320', 'F1385', 'F1909', 'F1915', 'F1920', 'F1975', 'F1977', 'F2113', 'F2118', 'F2281', 'F2366', 'F2395', 'F2501', 'F2504', 'F2507', 'F2516', 'F2575', 'F2716', 'F2881', 'F2953', 'F2956', 'F2971', 'F3188', 'F3211', 'F3263', 'F3344', 'F3395', 'F3436', 'F3484', 'F3560', 'F3580', 'F3627', 'F3645', 'F3701', 'F3784', 'F3856', 'F3936', 'F3968', 'F4013', 'F4398', 'F4446', 'F4455', 'F4466', 'F4547', 'F4557', 'F4637', 'F4678', 'F4706', 'F4730', 'F4876', 'F5029', 'F5210', 'F5224', 'F5292', 'F5326', 'F5352', 'F5374', 'F5423', 'F5430', 'F5457', 'F5467', 'F5502', 'F5620', 'F5995', 'F6014', 'F6161', 'F6222', 'F6413', 'F6490', 'F6503', 'F6736', 'F6893', 'F6970', 'F7160', 'F7192', 'F7230', 'F7322', 'F7356', 'F7366', 'F7381', 'F7390', 'F7411', 'F7413', 'F7437', 'F7442', 'F7468', 'F7500', 'F7536', 'F7867', 'F7901', 'F7927', 'F7947', 'F8096', 'F8102', 'F8103', 'F8334', 'F8394', 'F8577', 'F8631', 'F8644', 'F8790', 'F8819', 'F8922', 'F8948', 'F8958', 'F9004', 'F9080', 'F9125', 'F9138', 'F9410', 'F9456', 'F9476', 'F9484', 'F9493', 'F9549', 'F9556', 'F9598', 'F9653', 'F9747', 'F9751', 'F9817', 'F9856', 'F9921', 'F9938', 'F9952', 'F9959'
```

1. Filter data to only PD patient with 'low' or 'high' bdi_category. Also remove any mitochondira/chloroplast DNA

    ```
    qiime feature-table filter-samples \
      --i-table control-table.qza \
      --m-metadata-file pd_metadata_bdi_categorized.txt \
      --p-where "[bdi_category] IN ('low', 'high')" \
      --o-filtered-table pd-low-high-bdi-table.qza

    qiime taxa filter-table \
      --i-table pd-low-high-bdi-table.qza \
      --i-taxonomy taxonomy.qza \
      --p-exclude mitochondria,chloroplast \
      --o-filtered-table table-no-mitochondria-no-chloroplast.qza
    ```
  
2. Visualize data after denoising and filtering

    ```
    qiime feature-table summarize \
      --i-table table-no-mitochondria-no-chloroplast.qza \
      --o-visualization table-no-mitochondria-no-chloroplast.qzv \
      --m-sample-metadata-file pd_metadata_bdi_categorized.txt
    ```

3. Alpha-rarefaction to normalize sampling depth based on BDI_category metadata for downstream alpha/beta diversity analysis
-  Chose a sampling depth 5696. This allows us to retain 296,192 (61.69%) features in 52 (78.79%) samples with:
    - 14 samples for high BDI category
    - 38 samples for low BDI category

    ```
    qiime diversity alpha-rarefaction \
      --i-table table-no-mitochondria-no-chloroplast.qza \
      --i-phylogeny rooted-tree.qza \
      --p-max-depth 5696 \
      --m-metadata-file pd_metadata_bdi_categorized.txt \
      --o-visualization alpha-rarefaction.qzv
    ```
4. Export `table-no-mitochondria-no-chloroplast.qza` as a .txt file
    ```
    qiime tools export \
      --input-path table-no-mitochondria-no-chloroplast.qza \
      --output-path out 

    biom convert -i feature-table.biom --to-tsv -o feature-table.txt
    ```
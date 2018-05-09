#!/bin/bash

### EcoBioTec - QIIME2 tutorial

SECONDS=0

qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path source_info.csv --source-format PairedEndFastqManifestPhred33 --output-path demux-paired-end_OSD14.qza

echo 
echo 'qiime took '$SECONDS' seconds to import with success fastq files'
echo

SECONDS=0

qiime demux summarize --i-data demux-paired-end_OSD14.qza --o-visualization demux-paired-end_OSD14.qzv

echo 
echo 'qiime took '$SECONDS' seconds to summarize fastq files'
echo

SECONDS=0

for fastq in fastq/*; do /Applications/FastQC.app/Contents/MacOS/fastqc $fastq; done;
mkdir fastqc
mv fastq/*.zip fastq/*.html fastqc/

echo 
echo 'FastQC took '$SECONDS' seconds to summarize fastq files'
echo

SECONDS=0

qiime dada2 denoise-paired --i-demultiplexed-seqs demux-paired-end_OSD14.qza --p-trunc-len-f 225 --p-trunc-len-r 225 --p-n-reads-learn 30000 --p-n-threads 4 --o-representative-sequences rep-seqs-dada2_OSD14.qza --o-table table-dada2_OSD14.qza --output-dir dada2-dmx-pe_OSD14

echo
echo 'qiime took '$SECONDS' seconds to denoise the paired-end fastq files'
echo

SECONDS=0

qiime feature-table summarize --i-table table-dada2_OSD14.qza --o-visualization table-dada2_OSD14.qzv --m-sample-metadata-file metadata/osd14_metadata.tsv

qiime feature-table tabulate-seqs --i-data rep-seqs-dada2_OSD14.qza --o-visualization rep-seqs-dada2_OSD14.qzv

echo
echo 'qiime took '$SECONDS' seconds to tabulate ASVs and rep. seqs.'
echo

SECONDS=0

wget -O "gg-13-8-99-515-806-nb-classifier.qza" "https://data.qiime2.org/2018.2/common/gg-13-8-99-515-806-nb-classifier.qza"

echo
echo 'It took '$SECONDS' seconds to download Greengenes database.'
echo

SECONDS=0

qiime feature-classifier classify-sklearn --i-classifier gg-13-8-99-515-806-nb-classifier.qza --i-reads rep-seqs-dada2_OSD14.qza --o-classification taxonomy-rep-seqs-dada2_OSD14.qza

echo
echo 'qiime took '$SECONDS' seconds to classify ASVs to Greengenes database.'
echo

SECONDS=0

qiime metadata tabulate --m-input-file taxonomy-rep-seqs-dada2_OSD14.qza --o-visualization taxonomy-rep-seqs-dada2_OSD14.qzv 
qiime taxa barplot --i-table table-dada2_OSD14.qza --i-taxonomy taxonomy-rep-seqs-dada2_OSD14.qza --m-metadata-file metadata/osd14_metadata.tsv --o-visualization taxa-bar-plots_OSD14.qzv

echo
echo 'qiime took '$SECONDS' seconds to tabulate ASVs and plot them.'
echo

SECONDS=0

qiime alignment mafft --i-sequences rep-seqs-dada2_OSD14.qza --o-alignment mafft-rep-seqs-dada2_OSD14.qza

echo
echo 'mafft took '$SECONDS' seconds to align rep. seqs of ASVs.'
echo

SECONDS=0

qiime alignment mask --i-alignment mafft-rep-seqs-dada2_OSD14.qza --o-masked-alignment masked-msa-rep-seqs-dada2_OSD14.qza
qiime phylogeny fasttree --i-alignment masked-msa-rep-seqs-dada2_OSD14.qza --o-tree unroot-ml-tree-masked_OSD14.qza
qiime phylogeny midpoint-root --i-tree unroot-ml-tree-masked_OSD14.qza --o-rooted-tree root-ml-tree_OSD14.qza

echo
echo 'qiime took '$SECONDS' seconds to mask highly variable regions, build a ML-tree with FastTree, and root it through the midpoint rooting method.'
echo

SECONDS=0

qiime diversity core-metrics-phylogenetic --i-phylogeny root-ml-tree_OSD14.qza --i-table table-dada2_OSD14.qza --p-sampling-depth 4300 --m-metadata-file metadata/osd14_metadata.tsv --output-dir core-metrics-results

echo
echo 'qiime took '$SECONDS' seconds to perform core diversity analyses.'
echo
echo 'qiime concluded with success all the steps!'
echo
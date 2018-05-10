#!/bin/bash

### EcoBioTec - QIIME2 tutorial


echo 
echo 'Starting QIIME2 tutorial!'
echo

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

wget -O "gg-13-8-99-classifier.qza" "https://data.qiime2.org/2018.4/common/gg-13-8-99-nb-classifier.qza"

echo
echo 'It took '$SECONDS' seconds to download Greengenes database.'
echo

md5_checksum=[ md5 gg-13-8-99-classifier.qza | cut -d ' ' -f 4;
if [ md5_checksum==bb72a9e3f1a4c810dd50bceef3508105 ]; then echo; echo 'Greengenes DB is fine!'; echo; else echo; echo 'Greengenes DB is corrupted!'; echo; fi;

SECONDS=0

qiime feature-classifier classify-sklearn --i-classifier gg-13-8-99-classifier.qza --i-reads rep-seqs-dada2_OSD14.qza --o-classification taxonomy-rep-seqs-dada2_OSD14.qza

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

qiime diversity core-metrics-phylogenetic --i-phylogeny root-ml-tree_OSD14.qza --i-table table-dada2_OSD14.qza --p-sampling-depth 20000 --m-metadata-file metadata/osd14_metadata.tsv --output-dir core-metrics-results

echo
echo 'qiime took '$SECONDS' seconds to perform core diversity analyses.'
echo
echo 'QIIME2 tutorial concluded with success all the steps!'
echo

#### EcoBioTec PICRUSt tutorial

echo 
echo 'Starting PICRUSt tutorial!'
echo

mkdir PICRUSt
cd PICRUSt
cp ../demux-paired-end_OSD14.qza .

SECONDS=0

qiime vsearch join-pairs --i-demultiplexed-seqs demux-paired-end_OSD14.qza --p-allowmergestagger --o-joined-sequences dmx-jpe_OSD14.qza

echo
echo 'vsearch took '$SECONDS' seconds to join paired-end reads.'
echo

SECONDS=0

qiime quality-filter q-score-joined --i-demux dmx-jpe_OSD14.qza --o-filtered-sequences dmx-jpe-filter_OSD14.qza --o-filter-stats dmx-jpe-filter-stats.qza

echo
echo 'qiime took '$SECONDS' seconds to filter the paired-end joined reads.'
echo

SECONDS=0

qiime vsearch dereplicate-sequences --i-sequences dmx-jpe-filter_OSD14.qza --o-dereplicated-table drpl-tbl_OSD14.qza --o-dereplicated-sequences drpl-seqs.qza

echo
echo 'vsearch took '$SECONDS' seconds to dereplicate sequences.'
echo

SECONDS=0

wget ftp://greengenes.microbio.me/greengenes_release/gg_13_5/gg_13_5_otus.tar.gz
tar -xvzf gg_13_5_otus.tar.gz # uncompress DB
qiime tools import --type 'FeatureData[Sequence]' --input-path gg_13_5_otus/rep_set/97_otus.fasta --output-path 97_otus-GG.qza # import fasta seqs
qiime tools import --type 'FeatureData[Taxonomy]' --source-format HeaderlessTSVTaxonomyFormat --input-path gg_13_5_otus/taxonomy/97_otu_taxonomy.txt --output-path 97_otu-ref-taxonomy-GG.qza

echo
echo 'It took '$SECONDS' seconds to download Greengenes DB (v.13.5), uncompress it and import it to QIIME2.'
echo

SECONDS=0

qiime vsearch cluster-features-closed-reference --i-table drpl-tbl_OSD14.qza --i-sequences drpl-seqs.qza --i-reference-sequences 97_otus-GG.qza --p-perc-identity 0.97 --o-clustered-table tbl-cr-97_OSD14.qza --o-clustered-sequences rep-seqs-cr-97_OSD14.qza --o-unmatched-sequences unmatched-cr-97_OSD14.qza

echo
echo 'vsearch took '$SECONDS' seconds to cluster OTUs.'
echo
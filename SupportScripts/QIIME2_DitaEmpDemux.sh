"""Script to demultiplex Dita's EMP 16s data"""

conda activate qiime2-2022.2

qiime tools import \
   --type EMPPairedEndSequences \
   --input-path '/home/j/Dropbox/BioinformaticPipeline_Env/MultiplexedFastqs/Dita16sMultiplexed/Run1' \
   --output-path '/home/j/Dropbox/QIIME2/Dita16s_demultiplex/Run1.qza'


qiime demux emp-paired \
  --m-barcodes-file '/home/j/Dropbox/BioinformaticPipeline_Env/MultiplexedFastqs/Dita16sMultiplexed/Barcodes/Run1.tsv' \
  --m-barcodes-column barcode-sequence \
  --i-seqs  '/home/j/Dropbox/QIIME2/Dita16s_demultiplex/Run1.qza' \
  --p-no-golay-error-correction \
  --o-per-sample-sequences '/home/j/Dropbox/QIIME2/Dita16s_demultiplex/Run1_demux-full.qza' \
  --o-error-correction-details '/home/j/Dropbox/QIIME2/Dita16s_demultiplex/Run1_demux-details.qza'

qiime tools export \
    --input-path '/home/j/Dropbox/QIIME2/Dita16s_demultiplex/Run1_demux-full.qza'  \
    --output-path '/home/j/Dropbox/BioinformaticPipeline_Env/FASTQs/Dita16s/Run1'





qiime tools import \
   --type EMPPairedEndSequences \
   --input-path '/home/j/Dropbox/BioinformaticPipeline_Env/MultiplexedFastqs/Dita16sMultiplexed/Run2' \
   --output-path '/home/j/Dropbox/QIIME2/Dita16s_demultiplex/Run2.qza'


qiime demux emp-paired \
  --m-barcodes-file '/home/j/Dropbox/BioinformaticPipeline_Env/MultiplexedFastqs/Dita16sMultiplexed/Barcodes/Run2.tsv' \
  --m-barcodes-column barcode-sequence \
  --i-seqs  '/home/j/Dropbox/QIIME2/Dita16s_demultiplex/Run2.qza' \
  --p-no-golay-error-correction \
  --o-per-sample-sequences '/home/j/Dropbox/QIIME2/Dita16s_demultiplex/Run2_demux-full.qza' \
  --o-error-correction-details '/home/j/Dropbox/QIIME2/Dita16s_demultiplex/Run2_demux-details.qza'

qiime tools export \
    --input-path '/home/j/Dropbox/QIIME2/Dita16s_demultiplex/Run2_demux-full.qza'  \
    --output-path '/home/j/Dropbox/BioinformaticPipeline_Env/FASTQs/Dita16s/Run2'





qiime tools import \
   --type EMPPairedEndSequences \
   --input-path '/home/j/Dropbox/BioinformaticPipeline_Env/MultiplexedFastqs/Dita16sMultiplexed/Run3' \
   --output-path '/home/j/Dropbox/QIIME2/Dita16s_demultiplex/Run3.qza'


qiime demux emp-paired \
  --m-barcodes-file '/home/j/Dropbox/BioinformaticPipeline_Env/MultiplexedFastqs/Dita16sMultiplexed/Barcodes/Run3.tsv' \
  --m-barcodes-column barcode-sequence \
  --i-seqs  '/home/j/Dropbox/QIIME2/Dita16s_demultiplex/Run3.qza' \
  --p-no-golay-error-correction \
  --o-per-sample-sequences '/home/j/Dropbox/QIIME2/Dita16s_demultiplex/Run3_demux-full.qza' \
  --o-error-correction-details '/home/j/Dropbox/QIIME2/Dita16s_demultiplex/Run3_demux-details.qza'

qiime tools export \
    --input-path '/home/j/Dropbox/QIIME2/Dita16s_demultiplex/Run3_demux-full.qza'  \
    --output-path '/home/j/Dropbox/BioinformaticPipeline_Env/FASTQs/Dita16s/Run3'




qiime tools import \
   --type EMPPairedEndSequences \
   --input-path '/home/j/Dropbox/BioinformaticPipeline_Env/MultiplexedFastqs/Dita16sMultiplexed/Run4' \
   --output-path '/home/j/Dropbox/QIIME2/Dita16s_demultiplex/Run4.qza'


qiime demux emp-paired \
  --m-barcodes-file '/home/j/Dropbox/BioinformaticPipeline_Env/MultiplexedFastqs/Dita16sMultiplexed/Barcodes/Run4.tsv' \
  --m-barcodes-column barcode-sequence \
  --i-seqs  '/home/j/Dropbox/QIIME2/Dita16s_demultiplex/Run4.qza' \
  --p-no-golay-error-correction \
  --o-per-sample-sequences '/home/j/Dropbox/QIIME2/Dita16s_demultiplex/Run4_demux-full.qza' \
  --o-error-correction-details '/home/j/Dropbox/QIIME2/Dita16s_demultiplex/Run4_demux-details.qza'

qiime tools export \
    --input-path '/home/j/Dropbox/QIIME2/Dita16s_demultiplex/Run4_demux-full.qza'  \
    --output-path '/home/j/Dropbox/BioinformaticPipeline_Env/FASTQs/Dita16s/Run4'



qiime tools import \
   --type EMPPairedEndSequences \
   --input-path '/home/j/Dropbox/BioinformaticPipeline_Env/MultiplexedFastqs/Dita16sMultiplexed/Run5' \
   --output-path '/home/j/Dropbox/QIIME2/Dita16s_demultiplex/Run5.qza'


qiime demux emp-paired \
  --m-barcodes-file '/home/j/Dropbox/BioinformaticPipeline_Env/MultiplexedFastqs/Dita16sMultiplexed/Barcodes/Run5.tsv' \
  --m-barcodes-column barcode-sequence \
  --i-seqs  '/home/j/Dropbox/QIIME2/Dita16s_demultiplex/Run5.qza' \
  --p-no-golay-error-correction \
  --o-per-sample-sequences '/home/j/Dropbox/QIIME2/Dita16s_demultiplex/Run5_demux-full.qza' \
  --o-error-correction-details '/home/j/Dropbox/QIIME2/Dita16s_demultiplex/Run5_demux-details.qza'

qiime tools export \
    --input-path '/home/j/Dropbox/QIIME2/Dita16s_demultiplex/Run5_demux-full.qza'  \
    --output-path '/home/j/Dropbox/BioinformaticPipeline_Env/FASTQs/Dita16s/Run5'




qiime tools import \
   --type EMPPairedEndSequences \
   --input-path '/home/j/Dropbox/BioinformaticPipeline_Env/MultiplexedFastqs/Dita16sMultiplexed/Run6' \
   --output-path '/home/j/Dropbox/QIIME2/Dita16s_demultiplex/Run6.qza'


qiime demux emp-paired \
  --m-barcodes-file '/home/j/Dropbox/BioinformaticPipeline_Env/MultiplexedFastqs/Dita16sMultiplexed/Barcodes/Run6.tsv' \
  --m-barcodes-column barcode-sequence \
  --i-seqs  '/home/j/Dropbox/QIIME2/Dita16s_demultiplex/Run6.qza' \
  --p-no-golay-error-correction \
  --o-per-sample-sequences '/home/j/Dropbox/QIIME2/Dita16s_demultiplex/Run6_demux-full.qza' \
  --o-error-correction-details '/home/j/Dropbox/QIIME2/Dita16s_demultiplex/Run6_demux-details.qza'

qiime tools export \
    --input-path '/home/j/Dropbox/QIIME2/Dita16s_demultiplex/Run6_demux-full.qza'  \
    --output-path '/home/j/Dropbox/BioinformaticPipeline_Env/FASTQs/Dita16s/Run6'





qiime tools import \
   --type EMPPairedEndSequences \
   --input-path '/home/j/Dropbox/BioinformaticPipeline_Env/MultiplexedFastqs/Dita16sMultiplexed/Run7' \
   --output-path '/home/j/Dropbox/QIIME2/Dita16s_demultiplex/Run7.qza'


qiime demux emp-paired \
  --m-barcodes-file '/home/j/Dropbox/BioinformaticPipeline_Env/MultiplexedFastqs/Dita16sMultiplexed/Barcodes/Run7.tsv' \
  --m-barcodes-column barcode-sequence \
  --i-seqs  '/home/j/Dropbox/QIIME2/Dita16s_demultiplex/Run7.qza' \
  --p-no-golay-error-correction \
  --o-per-sample-sequences '/home/j/Dropbox/QIIME2/Dita16s_demultiplex/Run7_demux-full.qza' \
  --o-error-correction-details '/home/j/Dropbox/QIIME2/Dita16s_demultiplex/Run7_demux-details.qza'

qiime tools export \
    --input-path '/home/j/Dropbox/QIIME2/Dita16s_demultiplex/Run7_demux-full.qza'  \
    --output-path '/home/j/Dropbox/BioinformaticPipeline_Env/FASTQs/Dita16s/Run7'



qiime tools import \
   --type EMPPairedEndSequences \
   --input-path '/home/j/Dropbox/BioinformaticPipeline_Env/MultiplexedFastqs/Dita16sMultiplexed/Run8' \
   --output-path '/home/j/Dropbox/QIIME2/Dita16s_demultiplex/Run8.qza'


qiime demux emp-paired \
  --m-barcodes-file '/home/j/Dropbox/BioinformaticPipeline_Env/MultiplexedFastqs/Dita16sMultiplexed/Barcodes/Run8.tsv' \
  --m-barcodes-column barcode-sequence \
  --i-seqs  '/home/j/Dropbox/QIIME2/Dita16s_demultiplex/Run8.qza' \
  --p-no-golay-error-correction \
  --o-per-sample-sequences '/home/j/Dropbox/QIIME2/Dita16s_demultiplex/Run8_demux-full.qza' \
  --o-error-correction-details '/home/j/Dropbox/QIIME2/Dita16s_demultiplex/Run8_demux-details.qza'

qiime tools export \
    --input-path '/home/j/Dropbox/QIIME2/Dita16s_demultiplex/Run8_demux-full.qza'  \
    --output-path '/home/j/Dropbox/BioinformaticPipeline_Env/FASTQs/Dita16s/Run8'



"""create visualisations of demuxed runs"""

qiime demux summarize \
  --i-data '/home/j/Dropbox/QIIME2/Dita16s_demultiplex/Run1_demux-full.qza'   \
  --o-visualization '/home/j/Dropbox/QIIME2/Dita16s_demultiplex/Run1_demux-full.qzv'


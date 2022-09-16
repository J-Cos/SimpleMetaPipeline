"""Script to demultiplex CrossPacific EMP 16s data"""

conda activate qiime2-2022.2

qiime tools import \
   --type EMPPairedEndSequences \
   --input-path '/home/j/Dropbox/BioinformaticPipeline_Env/FASTQs/CrossPacific/150718_M03292_0022_000000000-AGKWR' \
   --output-path '/home/j/Dropbox/QIIME2/CrossPacific16s_demultiplex/150718_M03292_0022_000000000-AGKWR.qza'


qiime demux emp-paired \
  --m-barcodes-file '/home/j/Dropbox/BioinformaticPipeline_Env/FASTQs/CrossPacific/Barcodes/Run1.tsv' \
  --m-barcodes-column barcode-sequence \
  --i-seqs  '/home/j/Dropbox/QIIME2/CrossPacific16s_demultiplex/150718_M03292_0022_000000000-AGKWR.qza' \
  --p-no-golay-error-correction \
  --o-per-sample-sequences /home/j/Dropbox/QIIME2/CrossPacific16s_demultiplex/150718_M03292_0022_000000000-AGKWR_demux-full.qza \
  --o-error-correction-details /home/j/Dropbox/QIIME2/CrossPacific16s_demultiplex/150718_M03292_0022_000000000-AGKWR_demux-details.qza

qiime tools export \
    --input-path /home/j/Dropbox/QIIME2/CrossPacific16s_demultiplex/150718_M03292_0022_000000000-AGKWR_demux-full.qza \
    --output-path /home/j/Dropbox/QIIME2/CrossPacific16s_demultiplex/DemultiplexedFastqs






qiime tools import \
   --type EMPPairedEndSequences \
   --input-path '/home/j/Dropbox/BioinformaticPipeline_Env/FASTQs/CrossPacific/151110_M03292_0039_000000000-AHLEH' \
   --output-path '/home/j/Dropbox/QIIME2/CrossPacific16s_demultiplex/151110_M03292_0039_000000000-AHLEH.qza'


qiime demux emp-paired \
  --m-barcodes-file '/home/j/Dropbox/BioinformaticPipeline_Env/FASTQs/CrossPacific/Barcodes/Run2.tsv' \
  --m-barcodes-column barcode.sequence \
  --i-seqs  '/home/j/Dropbox/QIIME2/CrossPacific16s_demultiplex/151110_M03292_0039_000000000-AHLEH.qza' \
  --p-no-golay-error-correction \
  --o-per-sample-sequences /home/j/Dropbox/QIIME2/CrossPacific16s_demultiplex/151110_M03292_0039_000000000-AHLEH_demux-full.qza \
  --o-error-correction-details /home/j/Dropbox/QIIME2/CrossPacific16s_demultiplex/151110_M03292_0039_000000000-AHLEH_demux-details.qza

qiime tools export \
    --input-path /home/j/Dropbox/QIIME2/CrossPacific16s_demultiplex/151110_M03292_0039_000000000-AHLEH_demux-full.qza \
    --output-path /home/j/Dropbox/QIIME2/CrossPacific16s_demultiplex/DemultiplexedFastqs






qiime tools import \
   --type EMPPairedEndSequences \
   --input-path '/home/j/Dropbox/BioinformaticPipeline_Env/FASTQs/CrossPacific/151116_M03292_0040_000000000-AKEDL' \
   --output-path '/home/j/Dropbox/QIIME2/CrossPacific16s_demultiplex/151116_M03292_0040_000000000-AKEDL.qza'


qiime demux emp-paired \
  --m-barcodes-file '/home/j/Dropbox/BioinformaticPipeline_Env/FASTQs/CrossPacific/Barcodes/Run3.tsv' \
  --m-barcodes-column barcode.sequence \
  --i-seqs  '/home/j/Dropbox/QIIME2/CrossPacific16s_demultiplex/151116_M03292_0040_000000000-AKEDL.qza' \
  --p-no-golay-error-correction \
  --o-per-sample-sequences /home/j/Dropbox/QIIME2/CrossPacific16s_demultiplex/151116_M03292_0040_000000000-AKEDL_demux-full.qza \
  --o-error-correction-details /home/j/Dropbox/QIIME2/CrossPacific16s_demultiplex/151116_M03292_0040_000000000-AKEDL_demux-details.qza

qiime tools export \
    --input-path /home/j/Dropbox/QIIME2/CrossPacific16s_demultiplex/151116_M03292_0040_000000000-AKEDL_demux-full.qza \
    --output-path /home/j/Dropbox/QIIME2/CrossPacific16s_demultiplex/DemultiplexedFastqs






qiime tools import \
   --type EMPPairedEndSequences \
   --input-path '/home/j/Dropbox/BioinformaticPipeline_Env/FASTQs/CrossPacific/151210_M03292_0045_000000000-AKE7C' \
   --output-path '/home/j/Dropbox/QIIME2/CrossPacific16s_demultiplex/151210_M03292_0045_000000000-AKE7C.qza'


qiime demux emp-paired \
  --m-barcodes-file '/home/j/Dropbox/BioinformaticPipeline_Env/FASTQs/CrossPacific/Barcodes/Run4.tsv' \
  --m-barcodes-column barcode.sequence \
  --i-seqs  '/home/j/Dropbox/QIIME2/CrossPacific16s_demultiplex/151210_M03292_0045_000000000-AKE7C.qza' \
  --p-no-golay-error-correction \
  --o-per-sample-sequences /home/j/Dropbox/QIIME2/CrossPacific16s_demultiplex/151210_M03292_0045_000000000-AKE7C_demux-full.qza \
  --o-error-correction-details /home/j/Dropbox/QIIME2/CrossPacific16s_demultiplex/151210_M03292_0045_000000000-AKE7C_demux-details.qza

qiime tools export \
    --input-path /home/j/Dropbox/QIIME2/CrossPacific16s_demultiplex/151210_M03292_0045_000000000-AKE7C_demux-full.qza \
    --output-path /home/j/Dropbox/QIIME2/CrossPacific16s_demultiplex/DemultiplexedFastqs






qiime tools import \
   --type EMPPairedEndSequences \
   --input-path '/home/j/Dropbox/BioinformaticPipeline_Env/FASTQs/CrossPacific/151222_M03292_0048_000000000-AL8TW' \
   --output-path '/home/j/Dropbox/QIIME2/CrossPacific16s_demultiplex/151222_M03292_0048_000000000-AL8TW.qza'


qiime demux emp-paired \
  --m-barcodes-file '/home/j/Dropbox/BioinformaticPipeline_Env/FASTQs/CrossPacific/Barcodes/Run5.tsv' \
  --m-barcodes-column barcode.sequence \
  --i-seqs  '/home/j/Dropbox/QIIME2/CrossPacific16s_demultiplex/151222_M03292_0048_000000000-AL8TW.qza' \
  --p-no-golay-error-correction \
  --o-per-sample-sequences /home/j/Dropbox/QIIME2/CrossPacific16s_demultiplex/151222_M03292_0048_000000000-AL8TW_demux-full.qza \
  --o-error-correction-details /home/j/Dropbox/QIIME2/CrossPacific16s_demultiplex/151222_M03292_0048_000000000-AL8TW_demux-details.qza

qiime tools export \
    --input-path /home/j/Dropbox/QIIME2/CrossPacific16s_demultiplex/151222_M03292_0048_000000000-AL8TW_demux-full.qza \
    --output-path /home/j/Dropbox/QIIME2/CrossPacific16s_demultiplex/DemultiplexedFastqs






qiime tools import \
   --type EMPPairedEndSequences \
   --input-path '/home/j/Dropbox/BioinformaticPipeline_Env/FASTQs/CrossPacific/160322_M03292_0061_000000000-AL86W' \
   --output-path '/home/j/Dropbox/QIIME2/CrossPacific16s_demultiplex/160322_M03292_0061_000000000-AL86W.qza'


qiime demux emp-paired \
  --m-barcodes-file '/home/j/Dropbox/BioinformaticPipeline_Env/FASTQs/CrossPacific/Barcodes/Run6.tsv' \
  --m-barcodes-column barcode.sequence \
  --i-seqs  '/home/j/Dropbox/QIIME2/CrossPacific16s_demultiplex/160322_M03292_0061_000000000-AL86W.qza' \
  --p-no-golay-error-correction \
  --o-per-sample-sequences /home/j/Dropbox/QIIME2/CrossPacific16s_demultiplex/160322_M03292_0061_000000000-AL86W_demux-full.qza \
  --o-error-correction-details /home/j/Dropbox/QIIME2/CrossPacific16s_demultiplex/160322_M03292_0061_000000000-AL86W_demux-details.qza

qiime tools export \
    --input-path /home/j/Dropbox/QIIME2/CrossPacific16s_demultiplex/160322_M03292_0061_000000000-AL86W_demux-full.qza \
    --output-path /home/j/Dropbox/QIIME2/CrossPacific16s_demultiplex/DemultiplexedFastqs






qiime tools import \
   --type EMPPairedEndSequences \
   --input-path '/home/j/Dropbox/BioinformaticPipeline_Env/FASTQs/CrossPacific/160328_M03292_0063_000000000-ABD00' \
   --output-path '/home/j/Dropbox/QIIME2/CrossPacific16s_demultiplex/160328_M03292_0063_000000000-ABD00.qza'


qiime demux emp-paired \
  --m-barcodes-file '/home/j/Dropbox/BioinformaticPipeline_Env/FASTQs/CrossPacific/Barcodes/Run7.tsv' \
  --m-barcodes-column barcode.sequence \
  --i-seqs  '/home/j/Dropbox/QIIME2/CrossPacific16s_demultiplex/160328_M03292_0063_000000000-ABD00.qza' \
  --p-no-golay-error-correction \
  --o-per-sample-sequences /home/j/Dropbox/QIIME2/CrossPacific16s_demultiplex/160328_M03292_0063_000000000-ABD00_demux-full.qza \
  --o-error-correction-details /home/j/Dropbox/QIIME2/CrossPacific16s_demultiplex/160328_M03292_0063_000000000-ABD00_demux-details.qza

qiime tools export \
    --input-path /home/j/Dropbox/QIIME2/CrossPacific16s_demultiplex/160328_M03292_0063_000000000-ABD00_demux-full.qza \
    --output-path /home/j/Dropbox/QIIME2/CrossPacific16s_demultiplex/DemultiplexedFastqs/ABD00



qiime tools import \
   --type EMPPairedEndSequences \
   --input-path '/home/j/Dropbox/BioinformaticPipeline_Env/FASTQs/CrossPacific/160422_M03292_0005_000000000-ANL7K' \
   --output-path '/home/j/Dropbox/QIIME2/CrossPacific16s_demultiplex/160422_M03292_0005_000000000-ANL7K.qza'


qiime demux emp-paired \
  --m-barcodes-file '/home/j/Dropbox/BioinformaticPipeline_Env/FASTQs/CrossPacific/Barcodes/Run7.tsv' \
  --m-barcodes-column barcode.sequence \
  --i-seqs  '/home/j/Dropbox/QIIME2/CrossPacific16s_demultiplex/160422_M03292_0005_000000000-ANL7K.qza' \
  --p-no-golay-error-correction \
  --o-per-sample-sequences /home/j/Dropbox/QIIME2/CrossPacific16s_demultiplex/160328_M03292_0063_000000000-ANL7K_demux-full.qza \
  --o-error-correction-details /home/j/Dropbox/QIIME2/CrossPacific16s_demultiplex/160328_M03292_0063_000000000-ANL7K_demux-details.qza

qiime tools export \
    --input-path /home/j/Dropbox/QIIME2/CrossPacific16s_demultiplex/160328_M03292_0063_000000000-ANL7K_demux-full.qza \
    --output-path /home/j/Dropbox/QIIME2/CrossPacific16s_demultiplex/DemultiplexedFastqs/Duplicate







"""create visualisations of demuxed runs"""

qiime demux summarize \
  --i-data /home/j/Dropbox/QIIME2/CrossPacific16s_demultiplex/160328_M03292_0063_000000000-ANL7K_demux-full.qza  \
  --o-visualization /home/j/Dropbox/QIIME2/CrossPacific16s_demultiplex/160328_M03292_0063_000000000-ANL7K_demux-full.qzv


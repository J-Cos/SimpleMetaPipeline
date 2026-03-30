    path<-'/rds/general/user/jcw120/home/BioinformaticPipeline_Env/RawFastqs/RansomeAll_COI/RunsSplitByAdapter'
    #path<-'/rds/general/user/jcw120/home/BioinformaticPipeline_Env/RawFastqs/DitaCOIMultiplexed/RunsSplitByAdapter'
    #path<-'/home/j/Dropbox/BioinformaticPipeline_Env/MultiplexedFastqs/EmmaCOIMultiplexed'

    #R1<-sort(list.files(file.path(path, "Cahyani2017_COI_Run1"), pattern=".fastq.gz", full.names = TRUE))
    #R2<-sort(list.files(file.path(path, "Cahyani2017_COI_Run2"), pattern=".fastq.gz", full.names = TRUE))
    #R3<-sort(list.files(file.path(path, "Cahyani2017_COI_Run3"), pattern=".fastq.gz", full.names = TRUE))
    #R4<-sort(list.files(file.path(path, "Cahyani2017_COI_Run4"), pattern=".fastq.gz", full.names = TRUE))

    #R1<-sort(list.files("Dunn_221118_Data_PLAS_09_12", pattern=".fastq.gz", full.names = TRUE))
    #R2<-sort(list.files("Dunn_221125_Data_PLAS_01_04", pattern=".fastq.gz", full.names = TRUE))
    #R3<-sort(list.files("Dunn_221202_Data_PLAS_05_08", pattern=".fastq.gz", full.names = TRUE))


    R1<-sort(list.files(file.path(path, "151203_M03292_0044_000000000-AKYV3"), pattern=".fastq.gz", full.names = TRUE))
    R2<-sort(list.files(file.path(path, "151216_M03292_0046_000000000-AJMWE"), pattern=".fastq.gz", full.names = TRUE))
    R3<-sort(list.files(file.path(path, "160508_M03292_0010_000000000-ANT9R"), pattern=".fastq.gz", full.names = TRUE))
    R4<-sort(list.files(file.path(path, "CO1_run6"), pattern=".fastq.gz", full.names = TRUE))
    R5<-sort(list.files(file.path(path, "EmmaCO1_23744943_R1_PLUS_ARMS_R2"), pattern=".fastq.gz", full.names = TRUE))
    R6<-sort(list.files(file.path(path, "ShannonCO1-24072105_R2"), pattern=".fastq.gz", full.names = TRUE))

#pdf(file.path('/rds/general/user/jcw120/home/RansomeAllCOI_quality.pdf'))
pdf(file.path('/rds/general/user/jcw120/home/CahyaniAllCOI_quality.pdf'))
pdf("DunnQuality2.pdf")
    dada2::plotQualityProfile(R1[1:16])
    dada2::plotQualityProfile(R2[1:16])
    dada2::plotQualityProfile(R3[1:16])
    #dada2::plotQualityProfile(R4[1:4])
    #dada2::plotQualityProfile(R5[1:4])
    #dada2::plotQualityProfile(R6[1:4])
    #dada2::plotQualityProfile(c(R7))
    #dada2::plotQualityProfile(c(R8))
    #dada2::plotQualityProfile(c(R9))
    #dada2::plotQualityProfile(c(R10))
    #dada2::plotQualityProfile(c(R11))
dev.off()
#quantify contents of fastq

#dependencies and parameters
library(ShortRead)
path<-'/rds/general/user/jcw120/home/BioinformaticPipeline_Env/RawFastqs/RansomeAll_COI/RunsSplitByAdapter'

#functions
CountFastqRecordsPerRun<-function(Run){
    R<-sort(list.files(file.path(path, Run), pattern=".fastq.gz", full.names = TRUE))
    table<-ShortRead::countFastq(c(R))
    write.csv(table, file=file.path(path, "../..", paste0(Run, "_FastqRecordsCount.csv")))
}

#code
for (Run in list.files(path)){
    CountFastqRecordsPerRun(Run)
    print(paste0(Run, " complete"))
}
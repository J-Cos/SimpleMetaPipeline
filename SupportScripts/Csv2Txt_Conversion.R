#csv to tsv conversion

library(tidyverse)

path<-"/home/j/Dropbox/BioinformaticPipeline_Env/MultiplexedFastqs/EmmaCOIMultiplexed/Barcodes/160508_M03292_0010_000000000-ANT9R"

barcodecsvs<-list.files(file.path(path), pattern=".csv")
barcodetsvs<-barcodecsvs %>% str_replace(".csv", replacement=".txt")

for (i in 1:length(barcodecsvs)) {
    read.csv(file.path(path, barcodecsvs[i]), header=FALSE) %>% 
    write.table(file.path(path,barcodetsvs[i]), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
}
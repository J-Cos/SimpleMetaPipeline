#csv to tsv conversion

library(tidyverse)

path<-"/home/j/Dropbox/BioinformaticPipeline_Env/MultiplexedFastqs/Dita16sMultiplexed/Barcodes"

barcodecsvs<-list.files(file.path(path), pattern=".csv")
barcodetsvs<-barcodecsvs %>% str_replace(".csv", replacement=".tsv")

for (i in 1:length(barcodecsvs)) {
    readr::read_csv(file.path(path, barcodecsvs[i])) %>% 
    '['(1:2) %>%
    rename(c("sampleID"="sample")) %>%
    write_tsv(file.path(path,barcodetsvs[i]))
}
        loadFastq<-function(FASTQ_folder, pattern){
            output<-sort(list.files(file.path(path, "FASTQs", FASTQ_folder), pattern=pattern, full.names = TRUE))
            return(output)
        }